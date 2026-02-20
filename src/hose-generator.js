/**
 * HOSE Code Generator
 *
 * Generates text-format HOSE codes (Bremser notation) from an openchemlib
 * Molecule object. Compatible with nmrshiftdb2 database format.
 *
 * Implements the CDK ExtendedHOSECodeGenerator two-pass algorithm:
 * Pass 1: BFS tree building - only skip parent atom, no visited tracking
 * Pass 2: Score calculation + code generation - detect ring closures via visited tracking
 *
 * Reference: W. Bremser, "HOSE - A Novel Substructure Code",
 *            Analytica Chimica Acta, 1978, 103, 355-365.
 */

// Bremser element substitutions
const BREMSER = { Si: 'Q', Cl: 'X', Br: 'Y' };

// Element priority for canonical ordering (higher = higher priority)
const ELEMENT_RANK = {
  C: 9000, O: 8900, N: 8800, S: 8700, P: 8600,
  Si: 8500, B: 8400, F: 8300, Cl: 8200, Br: 8100,
  I: 7900,
};

// Ring closure rank (from CDK: RING_RANK = 1100)
const RING_RANK = 1100;

// Bond type ranking (index: 0=unused, 1=single, 2=double, 3=triple, 4=aromatic)
const BOND_RANKINGS = [0, 0, 200000, 300000, 100000];

// Bond type symbols (index: 0=unused, 1=single, 2=double, 3=triple, 4=aromatic)
const BOND_SYMBOLS = ['<', '', '=', '%', '*'];

// Sphere delimiters
const SPHERE_DELIMITERS = ['(', '/', '/', ')', '/', '/', '/', '/', '/', '/'];

/**
 * Generate a text-format HOSE code for a given atom in a molecule.
 */
export function generateHoseCode(mol, atomIndex, options = {}) {
  const { maxSpheres = 4 } = options;

  mol.ensureHelperArrays(mol.constructor.cHelperSymmetrySimple || 63);

  // PASS 1: Build BFS tree without visited tracking (CDK approach)
  // Only skip the parent atom for each node
  const spheres = buildBfsTree(mol, atomIndex, maxSpheres);

  // PASS 2: Score, detect ring closures, sort, and generate code
  return createCode(mol, spheres, maxSpheres, atomIndex);
}

/**
 * Pass 1: Build BFS tree. Only skip parent atom, no visited tracking.
 * All non-parent neighbors are added as children, even if they appeared in earlier spheres.
 * Implicit H atoms are included as leaf nodes.
 * Nodes are sorted by canonical label (symmetry rank) at each level.
 */
function buildBfsTree(mol, centerIdx, maxSpheres) {
  const spheres = [];

  // Sphere 0: all neighbors of center atom
  const sphere0 = [];
  const connCount = mol.getConnAtoms(centerIdx);
  for (let i = 0; i < connCount; i++) {
    const neighborIdx = mol.getConnAtom(centerIdx, i);
    const element = mol.getAtomLabel(neighborIdx);
    const bondIdx = mol.getConnBond(centerIdx, i);
    const bondType = mol.isAromaticBond(bondIdx) ? 4 : mol.getBondOrder(bondIdx);
    const degree = getTotalBondCount(mol, neighborIdx);
    sphere0.push(makeNode(neighborIdx, element, bondType, null, centerIdx, degree));
  }
  // Add implicit H for center
  const centerImplH = mol.getImplicitHydrogens(centerIdx);
  for (let h = 0; h < centerImplH; h++) {
    sphere0.push(makeNode(-1, 'H', 1, null, centerIdx, 1));
  }
  // Sort by canonical label (CDK TreeNodeComparator sorts by canonical label ascending)
  sortByCanonicalLabel(sphere0, mol);
  spheres.push(sphere0);

  // Subsequent spheres
  for (let s = 0; s < maxSpheres - 1; s++) {
    const prevSphere = spheres[s];
    const nextSphereNodes = [];

    for (const node of prevSphere) {
      // Skip symbols that are stoppers/separators
      if (node.element === ',' || node.element === '&' || node.element === '#') continue;
      // Skip H atoms (not expanded)
      if (node.element === 'H') continue;

      const atomIdx = node.atomIdx;
      if (atomIdx < 0) continue;

      const connCount = mol.getConnAtoms(atomIdx);
      const parentAtomIdx = node.sourceAtomIdx;

      const implH = mol.getImplicitHydrogens(atomIdx);

      if (connCount === 1 && implH === 0) {
        // Truly terminal atom: only parent connection, no implicit H. Add comma separator.
        nextSphereNodes.push(makeNode(-1, ',', -1, node, atomIdx, 0));
      } else {
        for (let i = 0; i < connCount; i++) {
          const neighborIdx = mol.getConnAtom(atomIdx, i);
          if (neighborIdx === parentAtomIdx) continue; // Only skip parent

          const element = mol.getAtomLabel(neighborIdx);
          const bondIdx = mol.getConnBond(atomIdx, i);
          const bondType = mol.isAromaticBond(bondIdx) ? 4 : mol.getBondOrder(bondIdx);
          const degree = getTotalBondCount(mol, neighborIdx);
          nextSphereNodes.push(makeNode(neighborIdx, element, bondType, node, atomIdx, degree));
        }
        // Add implicit H children
        for (let h = 0; h < implH; h++) {
          nextSphereNodes.push(makeNode(-1, 'H', 1, node, atomIdx, 1));
        }
      }
    }

    // Sort by canonical label
    sortByCanonicalLabel(nextSphereNodes, mol);
    spheres.push(nextSphereNodes);
  }

  return spheres;
}

/**
 * Pass 2: Score, detect ring closures, sort, generate HOSE code string.
 * Follows CDK createCode() sequence:
 * 1. Degree accumulation (bottom-up)
 * 2. Calculate scores + sort (element rank or RING_RANK based on visited)
 * 3. Add ranking to score + sort
 * 4. Build stringscores + sort
 * 5. Backward propagation + sort
 * 6. Second forward pass rebuild stringscores + sort
 * 7. Generate code string with visited-based ring closure detection
 */
function createCode(mol, spheres, maxSpheres, centerIdx) {
  // Step 1: Degree accumulation (bottom-up from outer to inner)
  for (let f = 0; f < maxSpheres - 1; f++) {
    const sphere = spheres[maxSpheres - 1 - f];
    for (const node of sphere) {
      if (node.source !== null) {
        node.source.ranking += node.degree;
      }
    }
  }

  // Step 2: Calculate node scores and sort (element rank + bond rank, with visited tracking)
  // Visited tracking across spheres: atoms scored in earlier spheres are "visited"
  const visitedAtoms = new Set();
  visitedAtoms.add(centerIdx); // Center is visited

  for (let f = 0; f < maxSpheres; f++) {
    const sphere = spheres[f];
    calculateNodeScores(sphere, visitedAtoms);
    sortNodesByStringscore(sphere);
  }

  // Step 3: Add ranking to score and re-sort
  for (let f = 0; f < maxSpheres; f++) {
    const sphere = spheres[f];
    for (const node of sphere) {
      node.score += node.ranking;
    }
    sortNodesByStringscore(sphere);
  }

  // Step 4: Build stringscores (forward pass) and sort
  for (let f = 0; f < maxSpheres; f++) {
    const sphere = spheres[f];
    for (const node of sphere) {
      const localScore = padScore(node.score);
      node.stringscore = (node.source ? node.source.stringscore : '') + localScore;
    }
    sortNodesByStringscore(sphere);
  }

  // Step 5: Backward propagation - parent takes child's stringscore, re-sort
  for (let f = maxSpheres - 1; f > 0; f--) {
    const sphere = spheres[f];
    const parentSphere = spheres[f - 1];
    for (const node of sphere) {
      for (const parentNode of parentSphere) {
        if (parentNode === node.source) {
          parentNode.stringscore = node.stringscore;
        }
      }
    }
    sortNodesByStringscore(parentSphere);
  }

  // Step 6: Second forward pass - rebuild stringscores and sort
  for (let f = 0; f < maxSpheres; f++) {
    const sphere = spheres[f];
    for (const node of sphere) {
      const localScore = padScore(node.score);
      node.stringscore = (node.source ? node.source.stringscore : '') + localScore;
    }
    sortNodesByStringscore(sphere);
  }

  // Step 7: Generate code string with visited-based ring closure detection
  return generateCodeString(mol, spheres, maxSpheres, centerIdx);
}

/**
 * Calculate node scores: element rank (or RING_RANK for already-visited atoms) + bond rank.
 * Marks atoms as visited AFTER scoring all nodes in the sphere (batch-wise).
 */
function calculateNodeScores(sphere, visitedAtoms) {
  const toMark = [];

  for (const node of sphere) {
    // Ring closure detection: if atom already visited, use RING_RANK
    if (node.atomIdx >= 0 && visitedAtoms.has(node.atomIdx)) {
      node.score += RING_RANK;
      node.isRingClosure = true;
    } else {
      node.score += getElementRank(node.element);
    }

    // Bond ranking
    const bt = node.bondType;
    if (bt >= 0 && bt <= 4) {
      node.score += BOND_RANKINGS[bt];
    } else if (bt === -1) {
      node.score += 50000; // Comma/separator
    }

    if (node.atomIdx >= 0) {
      toMark.push(node.atomIdx);
    }
  }

  // Mark all atoms as visited AFTER scoring (batch-wise, CDK behavior)
  for (const idx of toMark) {
    visitedAtoms.add(idx);
  }
}

/**
 * Generate the HOSE code string from sorted spheres.
 * Ring closures are detected by tracking visited atoms during code output.
 */
function generateCodeString(mol, spheres, maxSpheres, centerIdx) {
  let code = '';
  const visited = new Set();
  visited.add(centerIdx);

  // Sphere 0: prefix before delimiters
  const s0 = spheres[0];
  let branch = s0.length > 0 ? s0[0].source : null; // source of first node
  // For sphere 0, track the first source atom for branch detection
  let branchAtom = s0.length > 0 ? centerIdx : -1;

  for (const node of s0) {
    const bondSym = getBondSymbol(node.bondType);
    if (node.atomIdx >= 0 && visited.has(node.atomIdx)) {
      code += bondSym + '&';
      node.stopper = true;
    } else if (node.atomIdx >= 0) {
      code += bondSym + bremserElement(node.element);
    } else if (node.element === 'H') {
      code += bondSym + 'H';
    } else if (node.element === ',') {
      // Comma separator node — skip in sphere 0
    }
    if (node.atomIdx >= 0) visited.add(node.atomIdx);
    if (node.source && node.source.stopper) node.stopper = true;
  }

  // Subsequent spheres
  for (let f = 0; f < maxSpheres - 1; f++) {
    const sphere = spheres[f + 1];
    const delimIdx = f;
    code += delimIdx < SPHERE_DELIMITERS.length ? SPHERE_DELIMITERS[delimIdx] : '/';

    if (sphere.length === 0) continue;

    // Branch tracking: CDK initializes branch to first node's source atom
    let currentBranch = sphere[0].source ? sphere[0].source.atomIdx : -1;

    for (let i = 0; i < sphere.length; i++) {
      const node = sphere[i];
      const sourceAtom = node.source ? node.source.atomIdx : -1;

      // Comma when source parent changes AND source is not a stopper
      if (node.source && !node.source.stopper && sourceAtom !== currentBranch) {
        currentBranch = sourceAtom;
        code += ',';
        // Check for chirality (@) — we don't implement this
      }

      // Only output code if source is not a stopper
      if (node.source && !node.source.stopper) {
        const bondSym = getBondSymbol(node.bondType);

        if (node.atomIdx >= 0 && visited.has(node.atomIdx)) {
          // Ring closure
          code += bondSym + '&';
          node.stopper = true;
        } else if (node.atomIdx >= 0) {
          code += bondSym + bremserElement(node.element);
        } else if (node.element === 'H') {
          code += bondSym + 'H';
        } else if (node.element === ',') {
          // Comma separator (from terminal atom) — acts as empty branch
          // The symbol is "," but getElementSymbol returns "" for it
          // This just occupies a slot, no output
        }
      }

      if (node.atomIdx >= 0) visited.add(node.atomIdx);
      if (node.source && node.source.stopper) node.stopper = true;
    }
  }

  // Fill up sphere delimiters (CDK fillUpSphereDelimiters)
  // Add remaining delimiters up to maxSpheres
  // The last sphere delimiter was already added in the loop above,
  // but we need the closing ones
  const lastDelimAdded = maxSpheres - 2; // last delimIdx added in loop
  for (let f = lastDelimAdded + 1; f < Math.min(maxSpheres, SPHERE_DELIMITERS.length); f++) {
    code += SPHERE_DELIMITERS[f];
  }

  return code;
}

// Helper functions

function makeNode(atomIdx, element, bondType, source, sourceAtomIdx, degree) {
  return {
    atomIdx,
    element,
    bondType,
    source,           // parent TreeNode
    sourceAtomIdx,    // parent atom index (for parent-skip during BFS)
    degree,           // total molecular bond count
    score: 0,
    ranking: 0,
    sortOrder: 0,
    stringscore: '',
    stopper: false,
    isRingClosure: false,
    children: [],
  };
}

function getTotalBondCount(mol, atomIdx) {
  if (atomIdx < 0) return 1;
  return mol.getConnAtoms(atomIdx) + mol.getImplicitHydrogens(atomIdx);
}

function getElementRank(element) {
  if (element === 'H') return 799999;
  if (element === ',') return 1000;
  if (element === '&') return RING_RANK;
  if (ELEMENT_RANK[element] !== undefined) return ELEMENT_RANK[element];
  return 800000 - (atomicMass(element) || 100);
}

function getBondSymbol(bondType) {
  if (bondType < 0 || bondType > 4) return '';
  return BOND_SYMBOLS[bondType] || '';
}

function padScore(score) {
  let s = String(score);
  while (s.length < 6) s = '0' + s;
  return s;
}

function sortByCanonicalLabel(nodes, mol) {
  nodes.sort((a, b) => {
    if (a.atomIdx < 0 && b.atomIdx < 0) return 0;
    if (a.atomIdx < 0 || a.element === ',') return 0;
    if (b.atomIdx < 0 || b.element === ',') return 0;
    const rankA = mol.getSymmetryRank(a.atomIdx);
    const rankB = mol.getSymmetryRank(b.atomIdx);
    if (rankA < rankB) return -1;
    if (rankA > rankB) return 1;
    return 0;
  });
}

function sortNodesByStringscore(sphere) {
  // CDK bubble sort by stringscore descending
  let changed = true;
  while (changed) {
    changed = false;
    for (let i = 0; i < sphere.length - 1; i++) {
      if (sphere[i + 1].stringscore > sphere[i].stringscore) {
        const tmp = sphere[i + 1];
        sphere[i + 1] = sphere[i];
        sphere[i] = tmp;
        changed = true;
      }
    }
  }
  for (let i = 0; i < sphere.length; i++) {
    sphere[i].sortOrder = sphere.length - i;
  }
}

function bremserElement(element) {
  return BREMSER[element] || element;
}

function atomicMass(element) {
  const masses = {
    H: 1, He: 4, Li: 7, Be: 9, B: 11, C: 12, N: 14, O: 16, F: 19,
    Ne: 20, Na: 23, Mg: 24, Al: 27, Si: 28, P: 31, S: 32, Cl: 35,
    Ar: 40, K: 39, Ca: 40, Fe: 56, Co: 59, Ni: 59, Cu: 64, Zn: 65,
    As: 75, Se: 79, Br: 80, Kr: 84, Ag: 108, Sn: 119, Sb: 122,
    Te: 128, I: 127, Xe: 131, Pt: 195, Au: 197, Hg: 201, Tl: 204,
    Pb: 207, Bi: 209, Ge: 73,
  };
  return masses[element] || 100;
}
