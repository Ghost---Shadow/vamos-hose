/**
 * HOSE Code Generator
 *
 * Generates text-format HOSE codes (Bremser notation) from an openchemlib
 * Molecule object. Compatible with nmrshiftdb2 database format.
 *
 * Based on the CDK HOSECodeGenerator algorithm.
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

// Bond type ranking (higher = higher priority)
const BOND_RANK = { 3: 300000, 2: 200000, 1.5: 100000, 1: 0 };

// Bond type symbols: single="", double="=", triple="%", aromatic="*"
const BOND_SYMBOL = { 1: '', 2: '=', 3: '%', 1.5: '*' };

// Sphere delimiters: sphere1="(", sphere2="/", sphere3="/", sphere4=")", sphere5+="/"
const SPHERE_DELIMITERS = ['(', '/', '/', ')', '/', '/', '/', '/', '/', '/'];

// Ring closure bonus for scoring
const RING_RANK = 1100;

/**
 * Generate a text-format HOSE code for a given atom in a molecule.
 *
 * @param {object} mol - openchemlib Molecule instance
 * @param {number} atomIndex - index of the central atom
 * @param {object} options - { maxSpheres: 4 }
 * @returns {string} HOSE code string (part after semicolon, matching DB keys)
 */
function generateHoseCode(mol, atomIndex, options = {}) {
  const { maxSpheres = 4 } = options;

  const numAtoms = mol.getAllAtoms();
  const visited = new Uint8Array(numAtoms);
  visited[atomIndex] = 1;

  // Build BFS tree
  const spheres = []; // spheres[i] = array of TreeNode
  const sphere0 = buildSphere0(mol, atomIndex, visited);
  spheres.push(sphere0);

  for (let s = 1; s < maxSpheres; s++) {
    const prevSphere = spheres[s - 1];
    const nextSphereNodes = [];
    for (const node of prevSphere) {
      if (node.stopper) continue;
      const children = expandNode(mol, node, visited, atomIndex);
      node.children = children;
      nextSphereNodes.push(...children);
    }
    if (nextSphereNodes.length === 0 && s < maxSpheres - 1) {
      spheres.push([]);
      continue;
    }
    spheres.push(nextSphereNodes);
  }

  // Calculate scores for canonical ordering (bottom-up)
  calculateScores(spheres);

  // Sort each sphere
  for (const sphere of spheres) {
    sphere.sort(compareNodes);
  }

  // Generate code string
  return buildHoseString(mol, spheres, maxSpheres, atomIndex);
}

/**
 * TreeNode representing an atom in the HOSE code BFS tree.
 */
function makeNode(atomIdx, element, bondOrder, parent, isRingClosure, implicitH = 0) {
  return {
    atomIdx,
    element,
    bondOrder,
    parent,        // parent TreeNode (null for sphere-0 nodes)
    stopper: isRingClosure || element === 'H', // don't expand ring closures or hydrogens
    implicitH,     // number of implicit hydrogens on this atom
    children: [],
    score: 0,
    subtreeRank: 0,
  };
}

/**
 * Build sphere 0: direct neighbors of the central atom.
 * Note: The center atom itself is NOT included - that's handled separately in buildHoseString.
 */
function buildSphere0(mol, centerIdx, visited) {
  const nodes = [];

  // Add heavy atom neighbors
  const connCount = mol.getConnAtoms(centerIdx);
  for (let i = 0; i < connCount; i++) {
    const neighborIdx = mol.getConnAtom(centerIdx, i);
    const element = mol.getAtomLabel(neighborIdx);

    // Skip explicit hydrogens (already counted as implicit)
    if (element === 'H') continue;

    const bondIdx = mol.getConnBond(centerIdx, i);
    const bondOrder = getBondOrder(mol, bondIdx);

    const isRing = visited[neighborIdx];
    if (!isRing) visited[neighborIdx] = 1;

    const neighborImplicitH = mol.getImplicitHydrogens(neighborIdx);
    nodes.push(makeNode(neighborIdx, element, bondOrder, null, !!isRing, neighborImplicitH));
  }

  return nodes;
}

/**
 * Expand a node: find its unvisited neighbors (skipping H and parent).
 */
function expandNode(mol, parentNode, visited, centerAtomIdx) {
  const atomIdx = parentNode.atomIdx;
  const connCount = mol.getConnAtoms(atomIdx);
  const children = [];

  // Find parent atom index (to skip it)
  // For sphere 0 nodes, parent is null, so we use centerAtomIdx
  const parentAtomIdx = parentNode.parent ? parentNode.parent.atomIdx : centerAtomIdx;

  for (let i = 0; i < connCount; i++) {
    const neighborIdx = mol.getConnAtom(atomIdx, i);

    // Skip the parent atom
    if (neighborIdx === parentAtomIdx) continue;

    const element = mol.getAtomLabel(neighborIdx);
    if (element === 'H') continue;

    const bondIdx = mol.getConnBond(atomIdx, i);
    const bondOrder = getBondOrder(mol, bondIdx);

    const isRingClosure = !!visited[neighborIdx];
    if (!isRingClosure) visited[neighborIdx] = 1;

    const neighborImplicitH = isRingClosure ? 0 : mol.getImplicitHydrogens(neighborIdx);
    children.push(makeNode(neighborIdx, element, bondOrder, parentNode, isRingClosure, neighborImplicitH));
  }

  return children;
}

/**
 * Get bond order, handling aromatic bonds.
 */
function getBondOrder(mol, bondIdx) {
  if (mol.isAromaticBond(bondIdx)) return 1.5;
  return mol.getBondOrder(bondIdx);
}

/**
 * Calculate canonical scores for all nodes (bottom-up).
 */
function calculateScores(spheres) {
  // Score each node: element rank + bond rank
  for (const sphere of spheres) {
    for (const node of sphere) {
      const elemRank = ELEMENT_RANK[node.element] || (800000 - (atomicMass(node.element) || 0));
      const bondRank = BOND_RANK[node.bondOrder] || 0;
      const ringBonus = node.stopper ? RING_RANK : 0;
      node.score = elemRank + bondRank + ringBonus;
    }
  }

  // Accumulate subtree ranks bottom-up
  for (let s = spheres.length - 1; s >= 0; s--) {
    for (const node of spheres[s]) {
      let childRank = 0;
      for (const child of node.children) {
        childRank += child.score + child.subtreeRank;
      }
      node.subtreeRank = childRank;
    }
  }
}

/**
 * Compare two TreeNodes for canonical ordering.
 * First by parent order, then by score + subtree rank.
 */
function compareNodes(a, b) {
  // Same parent? Sort by own score
  if (a.parent === b.parent) {
    return (b.score + b.subtreeRank) - (a.score + a.subtreeRank);
  }
  // Different parents: maintain parent order
  // (parents are already sorted, so we use their position)
  if (a.parent && b.parent) {
    return 0; // stable sort preserves parent order
  }
  return (b.score + b.subtreeRank) - (a.score + a.subtreeRank);
}

/**
 * Build the HOSE code string from sorted spheres.
 * The center atom is NOT in spheres[0] - we need to output it separately.
 */
function buildHoseString(mol, spheres, maxSpheres, centerIdx) {
  let code = '';

  // Output sphere 0: the center atom itself
  const centerImplicitH = mol.getImplicitHydrogens(centerIdx);
  const centerElement = mol.getAtomLabel(centerIdx);
  code += 'H'.repeat(centerImplicitH);
  code += bremserElement(centerElement);

  // Find last non-empty sphere
  let lastNonEmptySphere = -1;
  for (let s = 0; s < spheres.length && s < maxSpheres; s++) {
    if (spheres[s].length > 0) {
      lastNonEmptySphere = s;
    }
  }

  // Output spheres 1+ (spheres[0] = neighbors, spheres[1] = 2nd shell, etc.)
  for (let s = 0; s <= lastNonEmptySphere; s++) {
    // Output sphere delimiter BEFORE sphere content
    if (s < SPHERE_DELIMITERS.length) {
      code += SPHERE_DELIMITERS[s];
    } else {
      code += '/';
    }

    const sphere = spheres[s];
    let lastParent = null;
    let firstInSphere = true;

    for (const node of sphere) {
      // Insert comma between branches (different parents) in same sphere
      if (!firstInSphere && s > 0 && node.parent !== lastParent) {
        code += ',';
      }
      firstInSphere = false;
      lastParent = node.parent;

      // Bond symbol
      const bondSym = BOND_SYMBOL[node.bondOrder] || '';
      // Element symbol (with Bremser substitution)
      // Ring closures use &, hydrogens use H, others use their element
      const isRingClosure = node.stopper && node.element !== 'H';
      const elemSym = isRingClosure ? '&' : bremserElement(node.element);

      // Ring marker: prepend @ if atom is in a ring (sphere 0/1 only)
      if (s === 0 && !node.stopper && isRingAtomInSphere0(node, spheres)) {
        code += '@';
      }

      code += bondSym;

      // For non-hydrogen, non-ring-closure atoms, prepend implicit hydrogens
      if (!isRingClosure && node.element !== 'H') {
        code += 'H'.repeat(node.implicitH);
      }

      code += elemSym;
    }
  }

  // Add trailing delimiter if we output any spheres
  if (lastNonEmptySphere >= 0) {
    const delimiterIdx = lastNonEmptySphere + 1;
    if (delimiterIdx < SPHERE_DELIMITERS.length) {
      code += SPHERE_DELIMITERS[delimiterIdx];
    } else {
      code += '/';
    }
  }

  return code;
}

/**
 * Check if a sphere-0 atom leads to a ring closure in later spheres.
 */
function isRingAtomInSphere0(node, spheres) {
  // Walk down the tree from this node — if any descendant is a ring closure,
  // this sphere-0 atom is part of a ring path.
  const stack = [node];
  while (stack.length > 0) {
    const n = stack.pop();
    if (n.stopper && n !== node) return true;
    for (const child of n.children) {
      stack.push(child);
    }
  }
  return false;
}

/**
 * Convert element symbol to Bremser notation.
 */
function bremserElement(element) {
  return BREMSER[element] || element;
}

/**
 * Rough atomic mass for unknown element ranking.
 */
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

export { generateHoseCode };
