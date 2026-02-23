import OCL from 'openchemlib';

// Elements we may substitute between (atomic numbers)
const SUBSTITUTION_ELEMENTS = [6, 7, 8, 16, 9, 17, 35]; // C, N, O, S, F, Cl, Br

// Elements we may add to a molecule
const ADDABLE_ATOMS = [6, 7, 8]; // C, N, O

/**
 * Clone a molecule by round-tripping through SMILES.
 * @param {string} smiles
 * @returns {{ mol: OCL.Molecule, smiles: string } | null}
 */
export function cloneMolecule(smiles) {
  try {
    const mol = OCL.Molecule.fromSmiles(smiles);
    const canonical = mol.toSmiles();
    return { mol, smiles: canonical };
  } catch {
    return null;
  }
}

/**
 * Apply a mutation function to a molecule, validate, return new SMILES or null.
 * @param {string} smiles - current molecule SMILES
 * @param {function(OCL.Molecule): void} mutateFn - mutation callback
 * @returns {string | null} new canonical SMILES if valid, null if chemically invalid
 */
export function tryMutation(smiles, mutateFn) {
  try {
    const mol = OCL.Molecule.fromSmiles(smiles);
    mutateFn(mol);
    const newSmiles = mol.toSmiles();
    // Re-parse to ensure consistency
    OCL.Molecule.fromSmiles(newSmiles);
    return newSmiles;
  } catch {
    return null;
  }
}

/**
 * Enumerate atom substitution mutations for given atom indices.
 * @param {string} smiles
 * @param {number[]} atomIndices - atom indices to mutate
 * @returns {Array<{ smiles: string, description: string }>}
 */
export function enumerateSubstitutions(smiles, atomIndices) {
  const results = [];
  const mol = OCL.Molecule.fromSmiles(smiles);

  for (const idx of atomIndices) {
    if (idx >= mol.getAllAtoms()) continue;
    const currentNo = mol.getAtomicNo(idx);
    if (currentNo <= 1) continue; // skip H

    for (const targetNo of SUBSTITUTION_ELEMENTS) {
      if (targetNo === currentNo) continue;
      const newSmiles = tryMutation(smiles, (m) => {
        m.setAtomicNo(idx, targetNo);
      });
      if (newSmiles && newSmiles !== smiles) {
        results.push({
          smiles: newSmiles,
          description: `substitute atom ${idx}: Z=${currentNo} -> Z=${targetNo}`,
        });
      }
    }
  }
  return results;
}

/**
 * Enumerate bond order change mutations for bonds adjacent to given atoms.
 * @param {string} smiles
 * @param {number[]} atomIndices
 * @returns {Array<{ smiles: string, description: string }>}
 */
export function enumerateBondChanges(smiles, atomIndices) {
  const results = [];
  const mol = OCL.Molecule.fromSmiles(smiles);
  mol.ensureHelperArrays(OCL.Molecule.cHelperNeighbours || 1);

  const triedBonds = new Set();

  for (const idx of atomIndices) {
    if (idx >= mol.getAllAtoms()) continue;
    const connCount = mol.getConnAtoms(idx);
    for (let i = 0; i < connCount; i++) {
      const bondIdx = mol.getConnBond(idx, i);
      if (triedBonds.has(bondIdx)) continue;
      triedBonds.add(bondIdx);

      const currentOrder = mol.getBondOrder(bondIdx);
      for (const targetOrder of [1, 2, 3]) {
        if (targetOrder === currentOrder) continue;
        const newSmiles = tryMutation(smiles, (m) => {
          m.setBondOrder(bondIdx, targetOrder);
        });
        if (newSmiles && newSmiles !== smiles) {
          results.push({
            smiles: newSmiles,
            description: `bond ${bondIdx}: order ${currentOrder} -> ${targetOrder}`,
          });
        }
      }
    }
  }
  return results;
}

/**
 * Enumerate atom addition mutations. Adds a new atom bonded to atoms with free valence.
 * @param {string} smiles
 * @param {number[]} atomIndices
 * @returns {Array<{ smiles: string, description: string }>}
 */
export function enumerateAdditions(smiles, atomIndices) {
  const results = [];
  const mol = OCL.Molecule.fromSmiles(smiles);
  mol.ensureHelperArrays(OCL.Molecule.cHelperNeighbours || 1);

  for (const idx of atomIndices) {
    if (idx >= mol.getAllAtoms()) continue;
    if (mol.getAtomicNo(idx) <= 1) continue;

    // Check free valence
    const freeVal = mol.getFreeValence(idx);
    if (freeVal <= 0) continue;

    for (const atomicNo of ADDABLE_ATOMS) {
      const newSmiles = tryMutation(smiles, (m) => {
        const newIdx = m.addAtom(atomicNo);
        m.addBond(idx, newIdx);
      });
      if (newSmiles && newSmiles !== smiles) {
        results.push({
          smiles: newSmiles,
          description: `add Z=${atomicNo} to atom ${idx}`,
        });
      }
    }
  }
  return results;
}

/**
 * Enumerate atom removal mutations. Only removes terminal (degree-1) heavy atoms.
 * @param {string} smiles
 * @param {number[]} atomIndices
 * @returns {Array<{ smiles: string, description: string }>}
 */
export function enumerateRemovals(smiles, atomIndices) {
  const results = [];
  const mol = OCL.Molecule.fromSmiles(smiles);
  mol.ensureHelperArrays(OCL.Molecule.cHelperRings || 3);

  const tried = new Set();

  for (const idx of atomIndices) {
    if (idx >= mol.getAllAtoms()) continue;

    // Try removing this atom if terminal
    if (!tried.has(idx) && mol.getAtomicNo(idx) > 1 && mol.getConnAtoms(idx) === 1 && !mol.isRingAtom(idx)) {
      tried.add(idx);
      const newSmiles = tryMutation(smiles, (m) => { m.deleteAtom(idx); });
      if (newSmiles && newSmiles !== smiles) {
        results.push({ smiles: newSmiles, description: `remove terminal atom ${idx}` });
      }
    }

    // Also try removing terminal neighbors
    if (idx >= mol.getAllAtoms()) continue;
    const connCount = mol.getConnAtoms(idx);
    for (let i = 0; i < connCount; i++) {
      const neighbor = mol.getConnAtom(idx, i);
      if (tried.has(neighbor)) continue;
      tried.add(neighbor);
      if (mol.getAtomicNo(neighbor) > 1 && mol.getConnAtoms(neighbor) === 1 && !mol.isRingAtom(neighbor)) {
        const newSmiles = tryMutation(smiles, (m) => { m.deleteAtom(neighbor); });
        if (newSmiles && newSmiles !== smiles) {
          results.push({ smiles: newSmiles, description: `remove terminal neighbor ${neighbor} of atom ${idx}` });
        }
      }
    }
  }
  return results;
}

/**
 * Enumerate all single-step mutations for given atom indices, deduplicated by SMILES.
 * @param {string} smiles
 * @param {number[]} atomIndices - atoms to focus mutations on
 * @returns {Array<{ smiles: string, description: string }>}
 */
export function enumerateAllMutations(smiles, atomIndices) {
  const all = [
    ...enumerateSubstitutions(smiles, atomIndices),
    ...enumerateBondChanges(smiles, atomIndices),
    ...enumerateAdditions(smiles, atomIndices),
    ...enumerateRemovals(smiles, atomIndices),
  ];

  // Deduplicate by canonical SMILES
  const seen = new Set();
  seen.add(smiles); // exclude identity mutation
  const unique = [];
  for (const m of all) {
    if (!seen.has(m.smiles)) {
      seen.add(m.smiles);
      unique.push(m);
    }
  }
  return unique;
}

/**
 * Get all heavy atom indices in a molecule.
 * @param {string} smiles
 * @returns {number[]}
 */
export function getAllHeavyAtomIndices(smiles) {
  const mol = OCL.Molecule.fromSmiles(smiles);
  const indices = [];
  for (let i = 0; i < mol.getAllAtoms(); i++) {
    if (mol.getAtomicNo(i) > 1) indices.push(i);
  }
  return indices;
}

/**
 * Get all carbon atom indices in a molecule.
 * @param {string} smiles
 * @returns {number[]}
 */
export function getAllCarbonIndices(smiles) {
  const mol = OCL.Molecule.fromSmiles(smiles);
  const indices = [];
  for (let i = 0; i < mol.getAllAtoms(); i++) {
    if (mol.getAtomicNo(i) === 6) indices.push(i);
  }
  return indices;
}
