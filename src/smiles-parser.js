/**
 * SMILES Parser
 *
 * Parses a SMILES string into a molecule graph representation.
 *
 * @param {string} smiles - SMILES string
 * @returns {{ atoms: Array<{element: string, charge: number, hydrogens: number, aromatic: boolean}>, bonds: Array<{from: number, to: number, order: number}> }}
 */
function parseSMILES(smiles) {
  // TODO: Implement SMILES parser
  // This is a substantial piece of work. Options:
  // 1. Use openchemlib-js (preferred - battle-tested)
  // 2. Write our own (educational but error-prone)
  // 3. Use smiles-drawer or other npm package
  //
  // For now, throw so tests fail clearly at this layer.
  throw new Error('parseSMILES: Not implemented yet');
}

module.exports = { parseSMILES };
