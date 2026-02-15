import OCL from 'openchemlib';

/**
 * SMILES Parser
 *
 * Parses a SMILES string into a molecule graph representation.
 *
 * @param {string} smiles - SMILES string
 * @returns {OCL.Molecule} openchemlib Molecule object
 */
export function parseSMILES(smiles) {
  // Use openchemlib's built-in SMILES parser
  return OCL.Molecule.fromSmiles(smiles);
}
