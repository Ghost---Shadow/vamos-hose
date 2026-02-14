/**
 * NMR Shift Lookup
 *
 * Given a SMILES string, returns predicted 13C NMR chemical shifts
 * by converting the SMILES to HOSE codes and looking them up in
 * the preprocessed database.
 *
 * @param {string} smiles - SMILES string of the molecule
 * @param {object} options - { nucleus: '13C' }
 * @returns {Array<{shift: number, atom: string, hose: string, smiles: string}>}
 */
function lookupNmrShifts(smiles, options = {}) {
  // TODO: Implement
  // 1. Parse SMILES into a molecule
  // 2. Generate HOSE codes for each atom
  // 3. Look up each HOSE code in the database
  // 4. Return array of shift predictions
  throw new Error('Not implemented yet');
}

module.exports = { lookupNmrShifts };
