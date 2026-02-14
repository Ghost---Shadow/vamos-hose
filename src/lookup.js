const { smilesToHoseCodes } = require('./smiles-to-hose');
const { loadDatabase, queryHose } = require('./database');

/**
 * Given a SMILES string, returns predicted 13C NMR chemical shifts
 * by converting the SMILES to HOSE codes and looking them up in
 * the preprocessed database.
 *
 * @param {string} smiles - SMILES string of the molecule
 * @param {object} options - { nucleus: '13C' }
 * @returns {Array<{shift: number, atom: string, hose: string, smiles: string}>}
 */
function lookupNmrShifts(smiles, options = {}) {
  const { nucleus = '13C' } = options;

  // Step 1: SMILES -> per-atom HOSE codes
  const hoseCodes = smilesToHoseCodes(smiles, { nucleus });

  // Step 2: Load database (cached after first call)
  const db = loadDatabase();

  // Step 3: Look up each HOSE code and collect shifts
  const results = [];
  for (const entry of hoseCodes) {
    const hit = queryHose(db, entry.hose);
    if (hit) {
      results.push({
        shift: hit.avgShift,
        atom: entry.atom,
        hose: entry.hose,
        smiles: hit.smiles,
      });
    }
  }

  return results;
}

module.exports = { lookupNmrShifts };
