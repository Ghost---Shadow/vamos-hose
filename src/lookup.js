import { smilesToHoseCodes } from './smiles-to-hose.js';
import { loadDatabase, queryHose } from './database.js';

/**
 * Given a SMILES string, returns predicted 13C NMR chemical shifts
 * by converting the SMILES to HOSE codes and looking them up in
 * the preprocessed database.
 *
 * @param {string} smiles - SMILES string of the molecule
 * @param {object} options - { nucleus: '13C' }
 * @returns {Array<{shift: number, atom: string, hose: string, smiles: string}>}
 */
export function lookupNmrShifts(smiles, options = {}) {
  const { nucleus = '13C' } = options;

  // Step 1: SMILES -> per-atom HOSE codes
  const hoseCodes = smilesToHoseCodes(smiles, { nucleus });

  // Step 2: Load database (cached after first call)
  const db = loadDatabase();

  // Step 3: Look up each HOSE code and collect shifts
  // Try progressively shorter codes by truncating spheres
  const results = [];
  for (const entry of hoseCodes) {
    let hit = null;
    let hoseToUse = entry.hose;

    // Try exact match first
    hit = queryHose(db, hoseToUse);

    // If no match, try progressively truncating spheres from the end
    if (!hit) {
      let truncated = hoseToUse;
      // Try removing last sphere (everything from last delimiter onwards)
      for (let attempt = 0; attempt < 3 && !hit; attempt++) {
        const lastDelimIdx = Math.max(
          truncated.lastIndexOf('/'),
          truncated.lastIndexOf(')'),
          truncated.lastIndexOf('('),
        );
        if (lastDelimIdx <= 0) break;

        // Keep the delimiter but remove content after it
        const beforeDelim = truncated.substring(0, lastDelimIdx);
        const delim = truncated[lastDelimIdx];
        truncated = beforeDelim + delim;

        hit = queryHose(db, truncated);
        if (hit) {
          hoseToUse = truncated;
          break;
        }

        // Also try without the delimiter
        truncated = beforeDelim;
        hit = queryHose(db, truncated);
        if (hit) {
          hoseToUse = truncated;
          break;
        }
      }
    }

    // Also try without leading H's
    if (!hit && hoseToUse.match(/^H+/)) {
      const withoutH = hoseToUse.replace(/^H+/, '');
      hit = queryHose(db, withoutH);
      if (hit) hoseToUse = withoutH;
    }

    if (hit) {
      results.push({
        shift: hit.avgShift,
        atom: entry.atom,
        hose: hoseToUse,
        smiles: hit.smiles,
      });
    }
  }

  return results;
}
