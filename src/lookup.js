import { smilesToHoseCodes } from './smiles-to-hose.js';
import { queryHose, preloadChunks } from './database.js';

/**
 * Given a SMILES string, returns predicted 13C NMR chemical shifts
 * by converting the SMILES to HOSE codes and looking them up in
 * the sharded database (chunks loaded on demand).
 *
 * @param {string} smiles - SMILES string of the molecule
 * @param {object} options - { nucleus: '13C' }
 * @returns {Promise<Array<{shift: number, atom: string, hose: string, smiles: string}>>}
 */
export async function lookupNmrShifts(smiles, options = {}) {
  const { nucleus = '13C' } = options;

  // Step 1: SMILES -> per-atom HOSE codes (synchronous, no DB access)
  const hoseCodes = smilesToHoseCodes(smiles, { nucleus });

  // Step 2: Preload chunks for all exact-match HOSE codes in parallel
  await preloadChunks(hoseCodes.map((e) => e.hose));

  // Step 3: Look up each HOSE code with truncation fallback
  const results = [];
  for (const entry of hoseCodes) {
    let hit = null;
    let hoseToUse = entry.hose;

    // Try exact match first
    hit = await queryHose(hoseToUse);

    // If no match, try progressively truncating spheres from the end
    if (!hit) {
      let truncated = hoseToUse;
      for (let attempt = 0; attempt < 8 && !hit; attempt++) {
        const lastDelimIdx = Math.max(
          truncated.lastIndexOf('/'),
          truncated.lastIndexOf(')'),
          truncated.lastIndexOf('('),
        );
        if (lastDelimIdx <= 0) break;

        const beforeDelim = truncated.substring(0, lastDelimIdx);
        const delim = truncated[lastDelimIdx];
        truncated = beforeDelim + delim;

        hit = await queryHose(truncated);
        if (hit) {
          hoseToUse = truncated;
          break;
        }

        truncated = beforeDelim;
        hit = await queryHose(truncated);
        if (hit) {
          hoseToUse = truncated;
          break;
        }
      }
    }

    // Also try without leading H's
    if (!hit && hoseToUse.match(/^H+/)) {
      const withoutH = hoseToUse.replace(/^H+/, '');
      hit = await queryHose(withoutH);
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
