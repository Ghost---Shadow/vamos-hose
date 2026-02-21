import { loadDatabase, computeWeightedAvg } from './database.js';

/**
 * Reverse lookup: given observed NMR peaks, find candidate structures
 * in the database whose HOSE-predicted shifts match the peak pattern.
 *
 * Loads the entire sharded database and scans all entries. This is a
 * heavy operation (~210 MB) intended for structure elucidation workflows.
 *
 * @param {object} options
 * @param {string}   [options.nucleus='13C'] - target nucleus
 * @param {number[]} options.peaks           - observed chemical shifts (ppm)
 * @param {number}   [options.tolerance=2.0] - ppm tolerance per peak
 * @param {number}   [options.minMatches=1]  - minimum peaks that must match
 * @param {number}   [options.maxResults=50] - maximum candidates to return
 * @returns {Promise<Array<{smiles: string, hose: string, matchedPeaks: number, score: number}>>}
 */
export async function estimateFromSpectra(options = {}) {
  const {
    nucleus = '13C',
    peaks = [],
    tolerance = 2.0,
    minMatches = 1,
    maxResults = 50,
  } = options;

  if (peaks.length === 0) return [];

  const targetNucleus = nucleusToShort(nucleus);
  const db = await loadDatabase();

  // Score every HOSE entry against the observed peaks
  const candidates = scoreDatabase(db, peaks, tolerance, minMatches, targetNucleus);

  // Sort by score descending, then by matchedPeaks descending
  candidates.sort((a, b) => b.score - a.score || b.matchedPeaks - a.matchedPeaks);

  return candidates.slice(0, maxResults);
}

/**
 * Score all database entries against observed peaks.
 *
 * For each HOSE entry, compute its weighted-average shift and check
 * how many observed peaks it falls within tolerance of. Group results
 * by SMILES to aggregate matches across multiple HOSE codes from the
 * same molecule.
 *
 * Exported for unit-testing.
 *
 * @param {object} db - full merged database
 * @param {number[]} peaks
 * @param {number} tolerance
 * @param {number} minMatches
 * @param {string} targetNucleus - e.g. 'C'
 * @returns {Array<{smiles: string, hose: string, matchedPeaks: number, score: number}>}
 */
export function scoreDatabase(db, peaks, tolerance, minMatches, targetNucleus) {
  // Build per-SMILES accumulator: smiles -> { hoses[], matchedPeakSet, totalError }
  const bySmiles = new Map();

  for (const [hoseCode, entry] of Object.entries(db)) {
    if (entry.n !== targetNucleus) continue;

    const shift = computeWeightedAvg(entry);

    // Check which observed peaks this shift matches
    for (let i = 0; i < peaks.length; i++) {
      const err = Math.abs(shift - peaks[i]);
      if (err <= tolerance) {
        const smiles = entry.s;
        if (!bySmiles.has(smiles)) {
          bySmiles.set(smiles, {
            hoses: [],
            matchedPeaks: new Set(),
            totalError: 0,
          });
        }
        const acc = bySmiles.get(smiles);
        if (!acc.matchedPeaks.has(i)) {
          acc.matchedPeaks.add(i);
          acc.totalError += err;
          acc.hoses.push(hoseCode);
        }
      }
    }
  }

  // Convert to result array, filtering by minMatches
  const results = [];
  for (const [smiles, acc] of bySmiles) {
    const matched = acc.matchedPeaks.size;
    if (matched < minMatches) continue;

    const avgError = acc.totalError / matched;
    const score = Math.round((matched / peaks.length) * (1 - avgError / tolerance) * 1000) / 1000;

    results.push({
      smiles,
      hose: acc.hoses[0],
      matchedPeaks: matched,
      score,
    });
  }

  return results;
}

/**
 * Map nucleus string to short nucleus letter.
 * '13C' -> 'C', '1H' -> 'H'
 */
function nucleusToShort(nucleus) {
  const match = nucleus.match(/(\d+)([A-Z][a-z]?)/);
  if (match) return match[2];
  return nucleus.replace(/[0-9]/g, '');
}
