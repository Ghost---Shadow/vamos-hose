import { loadChunk, computeWeightedAvg } from './database.js';
import { MinHeap } from './min-heap.js';

const NUM_PPM_BUCKETS = 250;

/**
 * Resolve the base URL for inverted-index/ relative to this module.
 */
const INDEX_BASE = new URL('../inverted-index/', import.meta.url).href;

/** Cache of loaded PPM files: ppm -> array of [chunkIdx, shift] */
const _ppmCache = new Map();

/**
 * Load a single PPM index file by bucket number.
 * Tries dynamic import() first; falls back to fetch() + Function eval.
 *
 * @param {number} ppm - integer ppm bucket (0-249)
 * @returns {Promise<Array<[number, number]>>} array of [chunkIdx, shift]
 */
async function loadPpmFile(ppm) {
  if (_ppmCache.has(ppm)) return _ppmCache.get(ppm);

  const fileName = `ppm_${String(ppm).padStart(3, '0')}.js`;
  const fileUrl = INDEX_BASE + fileName;

  let data;
  try {
    const mod = await import(fileUrl);
    data = mod.default;
  } catch {
    const res = await fetch(fileUrl);
    if (!res.ok) throw new Error(`Failed to load ${fileName}: ${res.status}`);
    const code = await res.text();
    const match = code.match(/export\s+default\s+/);
    if (!match) throw new Error(`Invalid PPM file format: ${fileName}`);
    // eslint-disable-next-line no-new-func
    data = new Function(`return (${code.slice(match.index + match[0].length)})`)();
  }

  _ppmCache.set(ppm, data);
  return data;
}

/**
 * Clear the PPM index cache. Frees memory after reverse lookup.
 */
export function clearEstimateCache() {
  _ppmCache.clear();
}

/**
 * Reverse lookup: given observed NMR peaks, find the closest HOSE codes
 * in the database using the shift-keyed inverted index.
 *
 * Loads only the relevant PPM index files (~1-2 MB) and uses a min-heap
 * to keep only the top results. Then loads the corresponding chunks to
 * retrieve the actual HOSE codes.
 *
 * @param {object} options
 * @param {string}   [options.nucleus='13C'] - target nucleus (only 13C supported)
 * @param {number[]} options.peaks           - observed chemical shifts (ppm)
 * @param {number}   [options.tolerance=2.0] - ppm tolerance per peak
 * @param {number}   [options.maxResults=10] - max results per peak
 * @returns {Promise<Object<string, Array<{hose: string, shift: number, error: number}>>>}
 *          Map from peak value to array of matching HOSE entries
 */
export async function estimateFromSpectra(options = {}) {
  const {
    peaks = [],
    tolerance = 2.0,
    maxResults = 10,
  } = options;

  if (peaks.length === 0) return {};

  const result = {};

  for (const peak of peaks) {
    // Determine which PPM bucket files to load
    const lo = Math.max(0, Math.floor(peak - tolerance));
    const hi = Math.min(NUM_PPM_BUCKETS - 1, Math.floor(peak + tolerance));

    // Load relevant PPM files in parallel
    const ppmIndices = [];
    for (let p = lo; p <= hi; p++) ppmIndices.push(p);
    const ppmFiles = await Promise.all(ppmIndices.map(loadPpmFile));

    // Use a min-heap (max-heap internally) to keep top N by smallest error
    const heap = new MinHeap(maxResults);

    for (const entries of ppmFiles) {
      for (const [chunkIdx, shift] of entries) {
        const error = Math.abs(shift - peak);
        if (error <= tolerance) {
          heap.push({ key: error, chunkIdx, shift });
        }
      }
    }

    // Resolve top results: load chunks and find matching HOSE codes
    const topEntries = heap.toArray();
    const peakResults = [];

    for (const entry of topEntries) {
      const chunk = await loadChunk(entry.chunkIdx);
      // Find HOSE entries in this chunk whose weighted avg matches the shift
      const hoseCode = findHoseByShift(chunk, entry.shift);
      if (hoseCode) {
        peakResults.push({
          hose: hoseCode,
          shift: entry.shift,
          error: Math.round(entry.key * 1000) / 1000,
        });
      }
    }

    result[peak] = peakResults;
  }

  return result;
}

/**
 * Find the HOSE code in a chunk whose weighted avg shift matches the target.
 *
 * @param {object} chunk - chunk data (hoseCode -> entry)
 * @param {number} targetShift - the exact shift value to find
 * @returns {string|null} the matching HOSE code, or null
 */
function findHoseByShift(chunk, targetShift) {
  for (const [hoseCode, entry] of Object.entries(chunk)) {
    if (entry.n !== 'C') continue;
    const shift = computeWeightedAvg(entry);
    if (shift === targetShift) return hoseCode;
  }
  return null;
}

/**
 * Score all database entries against observed peaks.
 * Kept for backward compatibility with existing tests.
 *
 * @param {object} db - full merged database
 * @param {number[]} peaks
 * @param {number} tolerance
 * @param {number} minMatches
 * @param {string} targetNucleus
 * @returns {Array<{smiles: string, hose: string, matchedPeaks: number, score: number}>}
 */
export function scoreDatabase(db, peaks, tolerance, minMatches, targetNucleus) {
  const bySmiles = new Map();

  for (const [hoseCode, entry] of Object.entries(db)) {
    if (entry.n !== targetNucleus) continue;

    const shift = computeWeightedAvg(entry);

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
