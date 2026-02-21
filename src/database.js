const NUM_CHUNKS = 256;

/**
 * Resolve the base URL for chunks/ relative to this module.
 * Works in both Node.js (file:// URLs) and browsers (http:// URLs).
 */
const CHUNKS_BASE = new URL('../chunks/', import.meta.url).href;

/** Cache of loaded chunks: index -> chunk data object */
const _chunkCache = new Map();

/**
 * Hash a string to a 32-bit integer (must match shard_database.mjs).
 */
function hashCode(str) {
  let h = 0;
  for (let i = 0; i < str.length; i++) {
    h = ((h << 5) - h + str.charCodeAt(i)) | 0;
  }
  return Math.abs(h);
}

/**
 * Get the chunk index for a given HOSE code.
 */
function chunkIndex(hoseCode) {
  return hashCode(hoseCode) % NUM_CHUNKS;
}

/**
 * Load a single chunk by index via dynamic import().
 * Cached after first load.
 *
 * Works in Node.js (file:// URLs) and browsers (http:// URLs)
 * as long as the chunks/ directory is served alongside this module.
 *
 * @param {number} idx - chunk index (0-255)
 * @returns {Promise<object>} the chunk's HOSE-to-entry mapping
 */
async function loadChunk(idx) {
  if (_chunkCache.has(idx)) return _chunkCache.get(idx);

  const chunkName = `chunk_${String(idx).padStart(3, '0')}.js`;
  const chunkUrl = CHUNKS_BASE + chunkName;

  const mod = await import(chunkUrl);
  const data = mod.default;
  _chunkCache.set(idx, data);
  return data;
}

/**
 * Query a single HOSE code against the sharded database.
 *
 * Loads the required chunk on demand, then looks up the HOSE code.
 * Returns the best shift estimate by averaging across all solvents,
 * weighted by measurement count.
 *
 * @param {string} hoseCode - HOSE code to look up
 * @returns {Promise<{ avgShift: number, smiles: string, nucleus: string, solvents: object } | null>}
 */
export async function queryHose(hoseCode) {
  const idx = chunkIndex(hoseCode);
  const chunk = await loadChunk(idx);
  const entry = chunk[hoseCode];
  if (!entry) return null;

  const avgShift = computeWeightedAvg(entry);

  return {
    avgShift,
    smiles: entry.s,
    nucleus: entry.n,
    solvents: extractSolvents(entry),
  };
}

/**
 * Preload chunks for a batch of HOSE codes.
 * Loads all needed chunks in parallel for better performance.
 *
 * @param {string[]} hoseCodes - array of HOSE codes
 * @returns {Promise<void>}
 */
export async function preloadChunks(hoseCodes) {
  const indices = new Set(hoseCodes.map(chunkIndex));
  await Promise.all([...indices].map(loadChunk));
}

/**
 * Load the entire database by loading all 256 chunks.
 * Returns a merged object with all HOSE entries.
 * Provided for backward compatibility with tests.
 *
 * @returns {Promise<object>} the full database as a single object
 */
export async function loadDatabase() {
  const indices = Array.from({ length: NUM_CHUNKS }, (_, i) => i);
  const chunks = await Promise.all(indices.map(loadChunk));
  const merged = {};
  for (const chunk of chunks) {
    Object.assign(merged, chunk);
  }
  return merged;
}

/**
 * Compute weighted average shift across all solvents in a db entry.
 *
 * @param {object} entry - database entry { n, s, solvent1: {min,max,avg,cnt}, ... }
 * @returns {number} weighted average shift
 */
export function computeWeightedAvg(entry) {
  let totalWeight = 0;
  let weightedSum = 0;

  for (const [key, val] of Object.entries(entry)) {
    if (key === 'n' || key === 's') continue;
    weightedSum += val.avg * val.cnt;
    totalWeight += val.cnt;
  }

  if (totalWeight === 0) return 0;
  return Math.round((weightedSum / totalWeight) * 10) / 10;
}

/**
 * Extract solvent data from a db entry into a clean object.
 *
 * @param {object} entry
 * @returns {object} { solventName: { min, max, avg, cnt }, ... }
 */
export function extractSolvents(entry) {
  const solvents = {};
  for (const [key, val] of Object.entries(entry)) {
    if (key === 'n' || key === 's') continue;
    solvents[key] = val;
  }
  return solvents;
}

/**
 * Clear the chunk cache. Useful for testing.
 */
export function clearCache() {
  _chunkCache.clear();
}
