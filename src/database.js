const NUM_CHUNKS = 256;
const MAX_RETRIES = 3;
const RETRY_BASE_MS = 1000;
export const IMPORT_TIMEOUT_MS = 15000;

/**
 * Resolve the base URL for chunks/ relative to this module.
 * Works in both Node.js (file:// URLs) and browsers (http:// URLs).
 */
const CHUNKS_BASE = new URL('../chunks/', import.meta.url).href;

/**
 * Race a promise against a timeout.
 */
export function withTimeout(promise, ms) {
  return Promise.race([
    promise,
    new Promise((_, reject) => setTimeout(() => reject(new Error('timeout')), ms)),
  ]);
}

/**
 * Fetch with retry + exponential backoff for CDN 503s.
 */
export async function fetchWithRetry(url, retries = MAX_RETRIES) {
  for (let attempt = 0; attempt <= retries; attempt++) {
    const res = await fetch(url);
    if (res.ok) return res;
    if (res.status === 503 && attempt < retries) {
      await new Promise(r => setTimeout(r, RETRY_BASE_MS * 2 ** attempt));
      continue;
    }
    throw new Error(`Failed to load ${url}: ${res.status}`);
  }
}

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
export function chunkIndex(hoseCode) {
  return hashCode(hoseCode) % NUM_CHUNKS;
}

/**
 * Load a single chunk by index.
 * Tries dynamic import() with timeout first; falls back to fetch() + Function eval.
 * Retries the entire load up to MAX_RETRIES times on failure.
 * Cached after first load.
 *
 * @param {number} idx - chunk index (0-255)
 * @returns {Promise<object>} the chunk's HOSE-to-entry mapping
 */
export async function loadChunk(idx) {
  if (_chunkCache.has(idx)) return _chunkCache.get(idx);

  const chunkName = `chunk_${String(idx).padStart(3, '0')}.js`;
  const chunkUrl = CHUNKS_BASE + chunkName;

  let lastError;
  for (let attempt = 0; attempt <= MAX_RETRIES; attempt++) {
    if (attempt > 0) {
      await new Promise(r => setTimeout(r, RETRY_BASE_MS * 2 ** (attempt - 1)));
    }
    try {
      let data;
      try {
        const mod = await withTimeout(import(chunkUrl), IMPORT_TIMEOUT_MS);
        data = mod.default;
      } catch {
        // Fallback: fetch as text + parse (for CDN environments where import() fails)
        const res = await fetchWithRetry(chunkUrl, 0);
        const code = await res.text();
        const match = code.match(/export\s+default\s+/);
        if (!match) throw new Error(`Invalid chunk format: ${chunkName}`);
        // eslint-disable-next-line no-new-func
        data = new Function(`return (${code.slice(match.index + match[0].length)})`)();
      }
      _chunkCache.set(idx, data);
      return data;
    } catch (err) {
      lastError = err;
    }
  }
  throw lastError;
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
 * Loads in batches to avoid overwhelming CDN rate limits.
 * Returns a merged object with all HOSE entries.
 *
 * @param {object} [options]
 * @param {function} [options.onProgress] - callback(loaded, total) for progress reporting
 * @returns {Promise<object>} the full database as a single object
 */
export async function loadDatabase(options = {}) {
  const { onProgress } = options;
  // Use smaller batches for HTTP to avoid CDN rate limits; full parallel for local files
  const isHTTP = CHUNKS_BASE.startsWith('http');
  const BATCH_SIZE = isHTTP ? 16 : NUM_CHUNKS;
  const merged = {};

  for (let start = 0; start < NUM_CHUNKS; start += BATCH_SIZE) {
    const end = Math.min(start + BATCH_SIZE, NUM_CHUNKS);
    const batch = Array.from({ length: end - start }, (_, i) => start + i);
    const chunks = await Promise.all(batch.map(loadChunk));
    for (const chunk of chunks) {
      Object.assign(merged, chunk);
    }
    if (onProgress) onProgress(end, NUM_CHUNKS);
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
