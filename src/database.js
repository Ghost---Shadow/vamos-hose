import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let _cachedDb = null;

const DB_PATH = path.join(
  __dirname,
  '..',
  'cleaning-scripts',
  'hose_shift_lookup.json',
);

/**
 * Load the HOSE-keyed shift database.
 * Cached after first call.
 *
 * Database format:
 *   { hoseCode: { "n": nucleus, "s": smiles, solventName: { min, max, avg, cnt } } }
 *
 * @returns {object} the parsed database object
 */
export function loadDatabase() {
  if (_cachedDb) return _cachedDb;
  const raw = fs.readFileSync(DB_PATH, 'utf-8');
  _cachedDb = JSON.parse(raw);
  return _cachedDb;
}

/**
 * Query a single HOSE code against the database.
 *
 * Returns the best shift estimate by averaging across all solvents,
 * weighted by measurement count.
 *
 * @param {object} db - the loaded database
 * @param {string} hoseCode - HOSE code to look up
 * @returns {{ avgShift: number, smiles: string, solvents: object } | null}
 */
export function queryHose(db, hoseCode) {
  const entry = db[hoseCode];
  if (!entry) return null;

  // Weighted average across all solvents
  const avgShift = computeWeightedAvg(entry);

  return {
    avgShift,
    smiles: entry.s,
    nucleus: entry.n,
    solvents: extractSolvents(entry),
  };
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
