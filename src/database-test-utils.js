import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

let _cachedTestDb = null;

const TEST_DB_PATH = path.join(
  __dirname,
  '..',
  'cleaning-scripts',
  'hose_shift_sample.json',
);

/**
 * Load the sample HOSE-keyed shift database for testing.
 * Uses the smaller hose_shift_sample.json instead of the full database.
 *
 * @returns {object} the parsed test database object
 */
export function loadTestDatabase() {
  if (_cachedTestDb) return _cachedTestDb;
  const raw = fs.readFileSync(TEST_DB_PATH, 'utf-8');
  _cachedTestDb = JSON.parse(raw);
  return _cachedTestDb;
}

/**
 * Get a list of known HOSE codes from the test database.
 * Useful for testing without depending on specific molecules.
 *
 * @returns {string[]} array of HOSE codes
 */
export function getTestHoseCodes() {
  const db = loadTestDatabase();
  return Object.keys(db);
}

/**
 * Find a simple HOSE code from the test database.
 * Returns the first HOSE code that matches the criteria.
 *
 * @param {object} options - { maxLength?: number, nucleus?: string }
 * @returns {string | null} a HOSE code or null if none found
 */
export function findSimpleTestHoseCode(options = {}) {
  const { maxLength = 20, nucleus = 'C' } = options;
  const db = loadTestDatabase();

  for (const [hoseCode, entry] of Object.entries(db)) {
    if (entry.n === nucleus && hoseCode.length <= maxLength) {
      return hoseCode;
    }
  }

  return null;
}
