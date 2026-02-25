import { parseSMILES } from './smiles-parser.js';
import { generateHoseCode } from './hose-generator.js';

/**
 * Convert a SMILES string to a list of per-atom HOSE codes.
 *
 * @param {string} smiles - SMILES string
 * @param {object} options - { nucleus: '13C' }
 * @returns {Array<{atom: string, index: number, hose: string}>}
 */
// Module-level HOSE code cache: avoids redundant openchemlib parsing
const _hoseCache = new Map();

export function smilesToHoseCodes(smiles, options = {}) {
  const { nucleus = '13C' } = options;
  const cacheKey = `${smiles}\0${nucleus}`;
  if (_hoseCache.has(cacheKey)) return _hoseCache.get(cacheKey);

  const targetElement = nucleusToElement(nucleus);

  // Step 1: Parse SMILES into a molecule graph
  const molecule = parseSMILES(smiles);

  // Step 2: For each atom of the target element, generate a HOSE code
  const results = [];
  const numAtoms = molecule.getAllAtoms();

  for (let i = 0; i < numAtoms; i++) {
    const atomLabel = molecule.getAtomLabel(i);
    if (atomLabel !== targetElement) continue;

    const hose = generateHoseCode(molecule, i, { maxSpheres: 4 });
    results.push({
      atom: atomLabel,
      index: i,
      hose,
    });
  }

  _hoseCache.set(cacheKey, results);
  return results;
}

export function clearHoseCache() { _hoseCache.clear(); }

/**
 * Map nucleus string to element symbol.
 * '13C' -> 'C', '1H' -> 'H', '15N' -> 'N', etc.
 *
 * @param {string} nucleus
 * @returns {string}
 */
export function nucleusToElement(nucleus) {
  const match = nucleus.match(/(\d+)([A-Z][a-z]?)/);
  if (match) return match[2];
  // Fallback: strip digits
  return nucleus.replace(/[0-9]/g, '');
}
