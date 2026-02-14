const { parseSMILES } = require('./smiles-parser');
const { generateHoseCode } = require('./hose-generator');

/**
 * Convert a SMILES string to a list of per-atom HOSE codes.
 *
 * @param {string} smiles - SMILES string
 * @param {object} options - { nucleus: '13C' }
 * @returns {Array<{atom: string, index: number, hose: string}>}
 */
function smilesToHoseCodes(smiles, options = {}) {
  const { nucleus = '13C' } = options;
  const targetElement = nucleusToElement(nucleus);

  // Step 1: Parse SMILES into a molecule graph
  const molecule = parseSMILES(smiles);

  // Step 2: For each atom of the target element, generate a HOSE code
  const results = [];
  for (let i = 0; i < molecule.atoms.length; i++) {
    const atom = molecule.atoms[i];
    if (atom.element !== targetElement) continue;

    const hose = generateHoseCode(molecule, i, { maxSpheres: 4 });
    results.push({
      atom: atom.element,
      index: i,
      hose,
    });
  }

  return results;
}

/**
 * Map nucleus string to element symbol.
 * '13C' -> 'C', '1H' -> 'H', '15N' -> 'N', etc.
 *
 * @param {string} nucleus
 * @returns {string}
 */
function nucleusToElement(nucleus) {
  const match = nucleus.match(/(\d+)([A-Z][a-z]?)/);
  if (match) return match[2];
  // Fallback: strip digits
  return nucleus.replace(/[0-9]/g, '');
}

module.exports = { smilesToHoseCodes, nucleusToElement };
