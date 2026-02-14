/**
 * HOSE Code Generator
 *
 * Given a molecule graph and an atom index, generates the HOSE code
 * describing the spherical environment around that atom.
 *
 * HOSE = Hierarchically Ordered Spherical Environment
 * Spheres are separated by '/', branches within a sphere by ','
 * Bond types: (single implicit), = (double), # (triple), * (aromatic)
 *
 * @param {{ atoms: Array, bonds: Array }} molecule - molecule graph
 * @param {number} atomIndex - index of the central atom
 * @param {object} options - { maxSpheres: 4 }
 * @returns {string} HOSE code string
 */
function generateHoseCode(molecule, atomIndex, options = {}) {
  // TODO: Implement HOSE code generation
  // This is the core algorithm:
  // 1. BFS outward from atomIndex in concentric spheres
  // 2. At each sphere, collect atoms and their bond types
  // 3. Sort atoms within each sphere by priority (canonical ordering)
  // 4. Encode as HOSE string: bondType + element for each atom
  // 5. Separate spheres with '/'
  // 6. Mark ring closures with '@'
  //
  // For now, throw so tests fail clearly at this layer.
  throw new Error('generateHoseCode: Not implemented yet');
}

module.exports = { generateHoseCode };
