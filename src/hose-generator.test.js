const { generateHoseCode } = require('./hose-generator');

describe('generateHoseCode', () => {
  test('throws not implemented', () => {
    const molecule = {
      atoms: [{ element: 'C' }],
      bonds: [],
    };
    expect(() => generateHoseCode(molecule, 0)).toThrow('Not implemented');
  });

  // Tests below define the expected interface once implemented.
  // Unskip as implementation progresses.

  describe.skip('ethane: C-C', () => {
    const molecule = {
      atoms: [
        { element: 'C', hydrogens: 3, aromatic: false, charge: 0 },
        { element: 'C', hydrogens: 3, aromatic: false, charge: 0 },
      ],
      bonds: [{ from: 0, to: 1, order: 1 }],
    };

    test('HOSE for atom 0 starts with central atom neighbors', () => {
      const hose = generateHoseCode(molecule, 0);
      expect(typeof hose).toBe('string');
      expect(hose.length).toBeGreaterThan(0);
      // Sphere 0 should contain C with 3 H neighbors
      expect(hose).toMatch(/C/);
    });
  });

  describe.skip('propane: C-C-C, center atom', () => {
    const molecule = {
      atoms: [
        { element: 'C', hydrogens: 3, aromatic: false, charge: 0 },
        { element: 'C', hydrogens: 2, aromatic: false, charge: 0 },
        { element: 'C', hydrogens: 3, aromatic: false, charge: 0 },
      ],
      bonds: [
        { from: 0, to: 1, order: 1 },
        { from: 1, to: 2, order: 1 },
      ],
    };

    test('center carbon sees 2 carbons in sphere 0', () => {
      const hose = generateHoseCode(molecule, 1);
      expect(typeof hose).toBe('string');
      // Should have two C atoms in sphere 0 plus H atoms
      expect(hose).toMatch(/C/);
    });
  });

  describe.skip('ethene: C=C', () => {
    const molecule = {
      atoms: [
        { element: 'C', hydrogens: 2, aromatic: false, charge: 0 },
        { element: 'C', hydrogens: 2, aromatic: false, charge: 0 },
      ],
      bonds: [{ from: 0, to: 1, order: 2 }],
    };

    test('HOSE code contains double bond marker', () => {
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('=');
    });
  });

  describe.skip('benzene ring atom', () => {
    const molecule = {
      atoms: Array.from({ length: 6 }, () => ({
        element: 'C',
        hydrogens: 1,
        aromatic: true,
        charge: 0,
      })),
      bonds: [
        { from: 0, to: 1, order: 1.5 },
        { from: 1, to: 2, order: 1.5 },
        { from: 2, to: 3, order: 1.5 },
        { from: 3, to: 4, order: 1.5 },
        { from: 4, to: 5, order: 1.5 },
        { from: 5, to: 0, order: 1.5 },
      ],
    };

    test('HOSE code contains aromatic bond marker', () => {
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('*');
    });

    test('HOSE code contains ring marker', () => {
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('@');
    });
  });

  describe.skip('methanol: C-O', () => {
    const molecule = {
      atoms: [
        { element: 'C', hydrogens: 3, aromatic: false, charge: 0 },
        { element: 'O', hydrogens: 1, aromatic: false, charge: 0 },
      ],
      bonds: [{ from: 0, to: 1, order: 1 }],
    };

    test('carbon HOSE code references oxygen', () => {
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('O');
    });
  });
});
