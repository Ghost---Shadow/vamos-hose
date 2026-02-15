import { parseSMILES } from './smiles-parser.js';

describe('parseSMILES', () => {
  test('returns openchemlib Molecule object', () => {
    const mol = parseSMILES('C');
    expect(mol).toBeDefined();
    expect(typeof mol.getAllAtoms).toBe('function');
    expect(mol.getAllAtoms()).toBeGreaterThan(0);
  });

  // Tests below define the expected interface once implemented.
  // Unskip as implementation progresses.

  describe.skip('ethane: CC', () => {
    let mol;
    beforeAll(() => { mol = parseSMILES('CC'); });

    test('has 2 atoms', () => {
      expect(mol.atoms.length).toBe(2);
    });

    test('both atoms are carbon', () => {
      expect(mol.atoms[0].element).toBe('C');
      expect(mol.atoms[1].element).toBe('C');
    });

    test('has 1 bond', () => {
      expect(mol.bonds.length).toBe(1);
    });

    test('bond is single order', () => {
      expect(mol.bonds[0].order).toBe(1);
    });
  });

  describe.skip('ethene: C=C', () => {
    let mol;
    beforeAll(() => { mol = parseSMILES('C=C'); });

    test('has 2 atoms', () => {
      expect(mol.atoms.length).toBe(2);
    });

    test('has 1 double bond', () => {
      expect(mol.bonds[0].order).toBe(2);
    });
  });

  describe.skip('benzene: c1ccccc1', () => {
    let mol;
    beforeAll(() => { mol = parseSMILES('c1ccccc1'); });

    test('has 6 atoms', () => {
      expect(mol.atoms.length).toBe(6);
    });

    test('all atoms are aromatic', () => {
      for (const atom of mol.atoms) {
        expect(atom.aromatic).toBe(true);
      }
    });

    test('has 6 bonds', () => {
      expect(mol.bonds.length).toBe(6);
    });
  });

  describe.skip('methanol: CO', () => {
    let mol;
    beforeAll(() => { mol = parseSMILES('CO'); });

    test('has 2 atoms', () => {
      expect(mol.atoms.length).toBe(2);
    });

    test('first atom is C, second is O', () => {
      expect(mol.atoms[0].element).toBe('C');
      expect(mol.atoms[1].element).toBe('O');
    });
  });

  describe.skip('branching: CC(C)C', () => {
    let mol;
    beforeAll(() => { mol = parseSMILES('CC(C)C'); });

    test('has 4 atoms', () => {
      expect(mol.atoms.length).toBe(4);
    });

    test('has 3 bonds', () => {
      expect(mol.bonds.length).toBe(3);
    });

    test('central carbon connects to 3 others', () => {
      const bondsFromAtom1 = mol.bonds.filter(
        (b) => b.from === 1 || b.to === 1,
      );
      expect(bondsFromAtom1.length).toBe(3);
    });
  });

  describe.skip('bracket atom with charge: [NH4+]', () => {
    let mol;
    beforeAll(() => { mol = parseSMILES('[NH4+]'); });

    test('has 1 atom', () => {
      expect(mol.atoms.length).toBe(1);
    });

    test('atom is nitrogen with +1 charge', () => {
      expect(mol.atoms[0].element).toBe('N');
      expect(mol.atoms[0].charge).toBe(1);
    });
  });

  describe.skip('chlorine: CCl', () => {
    let mol;
    beforeAll(() => { mol = parseSMILES('CCl'); });

    test('has 2 atoms', () => {
      expect(mol.atoms.length).toBe(2);
    });

    test('second atom is Cl', () => {
      expect(mol.atoms[1].element).toBe('Cl');
    });
  });

  describe.skip('ring: C1CCCCC1 (cyclohexane)', () => {
    let mol;
    beforeAll(() => { mol = parseSMILES('C1CCCCC1'); });

    test('has 6 atoms', () => {
      expect(mol.atoms.length).toBe(6);
    });

    test('has 6 bonds (ring closed)', () => {
      expect(mol.bonds.length).toBe(6);
    });
  });
});
