import { generateHoseCode } from './hose-generator.js';
import OCL from 'openchemlib';
describe('generateHoseCode', () => {
  test('generates HOSE code for ethane', () => {
    const molecule = OCL.Molecule.fromSmiles('CC');
    const hose = generateHoseCode(molecule, 0);
    expect(typeof hose).toBe('string');
    expect(hose.length).toBeGreaterThan(0);
    expect(hose).toContain('C');
    expect(hose).toContain('(');
  });

  // Tests below define the expected interface once implemented.
  // Unskip as implementation progresses.

  describe('ethane: C-C', () => {
    test('HOSE for atom 0 contains carbon', () => {
      const molecule = OCL.Molecule.fromSmiles('CC');
      const hose = generateHoseCode(molecule, 0);
      expect(typeof hose).toBe('string');
      expect(hose.length).toBeGreaterThan(0);
      expect(hose).toMatch(/C/);
    });
  });

  describe('propane: C-C-C, center atom', () => {
    test('center carbon generates valid HOSE code', () => {
      const molecule = OCL.Molecule.fromSmiles('CCC');
      const hose = generateHoseCode(molecule, 1);
      expect(typeof hose).toBe('string');
      expect(hose).toMatch(/C/);
    });
  });

  describe('ethene: C=C', () => {
    test('HOSE code contains double bond marker', () => {
      const molecule = OCL.Molecule.fromSmiles('C=C');
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('=');
    });
  });

  describe('benzene ring atom', () => {
    test('HOSE code contains aromatic bond marker', () => {
      const molecule = OCL.Molecule.fromSmiles('c1ccccc1');
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('*');
    });

    test('HOSE code contains ring closure marker', () => {
      const molecule = OCL.Molecule.fromSmiles('c1ccccc1');
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('&');
    });
  });

  describe('methanol: C-O', () => {
    test('carbon HOSE code references oxygen', () => {
      const molecule = OCL.Molecule.fromSmiles('CO');
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('O');
    });
  });

  describe('exact HOSE code matching against nmrshiftdb2 format', () => {
    test('benzene produces canonical HOSE code', () => {
      const mol = OCL.Molecule.fromSmiles('c1ccccc1');
      const hose = generateHoseCode(mol, 0);
      expect(hose).toBe('H*C*C(H,H,*C,*C/H,H,*C,*&/H*&)');
    });

    test('benzene produces identical codes for all atoms', () => {
      const mol = OCL.Molecule.fromSmiles('c1ccccc1');
      const codes = [];
      for (let i = 0; i < mol.getAllAtoms(); i++) {
        codes.push(generateHoseCode(mol, i));
      }
      expect(new Set(codes).size).toBe(1);
    });

    test('propane produces correct HOSE codes', () => {
      const mol = OCL.Molecule.fromSmiles('CCC');
      expect(generateHoseCode(mol, 0)).toBe('HHHC(HHC/HHH/)');
      expect(generateHoseCode(mol, 1)).toBe('HHCC(HHH,HHH//)');
      expect(generateHoseCode(mol, 2)).toBe('HHHC(HHC/HHH/)');
    });

    test('acetone produces correct HOSE codes', () => {
      const mol = OCL.Molecule.fromSmiles('CC(=O)C');
      expect(generateHoseCode(mol, 0)).toBe('HHHC(=OC/,HHH/)');
      expect(generateHoseCode(mol, 1)).toBe('=OCC(,HHH,HHH//)');
    });

    test('toluene produces correct HOSE codes for all atoms', () => {
      const mol = OCL.Molecule.fromSmiles('Cc1ccccc1');
      expect(generateHoseCode(mol, 0)).toBe('HHHC(*C*C/H,H,*C,*C/H,H,*C,*&)');
      expect(generateHoseCode(mol, 4)).toBe('H*C*C(H,H,*C,*C/H,H,*C,*&/*&C)');
    });

    test('naphthalene produces matching HOSE codes', () => {
      const mol = OCL.Molecule.fromSmiles('c1ccc2ccccc2c1');
      const n = mol.getAllAtoms();
      let hits = 0;
      for (let i = 0; i < n; i++) {
        const code = generateHoseCode(mol, i);
        if (code.includes('&')) hits++;
      }
      expect(hits).toBe(n);
    });

    test('cyclohexane produces correct HOSE code', () => {
      const mol = OCL.Molecule.fromSmiles('C1CCCCC1');
      const hose = generateHoseCode(mol, 0);
      expect(hose).toBe('HHCC(HH,HH,C,C/HH,HH,C,&/HH&)');
    });

    test('symmetric atoms produce identical HOSE codes', () => {
      const mol = OCL.Molecule.fromSmiles('Cc1ccccc1');
      const code2 = generateHoseCode(mol, 2);
      const code6 = generateHoseCode(mol, 6);
      expect(code2).toBe(code6);
      const code3 = generateHoseCode(mol, 3);
      const code5 = generateHoseCode(mol, 5);
      expect(code3).toBe(code5);
    });
  });
});
