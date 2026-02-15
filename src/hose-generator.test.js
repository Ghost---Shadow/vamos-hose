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

    test('HOSE code contains ring marker', () => {
      const molecule = OCL.Molecule.fromSmiles('c1ccccc1');
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('@');
    });
  });

  describe('methanol: C-O', () => {
    test('carbon HOSE code references oxygen', () => {
      const molecule = OCL.Molecule.fromSmiles('CO');
      const hose = generateHoseCode(molecule, 0);
      expect(hose).toContain('O');
    });
  });
});
