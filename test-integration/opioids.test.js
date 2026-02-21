import { lookupNmrShifts } from '../src/lookup.js';
import {
  morphine,
  codeine,
  fentanyl,
  tramadol,
  methadone,
} from './opioids.smiles.js';

function expectShiftsToMatch(result, expectedShifts, tolerance) {
  const returned = result.map((r) => r.shift);
  const unmatched = [];
  for (const expected of expectedShifts) {
    const found = returned.some(
      (r) => Math.abs(r - expected) <= tolerance,
    );
    if (!found) unmatched.push(expected);
  }
  expect(unmatched).toEqual([]);
}

function expectShiftsInRange(result, minPpm, maxPpm) {
  for (const r of result) {
    expect(r.shift).toBeGreaterThanOrEqual(minPpm);
    expect(r.shift).toBeLessThanOrEqual(maxPpm);
  }
}

describe('Opioid Analgesics - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // MORPHINE
  // SMILES: CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O
  // 17 carbons: morphinan skeleton with phenanthrene ring system
  // ------------------------------------------------------------------
  describe('Morphine', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(morphine, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(17);
    });

    test('N-methyl carbon around 40-48 ppm', () => {
      const nme = result.filter((r) => r.shift >= 40 && r.shift <= 48);
      expect(nme.length).toBeGreaterThanOrEqual(1);
    });

    test('phenolic aromatic carbons in 115-165 ppm', () => {
      const aromatic = result.filter((r) => r.shift >= 115 && r.shift <= 165);
      expect(aromatic.length).toBeGreaterThanOrEqual(5);
    });

    test('oxygenated sp3 carbons around 75-90 ppm', () => {
      const oxy = result.filter((r) => r.shift >= 70 && r.shift <= 92);
      expect(oxy.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [24.4, 35.8, 42.9, 43, 45.7, 47.1, 58.7, 75.9, 86, 118.5, 118.6, 128.5, 130.3, 131, 133.4, 147.2, 160.4];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // CODEINE
  // SMILES: CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O
  // 18 carbons: morphinan + O-methyl
  // ------------------------------------------------------------------
  describe('Codeine', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(codeine, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(18);
    });

    test('O-methyl carbon around 55-65 ppm', () => {
      const ome = result.filter((r) => r.shift >= 55 && r.shift <= 65);
      expect(ome.length).toBeGreaterThanOrEqual(1);
    });

    test('one more carbon than morphine (OMe)', () => {
      // Codeine = morphine + OCH3
      expect(result.length).toBe(18);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [24.4, 35.8, 42.9, 43, 45.7, 47.1, 58.7, 60.4, 75.9, 86, 112.6, 119.3, 128.5, 130.3, 131, 133.4, 159, 160.4];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // FENTANYL
  // SMILES: CCC(=O)N(C1CCN(CC1)CCC2=CC=CC=C2)C3=CC=CC=C3
  // 22 carbons: piperidine (5C), phenethyl (8C), aniline (6C),
  // propanoyl (3C)
  // ------------------------------------------------------------------
  describe('Fentanyl', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(fentanyl, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(22);
    });

    test('amide carbonyl above 170 ppm', () => {
      const amide = result.filter((r) => r.shift >= 170);
      expect(amide.length).toBeGreaterThanOrEqual(1);
    });

    test('piperidine ring carbons around 25-65 ppm', () => {
      const piperidine = result.filter((r) => r.shift >= 25 && r.shift <= 65);
      expect(piperidine.length).toBeGreaterThanOrEqual(5);
    });

    test('two phenyl rings give aromatic carbons', () => {
      const aromatic = result.filter((r) => r.shift >= 120 && r.shift <= 145);
      expect(aromatic.length).toBeGreaterThanOrEqual(10);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [9.7, 32, 33.7, 36.5, 36.5, 60.6, 60.6, 61.3, 64.4, 124.3, 124.3, 128.8, 128.8, 130.5, 130.5, 130.8, 130.8, 135, 138.5, 139.6, 140.7, 178.2];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // TRAMADOL
  // SMILES: CN(C)CC1CCCCC1(C2=CC(=CC=C2)OC)O
  // 16 carbons: cyclohexane (6C), phenyl (6C), dimethylaminomethyl (3C), OMe
  // ------------------------------------------------------------------
  describe('Tramadol', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(tramadol, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(16);
    });

    test('cyclohexane carbons around 20-50 ppm', () => {
      const cyclo = result.filter((r) => r.shift >= 20 && r.shift <= 55);
      expect(cyclo.length).toBeGreaterThanOrEqual(5);
    });

    test('methoxy carbon around 55-65 ppm', () => {
      const ome = result.filter((r) => r.shift >= 55 && r.shift <= 68);
      expect(ome.length).toBeGreaterThanOrEqual(1);
    });

    test('quaternary C-OH around 70-85 ppm', () => {
      const coh = result.filter((r) => r.shift >= 68 && r.shift <= 88);
      expect(coh.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [25.5, 32.8, 38.8, 43.5, 48.4, 48.4, 60.8, 65.2, 73.5, 83.3, 110.8, 122.3, 128.4, 128.6, 150.7, 159.2];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // METHADONE
  // SMILES: CCC(=O)C(CC(C)N(C)C)(C1=CC=CC=C1)C2=CC=CC=C2
  // 21 carbons: diphenyl (12C), propanoyl (3C), dimethylamino (3C),
  // quaternary C, methyl, CH2
  // ------------------------------------------------------------------
  describe('Methadone', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(methadone, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(21);
    });

    test('ketone carbonyl above 210 ppm', () => {
      const ketone = result.filter((r) => r.shift >= 210);
      expect(ketone.length).toBeGreaterThanOrEqual(1);
    });

    test('two phenyl rings in aromatic region', () => {
      const aromatic = result.filter((r) => r.shift >= 125 && r.shift <= 145);
      expect(aromatic.length).toBeGreaterThanOrEqual(10);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [10.1, 11.5, 38.2, 41.1, 41.1, 45.8, 54.5, 82.3, 129.2, 129.2, 129.2, 129.2, 130.4, 130.4, 130.4, 130.4, 138.5, 138.5, 166.6, 166.6, 221.7];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 230);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug: morphine vs codeine structural comparison
  // ------------------------------------------------------------------
  describe('Morphine vs Codeine comparison', () => {
    let morphineResult, codeineResult;

    beforeAll(async () => {
      morphineResult = await lookupNmrShifts(morphine, { nucleus: '13C' });
      codeineResult = await lookupNmrShifts(codeine, { nucleus: '13C' });
    });

    test('codeine has exactly one more carbon than morphine', () => {
      expect(codeineResult.length - morphineResult.length).toBe(1);
    });

    test('both share morphinan ring carbons in 20-50 ppm region', () => {
      const morphAliphatic = morphineResult.filter((r) => r.shift >= 20 && r.shift <= 50);
      const codAliphatic = codeineResult.filter((r) => r.shift >= 20 && r.shift <= 50);
      expect(morphAliphatic.length).toBeGreaterThanOrEqual(4);
      expect(codAliphatic.length).toBeGreaterThanOrEqual(4);
    });
  });
});
