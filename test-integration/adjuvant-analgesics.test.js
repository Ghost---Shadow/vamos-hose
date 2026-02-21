import { lookupNmrShifts } from '../src/lookup.js';
import {
  gabapentin,
  pregabalin,
  carbamazepine,
  valproicAcid,
} from './adjuvant-analgesics.smiles.js';

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

describe('Adjuvant Analgesics - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // GABAPENTIN
  // SMILES: C1CCC(CC1)(CC(=O)O)CN
  // 9 carbons: cyclohexane (6C), aminomethyl, acetic acid (CH2 + COOH)
  // ------------------------------------------------------------------
  describe('Gabapentin', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(gabapentin, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(9);
    });

    test('cyclohexane carbons around 22-45 ppm', () => {
      const cyclo = result.filter((r) => r.shift >= 20 && r.shift <= 48);
      expect(cyclo.length).toBeGreaterThanOrEqual(5);
    });

    test('carboxylic acid carbon above 170 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 170);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('quaternary carbon around 75-90 ppm', () => {
      const quat = result.filter((r) => r.shift >= 75 && r.shift <= 90);
      expect(quat.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [24.4, 24.4, 27.5, 43.9, 43.9, 51, 62.1, 82.3, 179.7];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // PREGABALIN
  // SMILES: CC(C)CC(CC(=O)O)CN
  // 8 carbons: isobutyl (3C), aminomethyl, CH, acetic acid (CH2 + COOH)
  // ------------------------------------------------------------------
  describe('Pregabalin', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(pregabalin, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(8);
    });

    test('isobutyl methyls around 22-25 ppm', () => {
      const methyls = result.filter((r) => r.shift >= 20 && r.shift <= 28);
      expect(methyls.length).toBeGreaterThanOrEqual(2);
    });

    test('carboxylic acid carbon above 170 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 170);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [23.3, 23.3, 26.4, 33.8, 43.2, 48.1, 73.5, 180.8];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // CARBAMAZEPINE
  // SMILES: C1=CC=C2C(=C1)C=CC3=CC=CC=C3N2C(=O)N
  // 15 carbons: dibenz[b,f]azepine (14C) + carboxamide C=O
  // ------------------------------------------------------------------
  describe('Carbamazepine', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(carbamazepine, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(15);
    });

    test('carboxamide carbon above 155 ppm', () => {
      const amide = result.filter((r) => r.shift >= 155);
      expect(amide.length).toBeGreaterThanOrEqual(1);
    });

    test('aromatic carbons dominate spectrum', () => {
      const aromatic = result.filter((r) => r.shift >= 118 && r.shift <= 160);
      expect(aromatic.length).toBeGreaterThanOrEqual(13);
    });

    test('no aliphatic carbons (all sp2)', () => {
      const aliphatic = result.filter((r) => r.shift < 100);
      expect(aliphatic.length).toBe(0);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [120.9, 120.9, 126.4, 126.4, 128.9, 128.9, 129.1, 129.1, 133.4, 133.4, 138.1, 138.1, 157, 157, 159.6];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // VALPROIC ACID
  // SMILES: CCCC(CCC)C(=O)O
  // 8 carbons: 2-propylpentanoic acid (branched fatty acid)
  // ------------------------------------------------------------------
  describe('Valproic Acid', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(valproicAcid, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(8);
    });

    test('terminal methyl carbons around 14 ppm', () => {
      const methyls = result.filter((r) => r.shift >= 10 && r.shift <= 18);
      expect(methyls.length).toBeGreaterThanOrEqual(2);
    });

    test('carboxylic acid carbon above 178 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 178);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('entirely aliphatic except COOH', () => {
      const aliphatic = result.filter((r) => r.shift < 60);
      expect(aliphatic.length).toBe(7);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [14.6, 14.6, 20.6, 20.6, 34.4, 34.4, 45.2, 183.9];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug: gabapentinoids share GABA-like scaffold
  // ------------------------------------------------------------------
  describe('Gabapentinoid comparison', () => {
    let gabapentinResult, pregabalinResult;

    beforeAll(async () => {
      gabapentinResult = await lookupNmrShifts(gabapentin, { nucleus: '13C' });
      pregabalinResult = await lookupNmrShifts(pregabalin, { nucleus: '13C' });
    });

    test('both have aminomethyl carbon around 40-55 ppm', () => {
      for (const result of [gabapentinResult, pregabalinResult]) {
        const amino = result.filter((r) => r.shift >= 40 && r.shift <= 55);
        expect(amino.length).toBeGreaterThanOrEqual(1);
      }
    });

    test('both have carboxylic acid above 175 ppm', () => {
      for (const result of [gabapentinResult, pregabalinResult]) {
        const cooh = result.filter((r) => r.shift >= 175);
        expect(cooh.length).toBeGreaterThanOrEqual(1);
      }
    });
  });
});
