import { lookupNmrShifts } from '../src/lookup.js';
import {
  cortisone,
  hydrocortisone,
  prednisone,
  dexamethasone,
  budesonide,
} from './steroids.smiles.js';

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

describe('Corticosteroids - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // CORTISONE
  // SMILES: CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C(C(=O)CO)CCC12
  // 21 carbons: steroid skeleton (17C), 2 angular methyls, C=O, COCH2OH
  // ------------------------------------------------------------------
  describe('Cortisone', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(cortisone, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(21);
    });

    test('ketone carbonyls above 195 ppm', () => {
      const ketones = result.filter((r) => r.shift >= 195);
      expect(ketones.length).toBeGreaterThanOrEqual(2);
    });

    test('angular methyl carbons below 25 ppm', () => {
      const methyls = result.filter((r) => r.shift >= 14 && r.shift <= 25);
      expect(methyls.length).toBeGreaterThanOrEqual(2);
    });

    test('enone olefinic carbon around 120-175 ppm', () => {
      const olefinic = result.filter((r) => r.shift >= 118 && r.shift <= 178);
      expect(olefinic.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [17, 23.4, 33, 33.1, 33.9, 35.8, 38.8, 38.8, 39.3, 40.8, 45.5, 48.3, 50.2, 57.7, 60.5, 68.6, 71.9, 124, 171.2, 203.3, 218.1];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 230);
    });
  });

  // ------------------------------------------------------------------
  // HYDROCORTISONE
  // SMILES: CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C(C(=O)CO)CCC12O
  // 21 carbons: same skeleton as cortisone but C11-OH instead of C11-ketone
  // ------------------------------------------------------------------
  describe('Hydrocortisone', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(hydrocortisone, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(21);
    });

    test('C20 ketone above 200 ppm', () => {
      const ketone = result.filter((r) => r.shift >= 200);
      expect(ketone.length).toBeGreaterThanOrEqual(1);
    });

    test('C3 enone carbonyl above 195 ppm', () => {
      const enone = result.filter((r) => r.shift >= 195);
      expect(enone.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [16.7, 23.4, 33.9, 35.8, 36.8, 38.8, 38.8, 39.3, 40.8, 43.5, 49.5, 49.9, 50.2, 60.5, 68.6, 71.9, 85.6, 124, 171.2, 203.3, 218.1];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 230);
    });
  });

  // ------------------------------------------------------------------
  // PREDNISONE
  // SMILES: CC12CC(=O)C=CC1=CC(O)C1C2CCC2(C)C(C(=O)CO)CCC12
  // 21 carbons: steroid with 1,4-diene-3-one A ring
  // ------------------------------------------------------------------
  describe('Prednisone', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(prednisone, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(21);
    });

    test('conjugated diene carbons in 130-180 ppm', () => {
      const diene = result.filter((r) => r.shift >= 130 && r.shift <= 180);
      expect(diene.length).toBeGreaterThanOrEqual(2);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [13.9, 21.1, 23.4, 38.8, 38.8, 39, 40.6, 44.9, 47, 51.8, 55.5, 57.1, 60.5, 68.6, 74.3, 135.4, 138.8, 146, 178.6, 201, 218.1];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 230);
    });
  });

  // ------------------------------------------------------------------
  // DEXAMETHASONE
  // SMILES: CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO
  // 22 carbons: fluorinated steroid
  // ------------------------------------------------------------------
  describe('Dexamethasone', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(dexamethasone, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(22);
    });

    test('side-chain ketone above 200 ppm', () => {
      const ketone = result.filter((r) => r.shift >= 200);
      expect(ketone.length).toBeGreaterThanOrEqual(1);
    });

    test('A-ring enone carbonyl above 180 ppm', () => {
      const enone = result.filter((r) => r.shift >= 180);
      expect(enone.length).toBeGreaterThanOrEqual(1);
    });

    test('16-methyl carbon below 20 ppm', () => {
      const methyl = result.filter((r) => r.shift >= 13 && r.shift <= 20);
      expect(methyl.length).toBeGreaterThanOrEqual(2);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [15.4, 18.4, 33.8, 35.2, 38.8, 40.8, 42.3, 47.3, 47.8, 48.3, 57.1, 67.7, 73.5, 91.6, 100.9, 104.3, 124.3, 131.3, 158.9, 166.8, 188.4, 212.3];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 230);
    });
  });

  // ------------------------------------------------------------------
  // BUDESONIDE
  // SMILES: CCCC1OC2CC3C4CCC5=CC(=O)C=CC5(C)C4(F)C(O)CC3(C)C2(O1)C(=O)CO
  // 25 carbons: steroid with 16,17-acetonide
  // ------------------------------------------------------------------
  describe('Budesonide', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(budesonide, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(25);
    });

    test('more carbons than dexamethasone due to acetonide', () => {
      expect(result.length).toBeGreaterThan(22);
    });

    test('acetal/ketal carbons around 90-110 ppm', () => {
      const acetal = result.filter((r) => r.shift >= 88 && r.shift <= 112);
      expect(acetal.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [15.4, 15.8, 29.3, 32.5, 35.2, 38.8, 40.8, 44.2, 47.3, 47.8, 48.3, 57.1, 67.7, 73.5, 90.2, 91.6, 100.9, 103.9, 104.3, 124.3, 131.3, 158.9, 166.8, 188.4, 212.3];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 230);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug: all steroids share the 4-en-3-one motif
  // ------------------------------------------------------------------
  describe('Shared steroid A-ring enone', () => {
    let cortisoneResult, hydrocortisoneResult, prednisoneResult, dexamethasoneResult;

    beforeAll(async () => {
      cortisoneResult = await lookupNmrShifts(cortisone, { nucleus: '13C' });
      hydrocortisoneResult = await lookupNmrShifts(hydrocortisone, { nucleus: '13C' });
      prednisoneResult = await lookupNmrShifts(prednisone, { nucleus: '13C' });
      dexamethasoneResult = await lookupNmrShifts(dexamethasone, { nucleus: '13C' });
    });

    test('all steroids show olefinic/enone carbons in 118-170 ppm', () => {
      for (const result of [cortisoneResult, hydrocortisoneResult, prednisoneResult, dexamethasoneResult]) {
        const olefinic = result.filter((r) => r.shift >= 118 && r.shift <= 170);
        expect(olefinic.length).toBeGreaterThanOrEqual(1);
      }
    });

    test('all steroids show angular methyl carbons', () => {
      for (const result of [cortisoneResult, hydrocortisoneResult, prednisoneResult, dexamethasoneResult]) {
        const methyls = result.filter((r) => r.shift >= 12 && r.shift <= 25);
        expect(methyls.length).toBeGreaterThanOrEqual(1);
      }
    });
  });
});
