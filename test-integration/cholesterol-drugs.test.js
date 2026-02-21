import { lookupNmrShifts } from '../src/lookup.js';
import {
  cholesterol,
  fenofibrate,
  gemfibrozil,
} from './cholesterol-drugs.smiles.js';

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

describe('Cholesterol-Lowering Drugs - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // CHOLESTEROL
  // SMILES: CC(CCCC(C)C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C
  // 27 carbons: steroid skeleton (17C), side chain (8C), 2 angular methyls
  // ------------------------------------------------------------------
  describe('Cholesterol', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(cholesterol, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(27);
    });

    test('predominantly aliphatic spectrum', () => {
      const aliphatic = result.filter((r) => r.shift < 80);
      expect(aliphatic.length).toBeGreaterThanOrEqual(24);
    });

    test('olefinic C5=C6 double bond carbons around 120-142 ppm', () => {
      const olefinic = result.filter((r) => r.shift >= 118 && r.shift <= 145);
      expect(olefinic.length).toBe(2);
    });

    test('C3-OH carbon around 68-75 ppm', () => {
      const c3 = result.filter((r) => r.shift >= 66 && r.shift <= 76);
      expect(c3.length).toBeGreaterThanOrEqual(1);
    });

    test('angular methyl carbon around 12-20 ppm', () => {
      const angularMe = result.filter((r) => r.shift >= 12 && r.shift <= 20);
      expect(angularMe.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [19.7, 22, 23.4, 24, 24, 24.8, 28.8, 32.2, 35.2, 36.6, 36.7, 37.1, 38.6, 38.8, 38.8, 40, 40.3, 41.3, 42.4, 42.9, 43.9, 54, 56.9, 57.5, 71.3, 123, 139];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // FENOFIBRATE
  // SMILES: CC(C)OC(=O)C(C)(C)Oc1ccc(cc1)C(=O)c2ccc(Cl)cc2
  // 20 carbons: isopropyl ester (4C), gem-dimethyl (2C), two phenyl
  // rings (12C), ketone C=O, ester C=O
  // ------------------------------------------------------------------
  describe('Fenofibrate', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(fenofibrate, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(20);
    });

    test('ketone carbonyl above 190 ppm', () => {
      const ketone = result.filter((r) => r.shift >= 190);
      expect(ketone.length).toBeGreaterThanOrEqual(1);
    });

    test('ester carbonyl around 172-185 ppm', () => {
      const ester = result.filter((r) => r.shift >= 172 && r.shift <= 185);
      expect(ester.length).toBeGreaterThanOrEqual(1);
    });

    test('aromatic carbons from two phenyl rings', () => {
      const aromatic = result.filter((r) => r.shift >= 125 && r.shift <= 160);
      expect(aromatic.length).toBeGreaterThanOrEqual(8);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [22.9, 22.9, 24.9, 24.9, 71.7, 78.9, 128.7, 128.8, 128.8, 130, 130, 132.1, 132.1, 132.3, 132.3, 136.6, 143.6, 154.9, 178.8, 196.6];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // GEMFIBROZIL
  // SMILES: CC1=CC(C)=CC=C1CCCC(C)(C)C(=O)O
  // 15 carbons: 2,5-dimethylphenyl (8C), butyl chain (4C),
  // gem-dimethyl (2C), COOH
  // ------------------------------------------------------------------
  describe('Gemfibrozil', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(gemfibrozil, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(15);
    });

    test('aryl methyl carbons around 18-23 ppm', () => {
      const arylMe = result.filter((r) => r.shift >= 16 && r.shift <= 24);
      expect(arylMe.length).toBeGreaterThanOrEqual(2);
    });

    test('carboxylic acid carbon above 180 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 180);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('aromatic carbons around 128-148 ppm', () => {
      const aromatic = result.filter((r) => r.shift >= 128 && r.shift <= 148);
      expect(aromatic.length).toBeGreaterThanOrEqual(5);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [19.1, 22, 25, 25, 32, 36.6, 41.5, 42.1, 130.9, 135.5, 136.7, 136.7, 140.8, 144.4, 185.1];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });
});
