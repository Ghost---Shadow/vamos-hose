import { lookupNmrShifts } from '../src/lookup.js';
import {
  lidocaine,
  bupivacaine,
  benzocaine,
  procaine,
} from './local-anesthetics.smiles.js';

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

describe('Local Anesthetics - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // LIDOCAINE
  // SMILES: CCN(CC)CC(=O)NC1=C(C)C=CC=C1C
  // 14 carbons: diethylaminoacetamide (6C), 2,6-dimethylphenyl (8C)
  // ------------------------------------------------------------------
  describe('Lidocaine', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(lidocaine, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(14);
    });

    test('N-ethyl CH3 carbons around 10-15 ppm', () => {
      const ethylMe = result.filter((r) => r.shift >= 8 && r.shift <= 16);
      expect(ethylMe.length).toBeGreaterThanOrEqual(2);
    });

    test('aryl methyl carbons around 17-22 ppm', () => {
      const arylMe = result.filter((r) => r.shift >= 16 && r.shift <= 22);
      expect(arylMe.length).toBeGreaterThanOrEqual(2);
    });

    test('amide carbonyl above 160 ppm', () => {
      const amide = result.filter((r) => r.shift >= 158);
      expect(amide.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [12.6, 12.6, 18.9, 18.9, 51.9, 51.9, 62.2, 128.9, 135.4, 135.4, 135.5, 135.9, 135.9, 163.1];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // BUPIVACAINE
  // SMILES: CCCCN1CCCCC1C(=O)NC2=C(C)C=CC=C2C
  // 18 carbons: piperidine (5C), n-butyl (4C), amide, 2,6-dimethylphenyl (8C)
  // ------------------------------------------------------------------
  describe('Bupivacaine', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(bupivacaine, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(18);
    });

    test('butyl terminal methyl around 12-16 ppm', () => {
      const me = result.filter((r) => r.shift >= 10 && r.shift <= 16);
      expect(me.length).toBeGreaterThanOrEqual(1);
    });

    test('piperidine ring carbons around 20-35 ppm', () => {
      const pip = result.filter((r) => r.shift >= 18 && r.shift <= 35);
      expect(pip.length).toBeGreaterThanOrEqual(4);
    });

    test('amide carbonyl above 165 ppm', () => {
      const amide = result.filter((r) => r.shift >= 165);
      expect(amide.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [14.6, 18.9, 18.9, 20.9, 26.4, 26.5, 28.1, 32.6, 55.7, 60.6, 63.7, 128.9, 135.4, 135.4, 135.5, 135.9, 135.9, 170];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // BENZOCAINE
  // SMILES: CCOC(=O)C1=CC=C(C=C1)N
  // 9 carbons: ethyl ester (2C), phenyl (6C), ester C=O
  // ------------------------------------------------------------------
  describe('Benzocaine', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(benzocaine, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(9);
    });

    test('ester carbonyl above 160 ppm', () => {
      const ester = result.filter((r) => r.shift >= 160);
      expect(ester.length).toBeGreaterThanOrEqual(1);
    });

    test('para-aminophenyl C-N around 148-158 ppm', () => {
      const cn = result.filter((r) => r.shift >= 148 && r.shift <= 160);
      expect(cn.length).toBeGreaterThanOrEqual(1);
    });

    test('ethyl ester CH3 around 14-22 ppm', () => {
      const me = result.filter((r) => r.shift >= 14 && r.shift <= 24);
      expect(me.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [20.9, 65.6, 114.5, 114.5, 130.1, 132.3, 132.3, 155, 167.4];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // PROCAINE
  // SMILES: CCN(CC)CCOC(=O)C1=CC=C(C=C1)N
  // 13 carbons: diethylaminoethyl (6C), ester (1C), para-aminobenzoate (6C)
  // ------------------------------------------------------------------
  describe('Procaine', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(procaine, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(13);
    });

    test('diethylamine CH3 carbons around 10-15 ppm', () => {
      const me = result.filter((r) => r.shift >= 8 && r.shift <= 16);
      expect(me.length).toBeGreaterThanOrEqual(2);
    });

    test('ester carbonyl above 163 ppm', () => {
      const ester = result.filter((r) => r.shift >= 163);
      expect(ester.length).toBeGreaterThanOrEqual(1);
    });

    test('shares para-aminobenzoate pattern with benzocaine', () => {
      const aromatic = result.filter((r) => r.shift >= 110 && r.shift <= 160);
      expect(aromatic.length).toBeGreaterThanOrEqual(4);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [13.1, 13.1, 48.4, 48.4, 51.7, 62.6, 114.5, 114.5, 130.1, 132.3, 132.3, 155, 168.5];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug: amide vs ester anesthetics
  // ------------------------------------------------------------------
  describe('Amide vs Ester structural comparison', () => {
    let lidocaineResult, bupivacaineResult, benzocaineResult, procaineResult;

    beforeAll(async () => {
      lidocaineResult = await lookupNmrShifts(lidocaine, { nucleus: '13C' });
      bupivacaineResult = await lookupNmrShifts(bupivacaine, { nucleus: '13C' });
      benzocaineResult = await lookupNmrShifts(benzocaine, { nucleus: '13C' });
      procaineResult = await lookupNmrShifts(procaine, { nucleus: '13C' });
    });

    test('amide anesthetics (lidocaine, bupivacaine) share 2,6-dimethylphenyl', () => {
      for (const result of [lidocaineResult, bupivacaineResult]) {
        const arylMe = result.filter((r) => r.shift >= 16 && r.shift <= 22);
        expect(arylMe.length).toBeGreaterThanOrEqual(2);
      }
    });

    test('ester anesthetics (benzocaine, procaine) share para-aminobenzoate', () => {
      for (const result of [benzocaineResult, procaineResult]) {
        const paba = result.filter((r) => r.shift >= 148 && r.shift <= 158);
        expect(paba.length).toBeGreaterThanOrEqual(1);
      }
    });
  });
});
