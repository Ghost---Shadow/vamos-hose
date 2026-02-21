import { lookupNmrShifts } from '../src/lookup.js';
import {
  acetaminophen,
  phenacetin,
} from './acetaminophen.smiles.js';

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

describe('Acetaminophen and Related - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // ACETAMINOPHEN (PARACETAMOL)
  // SMILES: CC(=O)Nc1ccc(O)cc1
  // 8 carbons: acetyl (CH3 + C=O), para-aminophenol (6C)
  // ------------------------------------------------------------------
  describe('Acetaminophen', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(acetaminophen, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(8);
    });

    test('acetyl methyl around 22-28 ppm', () => {
      const methyl = result.filter((r) => r.shift >= 20 && r.shift <= 30);
      expect(methyl.length).toBeGreaterThanOrEqual(1);
    });

    test('amide carbonyl above 165 ppm', () => {
      const amide = result.filter((r) => r.shift >= 165);
      expect(amide.length).toBeGreaterThanOrEqual(1);
    });

    test('phenolic C-OH around 155-165 ppm', () => {
      const coh = result.filter((r) => r.shift >= 155 && r.shift <= 168);
      expect(coh.length).toBeGreaterThanOrEqual(1);
    });

    test('para-substituted aromatic carbons (two pairs)', () => {
      const aromatic = result.filter((r) => r.shift >= 112 && r.shift <= 140);
      expect(aromatic.length).toBeGreaterThanOrEqual(3);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [25.8, 116.3, 116.3, 122.8, 122.8, 135.1, 161.3, 171.7];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // PHENACETIN
  // SMILES: CC(=O)Nc1ccc(OCC)cc1
  // 10 carbons: acetyl (2C), para-ethoxyaniline (6C), ethoxy (2C)
  // ------------------------------------------------------------------
  describe('Phenacetin', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(phenacetin, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(10);
    });

    test('ethoxy CH3 around 14-18 ppm', () => {
      const ethMe = result.filter((r) => r.shift >= 12 && r.shift <= 20);
      expect(ethMe.length).toBeGreaterThanOrEqual(1);
    });

    test('ethoxy OCH2 around 62-68 ppm', () => {
      const och2 = result.filter((r) => r.shift >= 60 && r.shift <= 70);
      expect(och2.length).toBeGreaterThanOrEqual(1);
    });

    test('two more carbons than acetaminophen (ethyl group)', () => {
      expect(result.length).toBe(10);
    });

    test('amide carbonyl above 165 ppm', () => {
      const amide = result.filter((r) => r.shift >= 165);
      expect(amide.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [15.8, 25.8, 65.5, 115, 115, 123.8, 123.8, 135.1, 159, 171.7];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug: acetaminophen vs phenacetin structural comparison
  // ------------------------------------------------------------------
  describe('Acetaminophen vs Phenacetin comparison', () => {
    let acetResult, phenResult;

    beforeAll(async () => {
      acetResult = await lookupNmrShifts(acetaminophen, { nucleus: '13C' });
      phenResult = await lookupNmrShifts(phenacetin, { nucleus: '13C' });
    });

    test('phenacetin has exactly 2 more carbons (ethyl vs H)', () => {
      expect(phenResult.length - acetResult.length).toBe(2);
    });

    test('both share acetamido group with similar shifts', () => {
      const acetAmide = acetResult.filter((r) => r.shift >= 165);
      const phenAmide = phenResult.filter((r) => r.shift >= 165);
      expect(acetAmide.length).toBeGreaterThanOrEqual(1);
      expect(phenAmide.length).toBeGreaterThanOrEqual(1);
    });

    test('both share para-substituted aromatic pattern', () => {
      for (const result of [acetResult, phenResult]) {
        const aromatic = result.filter((r) => r.shift >= 112 && r.shift <= 140);
        expect(aromatic.length).toBeGreaterThanOrEqual(3);
      }
    });
  });
});
