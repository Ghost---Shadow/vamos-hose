import { lookupNmrShifts } from '../src/lookup.js';
import {
  celecoxib,
  meloxicam,
  piroxicam,
  etodolac,
  ketorolac,
  nabumetone,
  oxaprozin,
} from './nsaids-prescription.smiles.js';

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

describe('Prescription NSAIDs - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // CELECOXIB
  // SMILES: CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F
  // 17 carbons: tolyl (7C), pyrazole (3C), sulfonamidophenyl (6C), CF3 (1C)
  // ------------------------------------------------------------------
  describe('Celecoxib', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(celecoxib, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(17);
    });

    test('CF3 carbon around 118-125 ppm', () => {
      const cf3 = result.filter((r) => r.shift >= 115 && r.shift <= 128);
      expect(cf3.length).toBeGreaterThanOrEqual(1);
    });

    test('aromatic carbons dominate', () => {
      const aromatic = result.filter((r) => r.shift >= 105 && r.shift <= 150);
      expect(aromatic.length).toBeGreaterThanOrEqual(12);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [23.6, 119.4, 121.6, 129.7, 129.7, 129.9, 129.9, 130.5, 130.5, 131.5, 134.1, 141.9, 143.2, 143.6, 143.6, 144.9, 147.9];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // MELOXICAM
  // SMILES: CC1=C(N=C(S1)NC(=O)C2=C(C3=CC=CC=C3S(=O)(=O)N2C)O)C
  // 15 carbons: thiazole (2C + 2 methyls), benzothiazine ring system (7C),
  // amide C=O, N-CH3, enol C-OH
  // ------------------------------------------------------------------
  describe('Meloxicam', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(meloxicam, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(15);
    });

    test('thiazole methyl carbons below 20 ppm', () => {
      const methyls = result.filter((r) => r.shift >= 10 && r.shift <= 20);
      expect(methyls.length).toBeGreaterThanOrEqual(1);
    });

    test('aromatic/heteroaromatic carbons in 120-170 ppm', () => {
      const aromatic = result.filter((r) => r.shift >= 120 && r.shift <= 170);
      expect(aromatic.length).toBeGreaterThanOrEqual(8);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [13.2, 17.8, 37.7, 123.6, 126.2, 128, 132.4, 135.1, 135.7, 136.1, 138.8, 149.5, 163, 163.9, 164.9];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // PIROXICAM
  // SMILES: CN1C(=C(C2=CC=CC=C2S1(=O)=O)O)C(=O)NC3=CC=CC=N3
  // 15 carbons: N-CH3, benzothiazine (7C), amide C=O, pyridine ring (5C), enol
  // ------------------------------------------------------------------
  describe('Piroxicam', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(piroxicam, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(15);
    });

    test('N-methyl carbon around 35-40 ppm', () => {
      const nme = result.filter((r) => r.shift >= 33 && r.shift <= 42);
      expect(nme.length).toBeGreaterThanOrEqual(1);
    });

    test('pyridine carbons present', () => {
      const pyridine = result.filter((r) => r.shift >= 110 && r.shift <= 155);
      expect(pyridine.length).toBeGreaterThanOrEqual(5);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [37.7, 114.9, 123.6, 123.9, 126.2, 128, 132.4, 135.1, 135.7, 136.1, 138.7, 148.5, 149.5, 151.5, 159.6];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // ETODOLAC
  // SMILES: CCC1=CC2=C(C=C1CC(=O)O)NC3=C2CCOC3(CC)CC
  // 19 carbons: indole ring system (8C), dihydropyran (3C), acetic acid (2C),
  // ethyl groups (6C)
  // ------------------------------------------------------------------
  describe('Etodolac', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(etodolac, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(19);
    });

    test('aliphatic carbons from ethyl groups and ring', () => {
      const aliphatic = result.filter((r) => r.shift < 50);
      expect(aliphatic.length).toBeGreaterThanOrEqual(8);
    });

    test('quaternary spiro-like carbon around 70-80 ppm', () => {
      const quat = result.filter((r) => r.shift >= 68 && r.shift <= 82);
      expect(quat.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [12.2, 12.2, 16.8, 29.7, 39.1, 39.1, 39.2, 41.5, 74.9, 104.3, 117.1, 125, 130.2, 133.2, 142.9, 155.3, 155.3, 177.5, 186.9];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // KETOROLAC
  // SMILES: OC(=O)C1CCN2C1=CC=C2C(=O)C3=CC=CC=C3
  // 15 carbons: pyrrolizine (7C), benzoyl (7C), COOH
  // ------------------------------------------------------------------
  describe('Ketorolac', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(ketorolac, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(15);
    });

    test('ketone carbonyl above 188 ppm', () => {
      const ketone = result.filter((r) => r.shift >= 188);
      expect(ketone.length).toBeGreaterThanOrEqual(1);
    });

    test('carboxylic acid carbon above 175 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 175 && r.shift < 190);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [36.5, 48.9, 73.5, 113.4, 130.3, 130.3, 130.9, 131.8, 131.8, 135.9, 138.3, 138.5, 138.9, 185.6, 195.1];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // NABUMETONE
  // SMILES: COC1=CC2=CC(=CC=C2C=C1)CCC(=O)C
  // 15 carbons: methoxy (1C), naphthalene (10C), propanone chain (3C), methyl
  // ------------------------------------------------------------------
  describe('Nabumetone', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(nabumetone, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(15);
    });

    test('ketone carbonyl above 200 ppm', () => {
      const ketone = result.filter((r) => r.shift >= 200);
      expect(ketone.length).toBeGreaterThanOrEqual(1);
    });

    test('methoxy carbon around 55-65 ppm', () => {
      const ome = result.filter((r) => r.shift >= 55 && r.shift <= 65);
      expect(ome.length).toBeGreaterThanOrEqual(1);
    });

    test('naphthalene aromatic carbons', () => {
      const aromatic = result.filter((r) => r.shift >= 100 && r.shift <= 165);
      expect(aromatic.length).toBeGreaterThanOrEqual(8);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [30.5, 30.8, 45, 60.8, 106.7, 117.8, 125.9, 128.2, 134.9, 136.8, 138, 142.6, 148.7, 159.7, 208.3];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // OXAPROZIN
  // SMILES: OC(=O)CCC1=NC(=C(O1)C2=CC=CC=C2)C3=CC=CC=C3
  // 18 carbons: propionic acid (3C), oxazole (2C), two phenyl rings (12C), COOH
  // ------------------------------------------------------------------
  describe('Oxaprozin', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(oxaprozin, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(18);
    });

    test('aliphatic propionic acid carbons', () => {
      const aliphatic = result.filter((r) => r.shift >= 25 && r.shift <= 45);
      expect(aliphatic.length).toBeGreaterThanOrEqual(2);
    });

    test('aromatic carbons from two phenyl rings', () => {
      const aromatic = result.filter((r) => r.shift >= 125 && r.shift <= 150);
      expect(aromatic.length).toBeGreaterThanOrEqual(10);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [31, 40.3, 126.5, 126.5, 127.9, 127.9, 128, 130, 130, 132.4, 132.4, 132.6, 135.2, 138.5, 138.5, 145.3, 169.3, 183.1];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });
});
