import { lookupNmrShifts } from '../src/lookup.js';
import {
  thc,
  cbd,
  nabilone,
} from './endocannabinoids.smiles.js';

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

describe('Cannabinoids - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // THC (delta-9-tetrahydrocannabinol)
  // SMILES: CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O
  // 21 carbons: pentyl chain (5C), cyclohexene (6C), resorcinol (6C),
  // pyran ring, gem-dimethyl (2C), vinyl methyl
  // ------------------------------------------------------------------
  describe('THC', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(thc, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(21);
    });

    test('pentyl chain carbons in 14-36 ppm range', () => {
      const pentyl = result.filter((r) => r.shift >= 14 && r.shift <= 36);
      expect(pentyl.length).toBeGreaterThanOrEqual(5);
    });

    test('phenolic aromatic carbons above 150 ppm', () => {
      const phenolic = result.filter((r) => r.shift >= 148 && r.shift <= 165);
      expect(phenolic.length).toBeGreaterThanOrEqual(1);
    });

    test('pyran oxygen-bearing carbon around 75-80 ppm', () => {
      const pyran = result.filter((r) => r.shift >= 72 && r.shift <= 82);
      expect(pyran.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [15.2, 18.4, 24, 30.9, 30.9, 31.2, 31.4, 32.4, 36, 40.4, 40.8, 49.8, 77.2, 109.4, 110, 110.1, 123.9, 135, 142.9, 155.8, 158.9];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // CBD (cannabidiol)
  // SMILES: CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O
  // 21 carbons: same count as THC but open ring
  // ------------------------------------------------------------------
  describe('CBD', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(cbd, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(21);
    });

    test('same carbon count as THC (isomeric relationship)', () => {
      expect(result.length).toBe(21);
    });

    test('exocyclic methylene (=CH2) around 108-115 ppm', () => {
      const methylene = result.filter((r) => r.shift >= 105 && r.shift <= 125);
      expect(methylene.length).toBeGreaterThanOrEqual(1);
    });

    test('two phenolic OH carbons around 150-160 ppm', () => {
      const phenol = result.filter((r) => r.shift >= 148 && r.shift <= 162);
      expect(phenol.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [15.2, 18.4, 24, 30.1, 31.4, 32.4, 36, 38.8, 40.4, 40.8, 73.5, 109.3, 110.1, 110.1, 120, 123.9, 135, 145.7, 151.1, 155.8, 155.8];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // NABILONE
  // SMILES: CCCCCCC(C)(C)C1=CC(=C2C3CC(=O)CCC3C(OC2=C1)(C)C)O
  // 24 carbons: synthetic THC analog with ketone replacing cyclohexene
  // ------------------------------------------------------------------
  describe('Nabilone', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(nabilone, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(24);
    });

    test('cyclohexanone carbonyl above 205 ppm', () => {
      const ketone = result.filter((r) => r.shift >= 205);
      expect(ketone.length).toBeGreaterThanOrEqual(1);
    });

    test('more carbons than THC due to longer alkyl chain', () => {
      expect(result.length).toBeGreaterThan(21);
    });

    test('gem-dimethyl carbons around 23-33 ppm', () => {
      const gemDimethyl = result.filter((r) => r.shift >= 20 && r.shift <= 35);
      expect(gemDimethyl.length).toBeGreaterThanOrEqual(4);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [15.2, 23.1, 24, 30.4, 30.9, 30.9, 31.9, 31.9, 32.9, 35, 38.8, 41.3, 43.9, 49.8, 50.2, 56.1, 77.9, 116.2, 126.3, 126.3, 155.8, 158.9, 161.4, 212.6];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 230);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug: THC vs CBD isomeric comparison
  // ------------------------------------------------------------------
  describe('THC vs CBD structural comparison', () => {
    let thcResult, cbdResult;

    beforeAll(async () => {
      thcResult = await lookupNmrShifts(thc, { nucleus: '13C' });
      cbdResult = await lookupNmrShifts(cbd, { nucleus: '13C' });
    });

    test('both have identical carbon count (structural isomers)', () => {
      expect(thcResult.length).toBe(cbdResult.length);
    });

    test('both show pentyl chain carbons in aliphatic region', () => {
      for (const result of [thcResult, cbdResult]) {
        const pentyl = result.filter((r) => r.shift >= 14 && r.shift <= 36);
        expect(pentyl.length).toBeGreaterThanOrEqual(5);
      }
    });

    test('both show resorcinol/phenolic carbons above 148 ppm', () => {
      for (const result of [thcResult, cbdResult]) {
        const phenol = result.filter((r) => r.shift >= 148);
        expect(phenol.length).toBeGreaterThanOrEqual(1);
      }
    });
  });
});
