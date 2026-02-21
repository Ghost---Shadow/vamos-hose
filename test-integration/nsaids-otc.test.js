import { lookupNmrShifts } from '../src/lookup.js';
import {
  aspirin,
  ibuprofen,
  naproxen,
  ketoprofen,
  diclofenac,
} from './nsaids-otc.smiles.js';

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

describe('OTC NSAIDs - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // ASPIRIN
  // SMILES: CC(=O)Oc1ccccc1C(=O)O
  // 9 carbons: acetyl CH3 + C=O, phenyl ring (6C), COOH
  // ------------------------------------------------------------------
  describe('Aspirin', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(aspirin, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(9);
    });

    test('acetyl methyl around 20 ppm', () => {
      const methyl = result.filter((r) => r.shift >= 18 && r.shift <= 25);
      expect(methyl.length).toBeGreaterThanOrEqual(1);
    });

    test('aromatic carbons in 120-155 ppm range', () => {
      const aromatic = result.filter((r) => r.shift >= 118 && r.shift <= 155);
      expect(aromatic.length).toBeGreaterThanOrEqual(5);
    });

    test('carbonyl carbons above 165 ppm', () => {
      const carbonyls = result.filter((r) => r.shift >= 165);
      expect(carbonyls.length).toBeGreaterThanOrEqual(2);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [21.5, 123.6, 124.9, 125.7, 132.8, 139.2, 151.3, 169.8, 170.4];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // IBUPROFEN
  // SMILES: CC(C)Cc1ccc(cc1)C(C)C(=O)O
  // 13 carbons: isobutyl (4C), phenyl (6C), methyl + CH + COOH
  // ------------------------------------------------------------------
  describe('Ibuprofen', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(ibuprofen, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(13);
    });

    test('aliphatic carbons below 50 ppm', () => {
      const aliphatic = result.filter((r) => r.shift < 50);
      expect(aliphatic.length).toBeGreaterThanOrEqual(5);
    });

    test('aromatic carbons in 125-150 ppm range', () => {
      const aromatic = result.filter((r) => r.shift >= 125 && r.shift <= 150);
      expect(aromatic.length).toBeGreaterThanOrEqual(5);
    });

    test('carboxylic acid carbon above 175 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 175);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [18.8, 23, 23, 30.2, 45.4, 45.7, 127.4, 127.4, 129.7, 129.7, 137.7, 144.9, 181.3];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // NAPROXEN
  // SMILES: COc1ccc2cc(ccc2c1)C(C)C(=O)O
  // 14 carbons: methoxy (1C), naphthalene (10C), methyl + CH + COOH
  // ------------------------------------------------------------------
  describe('Naproxen', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(naproxen, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(14);
    });

    test('methoxy carbon around 55-65 ppm', () => {
      const ome = result.filter((r) => r.shift >= 55 && r.shift <= 65);
      expect(ome.length).toBeGreaterThanOrEqual(1);
    });

    test('naphthalene carbons in aromatic region', () => {
      const aromatic = result.filter((r) => r.shift >= 100 && r.shift <= 165);
      expect(aromatic.length).toBeGreaterThanOrEqual(8);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [18.8, 45.3, 60.8, 106.9, 117.8, 124.7, 127.3, 128.9, 133.8, 134.8, 134.9, 138, 159.7, 181.3];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // KETOPROFEN
  // SMILES: CC(c1cccc(c1)C(=O)c2ccccc2)C(=O)O
  // 16 carbons: methyl + CH (2C), two phenyl rings (12C), ketone + COOH
  // ------------------------------------------------------------------
  describe('Ketoprofen', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(ketoprofen, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(16);
    });

    test('ketone carbonyl above 190 ppm', () => {
      const ketone = result.filter((r) => r.shift >= 188);
      expect(ketone.length).toBeGreaterThanOrEqual(1);
    });

    test('carboxylic acid carbon above 175 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 175 && r.shift <= 190);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('aromatic carbons dominate spectrum', () => {
      const aromatic = result.filter((r) => r.shift >= 125 && r.shift <= 145);
      expect(aromatic.length).toBeGreaterThanOrEqual(10);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [18.8, 45.3, 129.3, 129.4, 130.3, 130.3, 130.7, 130.7, 132.5, 133.4, 138, 138.5, 140.2, 141.1, 181.3, 194.1];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // DICLOFENAC
  // SMILES: OC(=O)Cc1ccccc1Nc2c(Cl)cccc2Cl
  // 14 carbons: acetic acid CH2 + COOH, two phenyl rings (12C)
  // ------------------------------------------------------------------
  describe('Diclofenac', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(diclofenac, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(14);
    });

    test('acetic acid CH2 around 38-48 ppm', () => {
      const ch2 = result.filter((r) => r.shift >= 38 && r.shift <= 48);
      expect(ch2.length).toBeGreaterThanOrEqual(1);
    });

    test('aromatic carbons in 115-150 ppm range', () => {
      const aromatic = result.filter((r) => r.shift >= 115 && r.shift <= 150);
      expect(aromatic.length).toBeGreaterThanOrEqual(10);
    });

    test('carboxylic acid carbon above 175 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 175);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('matches HOSE predictions within 5 ppm', () => {
      const predicted = [43.9, 118.2, 129.2, 130, 130, 130.6, 130.6, 130.6, 131.2, 132, 138.2, 139.7, 144.3, 180.5];
      expectShiftsToMatch(result, predicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug: all OTC NSAIDs share a carboxylic acid group
  // ------------------------------------------------------------------
  describe('Shared COOH pharmacophore', () => {
    let aspirinResult, ibuprofenResult, naproxenResult, ketoprofenResult, diclofenacResult;

    beforeAll(async () => {
      aspirinResult = await lookupNmrShifts(aspirin, { nucleus: '13C' });
      ibuprofenResult = await lookupNmrShifts(ibuprofen, { nucleus: '13C' });
      naproxenResult = await lookupNmrShifts(naproxen, { nucleus: '13C' });
      ketoprofenResult = await lookupNmrShifts(ketoprofen, { nucleus: '13C' });
      diclofenacResult = await lookupNmrShifts(diclofenac, { nucleus: '13C' });
    });

    test('all OTC NSAIDs show COOH carbon in 165-185 ppm range', () => {
      for (const result of [aspirinResult, ibuprofenResult, naproxenResult, ketoprofenResult, diclofenacResult]) {
        const cooh = result.filter((r) => r.shift >= 165 && r.shift <= 185);
        expect(cooh.length).toBeGreaterThanOrEqual(1);
      }
    });

    test('all OTC NSAIDs have aromatic carbons', () => {
      for (const result of [aspirinResult, ibuprofenResult, naproxenResult, ketoprofenResult, diclofenacResult]) {
        const aromatic = result.filter((r) => r.shift >= 115 && r.shift <= 160);
        expect(aromatic.length).toBeGreaterThanOrEqual(4);
      }
    });
  });
});
