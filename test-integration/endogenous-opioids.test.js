import { lookupNmrShifts } from '../src/lookup.js';
import {
  metEnkephalin,
  leuEnkephalin,
  endomorphin1,
  endomorphin2,
} from './endogenous-opioids.smiles.js';

function expectShiftsInRange(result, minPpm, maxPpm) {
  for (const r of result) {
    expect(r.shift).toBeGreaterThanOrEqual(minPpm);
    expect(r.shift).toBeLessThanOrEqual(maxPpm);
  }
}

describe('Endogenous Opioid Peptides - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // MET-ENKEPHALIN
  // SMILES: CSCCC(C(=O)O)NC(=O)C(CC1=CC=CC=C1)NC(=O)CNC(=O)CNC(=O)C(CC2=CC=C(C=C2)O)N
  // Tyr-Gly-Gly-Phe-Met: 28 carbons
  // ------------------------------------------------------------------
  describe('Met-Enkephalin', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(metEnkephalin, { nucleus: '13C' });
    });

    test('returns expected carbon count for pentapeptide', () => {
      expect(result.length).toBe(27);
    });

    test('multiple amide carbonyls above 165 ppm', () => {
      const amides = result.filter((r) => r.shift >= 165);
      // 4 peptide bonds + 1 terminal COOH = 5
      expect(amides.length).toBeGreaterThanOrEqual(4);
    });

    test('aromatic carbons from Tyr and Phe', () => {
      const aromatic = result.filter((r) => r.shift >= 115 && r.shift <= 160);
      expect(aromatic.length).toBeGreaterThanOrEqual(8);
    });

    test('thioether methyl (S-CH3) around 15-20 ppm', () => {
      const sme = result.filter((r) => r.shift >= 12 && r.shift <= 22);
      expect(sme.length).toBeGreaterThanOrEqual(1);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // LEU-ENKEPHALIN
  // SMILES: CC(C)CC(C(=O)O)NC(=O)C(CC1=CC=CC=C1)NC(=O)CNC(=O)CNC(=O)C(CC2=CC=C(C=C2)O)N
  // Tyr-Gly-Gly-Phe-Leu: 29 carbons
  // ------------------------------------------------------------------
  describe('Leu-Enkephalin', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(leuEnkephalin, { nucleus: '13C' });
    });

    test('returns expected carbon count', () => {
      expect(result.length).toBe(28);
    });

    test('one more carbon than met-enkephalin (Leu vs Met)', () => {
      // Leu has isobutyl (4C) vs Met thioether (3C)
      expect(result.length).toBe(28);
    });

    test('leucine isopropyl methyls around 22-25 ppm', () => {
      const methyls = result.filter((r) => r.shift >= 18 && r.shift <= 28);
      expect(methyls.length).toBeGreaterThanOrEqual(2);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // ENDOMORPHIN-1
  // SMILES: C1CC(N(C1)C(=O)C(CC2=CC=C(C=C2)O)N)C(=O)NC(CC3=CNC4=CC=CC=C43)C(=O)NC(CC5=CC=CC=C5)C(=O)N
  // Tyr-Pro-Trp-Phe-NH2: 33 carbons
  // ------------------------------------------------------------------
  describe('Endomorphin-1', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(endomorphin1, { nucleus: '13C' });
    });

    test('returns expected carbon count for tetrapeptide', () => {
      expect(result.length).toBe(34);
    });

    test('tryptophan indole carbons in 105-140 ppm', () => {
      const indole = result.filter((r) => r.shift >= 105 && r.shift <= 140);
      expect(indole.length).toBeGreaterThanOrEqual(6);
    });

    test('proline ring carbons around 25-50 ppm', () => {
      const pro = result.filter((r) => r.shift >= 22 && r.shift <= 52);
      expect(pro.length).toBeGreaterThanOrEqual(3);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // ENDOMORPHIN-2
  // SMILES: C1CC(N(C1)C(=O)C(CC2=CC=C(C=C2)O)N)C(=O)NC(CC3=CC=CC=C3)C(=O)NC(CC4=CC=CC=C4)C(=O)N
  // Tyr-Pro-Phe-Phe-NH2: 31 carbons
  // ------------------------------------------------------------------
  describe('Endomorphin-2', () => {
    let result;

    beforeAll(async () => {
      result = await lookupNmrShifts(endomorphin2, { nucleus: '13C' });
    });

    test('returns expected carbon count', () => {
      expect(result.length).toBe(32);
    });

    test('fewer carbons than endomorphin-1 (Phe vs Trp)', () => {
      expect(result.length).toBeLessThan(34);
    });

    test('three phenyl rings give many aromatic carbons', () => {
      const aromatic = result.filter((r) => r.shift >= 115 && r.shift <= 160);
      expect(aromatic.length).toBeGreaterThanOrEqual(12);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug: enkephalin comparison
  // ------------------------------------------------------------------
  describe('Enkephalin comparison', () => {
    let metResult, leuResult;

    beforeAll(async () => {
      metResult = await lookupNmrShifts(metEnkephalin, { nucleus: '13C' });
      leuResult = await lookupNmrShifts(leuEnkephalin, { nucleus: '13C' });
    });

    test('both share Tyr-Gly-Gly-Phe core', () => {
      // Both should have similar aromatic region from Tyr + Phe
      const metAromatic = metResult.filter((r) => r.shift >= 115 && r.shift <= 160);
      const leuAromatic = leuResult.filter((r) => r.shift >= 115 && r.shift <= 160);
      expect(metAromatic.length).toBeGreaterThanOrEqual(8);
      expect(leuAromatic.length).toBeGreaterThanOrEqual(8);
    });

    test('both have multiple peptide bond carbonyls', () => {
      for (const result of [metResult, leuResult]) {
        const carbonyls = result.filter((r) => r.shift >= 165);
        expect(carbonyls.length).toBeGreaterThanOrEqual(4);
      }
    });
  });
});
