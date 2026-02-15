import { lookupNmrShifts } from '../src/lookup.js';
import {
  telmisartan,
  losartan,
  valsartan,
  irbesartan,
} from './hypertension-medication.smiles.js';

// Reference data sources:
// - nmrshiftdb2 HOSE code prediction (same engine as our database)
//   https://nmrshiftdb.nmr.uni-koeln.de/
// - Valsartan experimental: PMC2874312, 13C NMR (100 MHz, DMSO-d6)
//   https://pmc.ncbi.nlm.nih.gov/articles/PMC2874312/

// Helper: check that every expected shift has at least one match within
// tolerance in the returned shifts array
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

// Helper: check that no returned shift is wildly outside the expected
// range (sanity check that we're not returning garbage)
function expectShiftsInRange(result, minPpm, maxPpm) {
  for (const r of result) {
    expect(r.shift).toBeGreaterThanOrEqual(minPpm);
    expect(r.shift).toBeLessThanOrEqual(maxPpm);
  }
}

describe('Hypertension Medications - 13C NMR Lookup', () => {
  // ------------------------------------------------------------------
  // LOSARTAN
  // SMILES: CCCCC1=NC(Cl)=C(CO)N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4
  // 22 carbons: butyl chain (4C), imidazole ring (3C), hydroxymethyl (1C),
  // benzyl CH2 (1C), two phenyl rings (12C), tetrazole (1C)
  // ------------------------------------------------------------------
  describe('Losartan', () => {
    let result;

    beforeAll(() => {
      result = lookupNmrShifts(losartan, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(22);
    });

    test('aliphatic carbons in expected range', () => {
      // Butyl chain: ~14, ~23, ~26, ~30 ppm
      const aliphatic = result.filter((r) => r.shift < 70);
      expect(aliphatic.length).toBeGreaterThanOrEqual(6);
    });

    test('aromatic carbons in expected range', () => {
      // Biphenyl + imidazole aromatics: ~120-155 ppm
      const aromatic = result.filter(
        (r) => r.shift >= 115 && r.shift <= 160,
      );
      expect(aromatic.length).toBeGreaterThanOrEqual(12);
    });

    test('tetrazole carbon present around 155 ppm', () => {
      const tetrazole = result.filter(
        (r) => r.shift >= 150 && r.shift <= 165,
      );
      expect(tetrazole.length).toBeGreaterThanOrEqual(1);
    });

    test('hydroxymethyl carbon around 55-65 ppm', () => {
      const ch2oh = result.filter(
        (r) => r.shift >= 50 && r.shift <= 68,
      );
      expect(ch2oh.length).toBeGreaterThanOrEqual(1);
    });

    test('matches nmrshiftdb2 HOSE predictions within 5 ppm', () => {
      const nmrshiftdbPredicted = [
        13.9, 23.4, 25.7, 30.3, 49.1, 60.1, 121.4, 126.8, 127.5,
        127.5, 127.7, 127.7, 129.1, 130.4, 131.2, 135.8, 136.6,
        139.5, 139.6, 145.2, 150.2, 155.0,
      ];
      expectShiftsToMatch(result, nmrshiftdbPredicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // VALSARTAN
  // SMILES: CCCCC(=O)N(CC1=CC=C(C=C1)C2=CC=CC=C2C3=NNN=N3)C(C(C)C)C(=O)O
  // 24 carbons: pentanoyl chain (5C), N-CH2 (1C), two phenyl rings (12C),
  // tetrazole (1C), valine branch (3C+COOH), amide C=O
  // ------------------------------------------------------------------
  describe('Valsartan', () => {
    let result;

    beforeAll(() => {
      result = lookupNmrShifts(valsartan, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(24);
    });

    test('carboxylic acid carbon above 165 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 165);
      // Should have COOH (~170 ppm) and amide C=O (~172 ppm)
      expect(cooh.length).toBeGreaterThanOrEqual(2);
    });

    test('isopropyl methyls around 17-20 ppm', () => {
      const methyls = result.filter(
        (r) => r.shift >= 14 && r.shift <= 22,
      );
      // Two isopropyl methyls + terminal CH3 of pentanoyl
      expect(methyls.length).toBeGreaterThanOrEqual(3);
    });

    test('matches nmrshiftdb2 HOSE predictions within 5 ppm', () => {
      const nmrshiftdbPredicted = [
        13.8, 17.4, 17.4, 22.4, 27.2, 27.7, 35.7, 44.4, 62.9, 121.4,
        126.8, 127.5, 127.5, 127.7, 127.7, 129.1, 130.4, 131.2,
        137.0, 139.5, 139.6, 155.0, 170.2, 172.4,
      ];
      expectShiftsToMatch(result, nmrshiftdbPredicted, 5.0);
    });

    test('matches experimental data within 8 ppm', () => {
      // Experimental: PMC2874312, 13C NMR (100 MHz, DMSO-d6)
      const experimental = [
        174.0, 172.4, 171.8, 141.7, 138.2, 131.5, 131.1, 131.0,
        129.3, 128.8, 128.2, 127.4, 126.7, 70.3, 63.4, 49.9, 32.9,
        28.1, 27.3, 22.2, 20.6, 14.2,
      ];
      // Note: experimental has 22 peaks (some overlap), predicted has 24
      // We check that each experimental peak has a match
      expectShiftsToMatch(result, experimental, 8.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // IRBESARTAN
  // SMILES: CCCCC1=NC2(CCCC2)C(=O)N1CC3=CC=C(C=C3)C4=CC=CC=C4C5=NNN=N5
  // 25 carbons: butyl chain (4C), spirocyclopentane (4C), imidazolinone (3C),
  // benzyl CH2 (1C), two phenyl rings (12C), tetrazole (1C)
  // ------------------------------------------------------------------
  describe('Irbesartan', () => {
    let result;

    beforeAll(() => {
      result = lookupNmrShifts(irbesartan, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(25);
    });

    test('carbonyl carbon above 165 ppm', () => {
      // Imidazolinone C=O (~172 ppm) and C=N (~165 ppm)
      const carbonyl = result.filter((r) => r.shift >= 160);
      expect(carbonyl.length).toBeGreaterThanOrEqual(2);
    });

    test('spirocyclopentane carbons around 25-35 ppm', () => {
      const cyclopentane = result.filter(
        (r) => r.shift >= 20 && r.shift <= 40,
      );
      // 4 cyclopentane CH2 + butyl chain CH2s
      expect(cyclopentane.length).toBeGreaterThanOrEqual(4);
    });

    test('spiro quaternary carbon around 65-75 ppm', () => {
      const spiro = result.filter(
        (r) => r.shift >= 60 && r.shift <= 80,
      );
      expect(spiro.length).toBeGreaterThanOrEqual(1);
    });

    test('matches nmrshiftdb2 HOSE predictions within 5 ppm', () => {
      const nmrshiftdbPredicted = [
        13.8, 22.5, 24.5, 24.5, 28.8, 29.8, 34.9, 34.9, 54.1, 71.3,
        121.4, 126.8, 127.5, 127.5, 127.7, 127.7, 129.1, 130.4,
        131.2, 137.0, 139.5, 139.6, 155.0, 165.4, 172.1,
      ];
      expectShiftsToMatch(result, nmrshiftdbPredicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // TELMISARTAN
  // SMILES: CCCC1=NC2=C(C=C(C=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C(=O)O)C5=NC6=CC=CC=C6N5C)C
  // 33 carbons: propyl chain (3C), bis-benzimidazole (16C), N-methyl (1C),
  // aryl methyl (1C), benzyl CH2 (1C), two phenyl rings (12C), COOH
  // ------------------------------------------------------------------
  describe('Telmisartan', () => {
    let result;

    beforeAll(() => {
      result = lookupNmrShifts(telmisartan, { nucleus: '13C' });
    });

    test('returns correct number of carbon environments', () => {
      expect(result.length).toBe(33);
    });

    test('carboxylic acid carbon above 165 ppm', () => {
      const cooh = result.filter((r) => r.shift >= 165);
      expect(cooh.length).toBeGreaterThanOrEqual(1);
    });

    test('benzimidazole carbons in aromatic region', () => {
      // Bis-benzimidazole: many carbons in 105-155 ppm range
      const aromatic = result.filter(
        (r) => r.shift >= 105 && r.shift <= 160,
      );
      // Benzimidazole (8C each x2) + phenyl rings, minus some overlap
      expect(aromatic.length).toBeGreaterThanOrEqual(20);
    });

    test('N-methyl and aryl methyl around 14-32 ppm', () => {
      const methyls = result.filter(
        (r) => r.shift >= 10 && r.shift <= 35,
      );
      // propyl CH3 (~14), aryl-CH3 (~18), propyl CH2 (~22), N-CH3 (~26, ~32)
      expect(methyls.length).toBeGreaterThanOrEqual(4);
    });

    test('matches nmrshiftdb2 HOSE predictions within 5 ppm', () => {
      const nmrshiftdbPredicted = [
        14.1, 17.8, 21.8, 25.8, 31.6, 49.1, 109.6, 117.0, 119.7,
        122.3, 122.7, 126.3, 127.3, 127.4, 127.5, 127.5, 127.7,
        127.7, 128.3, 129.4, 130.6, 131.2, 132.0, 133.1, 135.8,
        136.5, 141.6, 142.9, 143.4, 143.6, 150.2, 153.6, 174.1,
      ];
      expectShiftsToMatch(result, nmrshiftdbPredicted, 5.0);
    });

    test('all shifts within chemically valid range', () => {
      expectShiftsInRange(result, 0, 220);
    });
  });

  // ------------------------------------------------------------------
  // Cross-drug structural comparisons
  // All 4 drugs share a biphenyl-tetrazole pharmacophore
  // ------------------------------------------------------------------
  describe('Shared biphenyl-tetrazole pharmacophore', () => {
    let losartanResult, valsartanResult, irbesartanResult, telmisartanResult;

    beforeAll(() => {
      losartanResult = lookupNmrShifts(losartan, { nucleus: '13C' });
      valsartanResult = lookupNmrShifts(valsartan, { nucleus: '13C' });
      irbesartanResult = lookupNmrShifts(irbesartan, { nucleus: '13C' });
      // Telmisartan has COOH instead of tetrazole on the same biphenyl
      telmisartanResult = lookupNmrShifts(telmisartan, { nucleus: '13C' });
    });

    test('all ARBs with tetrazole show C5-tetrazole around 150-160 ppm', () => {
      // Losartan, valsartan, irbesartan all have tetrazole
      for (const result of [
        losartanResult,
        valsartanResult,
        irbesartanResult,
      ]) {
        const tetrazoleC = result.filter(
          (r) => r.shift >= 148 && r.shift <= 162,
        );
        expect(tetrazoleC.length).toBeGreaterThanOrEqual(1);
      }
    });

    test('biphenyl carbons cluster around 125-142 ppm in all drugs', () => {
      for (const result of [
        losartanResult,
        valsartanResult,
        irbesartanResult,
        telmisartanResult,
      ]) {
        const biphenyl = result.filter(
          (r) => r.shift >= 125 && r.shift <= 145,
        );
        // At least 8 carbons in the biphenyl region
        expect(biphenyl.length).toBeGreaterThanOrEqual(8);
      }
    });
  });
});
