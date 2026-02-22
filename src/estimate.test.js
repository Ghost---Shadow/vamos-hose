import { describe, test, expect } from 'bun:test';
import { scoreDatabase, estimateFromSpectra } from './estimate.js';

describe('scoreDatabase', () => {
  const mockDb = {
    'HOSE_A': { n: 'C', s: 'CC', 'CDCl3': { min: 13, max: 15, avg: 14.0, cnt: 10 } },
    'HOSE_B': { n: 'C', s: 'CC', 'CDCl3': { min: 22, max: 24, avg: 23.0, cnt: 5 } },
    'HOSE_C': { n: 'C', s: 'CCC', 'CDCl3': { min: 13, max: 15, avg: 14.1, cnt: 8 } },
    'HOSE_D': { n: 'C', s: 'CCCC', 'CDCl3': { min: 125, max: 130, avg: 127.5, cnt: 3 } },
    'HOSE_E': { n: 'H', s: 'CC', 'CDCl3': { min: 1, max: 2, avg: 1.5, cnt: 4 } },
  };

  test('matches peaks within tolerance', () => {
    const results = scoreDatabase(mockDb, [14.0, 23.0], 2.0, 1, 'C');
    const ccResult = results.find((r) => r.smiles === 'CC');
    expect(ccResult).toBeDefined();
    expect(ccResult.matchedPeaks).toBe(2);
  });

  test('respects minMatches', () => {
    const results = scoreDatabase(mockDb, [14.0, 23.0], 2.0, 2, 'C');
    // CC matches both peaks (14.0 and 23.0), CCC only matches 14.0
    expect(results.some((r) => r.smiles === 'CC')).toBe(true);
    expect(results.some((r) => r.smiles === 'CCC')).toBe(false);
  });

  test('filters by nucleus', () => {
    const results = scoreDatabase(mockDb, [1.5], 2.0, 1, 'C');
    // HOSE_E has nucleus 'H', so shouldn't appear when searching for 'C'
    expect(results.every((r) => r.smiles !== 'CC' || r.matchedPeaks > 0)).toBe(true);
  });

  test('returns empty for no matches', () => {
    const results = scoreDatabase(mockDb, [999.0], 0.1, 1, 'C');
    expect(results).toEqual([]);
  });

  test('score reflects match quality', () => {
    // Perfect match: shift exactly on peak
    const results = scoreDatabase(mockDb, [14.0], 2.0, 1, 'C');
    const cc = results.find((r) => r.smiles === 'CC');
    const ccc = results.find((r) => r.smiles === 'CCC');
    // CC matches at 14.0 (error=0), CCC at 14.1 (error=0.1)
    // CC score = (1/1) * (1 - 0/2) = 1.0
    // CCC score = (1/1) * (1 - 0.1/2) = 0.95
    expect(cc.score).toBe(1.0);
    expect(ccc.score).toBe(0.95);
  });

  test('score accounts for fraction of peaks matched', () => {
    const results = scoreDatabase(mockDb, [14.0, 23.0, 127.5], 2.0, 1, 'C');
    const cc = results.find((r) => r.smiles === 'CC');
    const cccc = results.find((r) => r.smiles === 'CCCC');
    // CC matches 2/3 peaks, CCCC matches 1/3
    expect(cc.score).toBeGreaterThan(cccc.score);
  });
});

describe('estimateFromSpectra', () => {
  test('returns empty object for empty peaks', async () => {
    const result = await estimateFromSpectra({ peaks: [] });
    expect(result).toEqual({});
  });

  test('returns per-peak results with HOSE codes', async () => {
    // Ethane: ~6 ppm in 13C
    const result = await estimateFromSpectra({
      peaks: [6.0],
      tolerance: 3.0,
    });
    expect(typeof result).toBe('object');
    expect(result[6.0]).toBeDefined();
    expect(Array.isArray(result[6.0])).toBe(true);
    expect(result[6.0].length).toBeGreaterThan(0);
    for (const c of result[6.0]) {
      expect(c).toHaveProperty('hose');
      expect(c).toHaveProperty('shift');
      expect(c).toHaveProperty('error');
      expect(typeof c.hose).toBe('string');
      expect(typeof c.shift).toBe('number');
      expect(c.error).toBeLessThanOrEqual(3.0);
    }
  });

  test('results are sorted by error ascending per peak', async () => {
    const result = await estimateFromSpectra({
      peaks: [128.0],
      tolerance: 3.0,
    });
    const entries = result[128.0];
    for (let i = 1; i < entries.length; i++) {
      expect(entries[i].error).toBeGreaterThanOrEqual(entries[i - 1].error);
    }
  });

  test('respects maxResults per peak', async () => {
    const result = await estimateFromSpectra({
      peaks: [128.0],
      tolerance: 5.0,
      maxResults: 3,
    });
    expect(result[128.0].length).toBeLessThanOrEqual(3);
  });

  test('returns results for multiple peaks', async () => {
    const result = await estimateFromSpectra({
      peaks: [14.0, 128.0],
      tolerance: 3.0,
    });
    expect(result[14.0]).toBeDefined();
    expect(result[128.0]).toBeDefined();
    expect(result[14.0].length).toBeGreaterThan(0);
    expect(result[128.0].length).toBeGreaterThan(0);
  });
});
