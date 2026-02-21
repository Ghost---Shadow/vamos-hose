import { describe, test, expect } from 'bun:test';
import {
  loadDatabase,
  queryHose,
  computeWeightedAvg,
  extractSolvents,
  clearCache,
} from './database.js';

describe('computeWeightedAvg', () => {
  test('single solvent', () => {
    const entry = {
      n: 'C',
      s: 'CC',
      'C_CDCl3': { min: 10.0, max: 12.0, avg: 11.0, cnt: 5 },
    };
    expect(computeWeightedAvg(entry)).toBe(11.0);
  });

  test('two solvents weighted by count', () => {
    const entry = {
      n: 'C',
      s: 'CC',
      'C_CDCl3': { min: 10.0, max: 12.0, avg: 10.0, cnt: 3 },
      'C_DMSO': { min: 11.0, max: 13.0, avg: 20.0, cnt: 7 },
    };
    // weighted = (10*3 + 20*7) / (3+7) = (30+140)/10 = 17.0
    expect(computeWeightedAvg(entry)).toBe(17.0);
  });

  test('zero total count returns 0', () => {
    const entry = { n: 'C', s: 'CC' };
    expect(computeWeightedAvg(entry)).toBe(0);
  });

  test('rounds to 1 decimal place', () => {
    const entry = {
      n: 'C',
      s: 'CC',
      'C_CDCl3': { min: 1.0, max: 2.0, avg: 10.0, cnt: 3 },
      'C_DMSO': { min: 1.0, max: 2.0, avg: 20.0, cnt: 3 },
    };
    // weighted = (10*3 + 20*3) / 6 = 90/6 = 15.0
    expect(computeWeightedAvg(entry)).toBe(15.0);
  });

  test('skips n and s fields', () => {
    const entry = {
      n: 'C',
      s: 'CCC',
      'C_CDCl3': { min: 5.0, max: 7.0, avg: 6.0, cnt: 2 },
    };
    expect(computeWeightedAvg(entry)).toBe(6.0);
  });
});

describe('extractSolvents', () => {
  test('extracts solvent entries, excludes n and s', () => {
    const entry = {
      n: 'C',
      s: 'CC',
      'C_CDCl3': { min: 10.0, max: 12.0, avg: 11.0, cnt: 5 },
      'C_DMSO': { min: 11.0, max: 13.0, avg: 12.0, cnt: 3 },
    };
    expect(extractSolvents(entry)).toEqual({
      'C_CDCl3': { min: 10.0, max: 12.0, avg: 11.0, cnt: 5 },
      'C_DMSO': { min: 11.0, max: 13.0, avg: 12.0, cnt: 3 },
    });
  });

  test('returns empty object when only metadata present', () => {
    const entry = { n: 'C', s: 'CC' };
    expect(extractSolvents(entry)).toEqual({});
  });
});

describe('queryHose', () => {
  test('returns shift data for known HOSE code', async () => {
    // Use a HOSE code we know exists in the sharded chunks
    const db = await loadDatabase();
    const testCode = Object.keys(db)[0];
    const result = await queryHose(testCode);
    expect(result).not.toBeNull();
    expect(result).toHaveProperty('avgShift');
    expect(result).toHaveProperty('smiles');
    expect(result).toHaveProperty('nucleus');
    expect(result).toHaveProperty('solvents');
    expect(typeof result.avgShift).toBe('number');
  });

  test('returns null for unknown HOSE code', async () => {
    expect(await queryHose('NONEXISTENT_CODE_XYZ')).toBeNull();
  });
});

describe('loadDatabase', () => {
  test('returns an object with HOSE code keys', async () => {
    const db = await loadDatabase();
    expect(typeof db).toBe('object');
    expect(Object.keys(db).length).toBeGreaterThan(0);
  });

  test('entries have n and s fields', async () => {
    const db = await loadDatabase();
    const firstKey = Object.keys(db)[0];
    expect(db[firstKey]).toHaveProperty('n');
    expect(db[firstKey]).toHaveProperty('s');
  });

  test('solvent entries have min, max, avg, cnt', async () => {
    const db = await loadDatabase();
    const firstKey = Object.keys(db)[0];
    const entry = db[firstKey];
    const solventKeys = Object.keys(entry).filter(
      (k) => k !== 'n' && k !== 's',
    );
    expect(solventKeys.length).toBeGreaterThan(0);
    const solvent = entry[solventKeys[0]];
    expect(solvent).toHaveProperty('min');
    expect(solvent).toHaveProperty('max');
    expect(solvent).toHaveProperty('avg');
    expect(solvent).toHaveProperty('cnt');
  });
});
