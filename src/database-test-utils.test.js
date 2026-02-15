import { describe, test, expect } from 'bun:test';
import {
  loadTestDatabase,
  getTestHoseCodes,
  findSimpleTestHoseCode,
} from './database-test-utils.js';
import { queryHose, computeWeightedAvg } from './database.js';

describe('Test Database Utils', () => {
  test('loadTestDatabase returns valid database', () => {
    const db = loadTestDatabase();
    expect(typeof db).toBe('object');
    expect(Object.keys(db).length).toBeGreaterThan(0);
  });

  test('database entries have correct structure', () => {
    const db = loadTestDatabase();
    const keys = Object.keys(db);
    expect(keys.length).toBeGreaterThan(0);

    const firstEntry = db[keys[0]];
    expect(firstEntry).toHaveProperty('n'); // nucleus
    expect(firstEntry).toHaveProperty('s'); // smiles
  });

  test('getTestHoseCodes returns array of codes', () => {
    const codes = getTestHoseCodes();
    expect(Array.isArray(codes)).toBe(true);
    expect(codes.length).toBeGreaterThan(0);
    codes.forEach((code) => {
      expect(typeof code).toBe('string');
    });
  });

  test('findSimpleTestHoseCode returns a valid code', () => {
    const code = findSimpleTestHoseCode();
    expect(typeof code).toBe('string');
    expect(code.length).toBeGreaterThan(0);
    expect(code.length).toBeLessThanOrEqual(20);
  });

  test('findSimpleTestHoseCode filters by nucleus', () => {
    const cCode = findSimpleTestHoseCode({ nucleus: 'C' });
    if (cCode) {
      const db = loadTestDatabase();
      expect(db[cCode].n).toBe('C');
    }
  });

  test('queryHose works with sample database', () => {
    const db = loadTestDatabase();
    const testCode = Object.keys(db)[0];

    const result = queryHose(db, testCode);
    expect(result).not.toBeNull();
    expect(result).toHaveProperty('avgShift');
    expect(result).toHaveProperty('smiles');
    expect(result).toHaveProperty('nucleus');
    expect(result).toHaveProperty('solvents');
    expect(typeof result.avgShift).toBe('number');
  });

  test('computeWeightedAvg works with sample entries', () => {
    const db = loadTestDatabase();
    const testCode = Object.keys(db)[0];
    const entry = db[testCode];

    const avg = computeWeightedAvg(entry);
    expect(typeof avg).toBe('number');
    expect(avg).toBeGreaterThanOrEqual(-100);
    expect(avg).toBeLessThan(300);
  });

  test('sample database is cached', () => {
    const db1 = loadTestDatabase();
    const db2 = loadTestDatabase();
    expect(db1).toBe(db2); // Same reference
  });
});
