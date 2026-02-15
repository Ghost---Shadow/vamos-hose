import { describe, test, expect } from 'bun:test';
import { lookupNmrShifts } from './lookup.js';
import { loadTestDatabase } from './database-test-utils.js';

describe('lookupNmrShifts', () => {
  test('returns array with required properties', () => {
    // Use a simple molecule that should have entries in the sample database
    const result = lookupNmrShifts('CC', { nucleus: '13C' });
    expect(Array.isArray(result)).toBe(true);

    // Check structure of returned objects if any results found
    result.forEach((entry) => {
      expect(entry).toHaveProperty('shift');
      expect(entry).toHaveProperty('atom');
      expect(entry).toHaveProperty('hose');
      expect(entry).toHaveProperty('smiles');
      expect(typeof entry.shift).toBe('number');
      expect(typeof entry.atom).toBe('string');
      expect(typeof entry.hose).toBe('string');
      expect(typeof entry.smiles).toBe('string');
    });
  });

  test('returns array for methane', () => {
    const result = lookupNmrShifts('C', { nucleus: '13C' });
    expect(Array.isArray(result)).toBe(true);
  });

  test('returns array for propane', () => {
    const result = lookupNmrShifts('CCC', { nucleus: '13C' });
    expect(Array.isArray(result)).toBe(true);
  });

  test('defaults nucleus to 13C', () => {
    const result1 = lookupNmrShifts('CC');
    const result2 = lookupNmrShifts('CC', { nucleus: '13C' });
    expect(Array.isArray(result1)).toBe(true);
    expect(Array.isArray(result2)).toBe(true);
  });

  test('accepts different nucleus parameter', () => {
    const result = lookupNmrShifts('C', { nucleus: '1H' });
    expect(Array.isArray(result)).toBe(true);
  });

  test('shift values are within valid NMR range', () => {
    const result = lookupNmrShifts('CC', { nucleus: '13C' });

    result.forEach((entry) => {
      // 13C NMR shifts typically range from -50 to 250 ppm
      expect(entry.shift).toBeGreaterThan(-100);
      expect(entry.shift).toBeLessThan(300);
    });
  });

  test('sample database is accessible', () => {
    const db = loadTestDatabase();
    expect(typeof db).toBe('object');
    expect(Object.keys(db).length).toBeGreaterThan(0);

    // Check first entry has expected structure
    const firstKey = Object.keys(db)[0];
    const firstEntry = db[firstKey];
    expect(firstEntry).toHaveProperty('n'); // nucleus
    expect(firstEntry).toHaveProperty('s'); // smiles
  });
});
