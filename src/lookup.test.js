import { describe, test, expect } from 'bun:test';
import { lookupNmrShifts } from './lookup.js';

describe('lookupNmrShifts', () => {
  test('returns array of shift predictions for methane', () => {
    const result = lookupNmrShifts('C', { nucleus: '13C' });
    expect(Array.isArray(result)).toBe(true);
    expect(result.length).toBeGreaterThan(0);
    if (result.length > 0) {
      expect(result[0]).toHaveProperty('shift');
      expect(result[0]).toHaveProperty('atom');
      expect(result[0]).toHaveProperty('hose');
      expect(result[0]).toHaveProperty('smiles');
    }
  });

  test('returns array of shift predictions for ethane', () => {
    const result = lookupNmrShifts('CC', { nucleus: '13C' });
    expect(Array.isArray(result)).toBe(true);
  });

  test('defaults nucleus to 13C', () => {
    const result = lookupNmrShifts('C');
    expect(Array.isArray(result)).toBe(true);
  });

  test('handles different nuclei', () => {
    const result = lookupNmrShifts('C', { nucleus: '1H' });
    expect(Array.isArray(result)).toBe(true);
  });
});
