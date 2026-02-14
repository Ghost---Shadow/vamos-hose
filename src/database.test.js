const {
  loadDatabase,
  queryHose,
  computeWeightedAvg,
  extractSolvents,
} = require('./database');

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
  const mockDb = {
    'HCC(': {
      n: 'C',
      s: 'CCC',
      'C_CDCl3': { min: 25.0, max: 27.0, avg: 26.0, cnt: 10 },
      'C_DMSO': { min: 26.0, max: 28.0, avg: 27.0, cnt: 5 },
    },
  };

  test('returns shift data for known HOSE code', () => {
    const result = queryHose(mockDb, 'HCC(');
    expect(result).not.toBeNull();
    expect(result.smiles).toBe('CCC');
    expect(result.nucleus).toBe('C');
    // weighted avg: (26*10 + 27*5) / 15 = (260+135)/15 = 26.333... -> 26.3
    expect(result.avgShift).toBe(26.3);
    expect(result.solvents).toEqual({
      'C_CDCl3': { min: 25.0, max: 27.0, avg: 26.0, cnt: 10 },
      'C_DMSO': { min: 26.0, max: 28.0, avg: 27.0, cnt: 5 },
    });
  });

  test('returns null for unknown HOSE code', () => {
    expect(queryHose(mockDb, 'NONEXISTENT')).toBeNull();
  });
});

describe('loadDatabase', () => {
  test('returns an object with HOSE code keys', () => {
    const db = loadDatabase();
    expect(typeof db).toBe('object');
    expect(Object.keys(db).length).toBeGreaterThan(0);
  });

  test('entries have n and s fields', () => {
    const db = loadDatabase();
    const firstKey = Object.keys(db)[0];
    expect(db[firstKey]).toHaveProperty('n');
    expect(db[firstKey]).toHaveProperty('s');
  });

  test('solvent entries have min, max, avg, cnt', () => {
    const db = loadDatabase();
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

  test('returns same reference on second call (cached)', () => {
    const db1 = loadDatabase();
    const db2 = loadDatabase();
    expect(db1).toBe(db2);
  });
});
