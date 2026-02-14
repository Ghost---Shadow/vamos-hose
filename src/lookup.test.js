const { lookupNmrShifts } = require('./lookup');

jest.mock('./smiles-to-hose', () => ({
  smilesToHoseCodes: jest.fn(),
}));

jest.mock('./database', () => ({
  loadDatabase: jest.fn(),
  queryHose: jest.fn(),
}));

const { smilesToHoseCodes } = require('./smiles-to-hose');
const { loadDatabase, queryHose } = require('./database');

describe('lookupNmrShifts', () => {
  beforeEach(() => {
    jest.clearAllMocks();
  });

  test('wires smilesToHoseCodes -> queryHose -> result', () => {
    smilesToHoseCodes.mockReturnValue([
      { atom: 'C', index: 0, hose: 'HOSE_A' },
      { atom: 'C', index: 1, hose: 'HOSE_B' },
    ]);

    const mockDb = {};
    loadDatabase.mockReturnValue(mockDb);

    queryHose
      .mockReturnValueOnce({ avgShift: 25.3, smiles: 'CC', nucleus: 'C', solvents: {} })
      .mockReturnValueOnce({ avgShift: 130.5, smiles: 'C=C', nucleus: 'C', solvents: {} });

    const result = lookupNmrShifts('CC', { nucleus: '13C' });

    expect(smilesToHoseCodes).toHaveBeenCalledWith('CC', { nucleus: '13C' });
    expect(loadDatabase).toHaveBeenCalled();
    expect(queryHose).toHaveBeenCalledWith(mockDb, 'HOSE_A');
    expect(queryHose).toHaveBeenCalledWith(mockDb, 'HOSE_B');
    expect(result).toEqual([
      { shift: 25.3, atom: 'C', hose: 'HOSE_A', smiles: 'CC' },
      { shift: 130.5, atom: 'C', hose: 'HOSE_B', smiles: 'C=C' },
    ]);
  });

  test('skips atoms with no database hit', () => {
    smilesToHoseCodes.mockReturnValue([
      { atom: 'C', index: 0, hose: 'KNOWN' },
      { atom: 'C', index: 1, hose: 'UNKNOWN' },
      { atom: 'C', index: 2, hose: 'KNOWN2' },
    ]);

    loadDatabase.mockReturnValue({});

    queryHose
      .mockReturnValueOnce({ avgShift: 40.0, smiles: 'C', nucleus: 'C', solvents: {} })
      .mockReturnValueOnce(null)
      .mockReturnValueOnce({ avgShift: 77.0, smiles: 'CO', nucleus: 'C', solvents: {} });

    const result = lookupNmrShifts('CCC');

    expect(result.length).toBe(2);
    expect(result[0].shift).toBe(40.0);
    expect(result[1].shift).toBe(77.0);
  });

  test('returns empty array when no HOSE codes match', () => {
    smilesToHoseCodes.mockReturnValue([
      { atom: 'C', index: 0, hose: 'X' },
    ]);

    loadDatabase.mockReturnValue({});
    queryHose.mockReturnValue(null);

    const result = lookupNmrShifts('C');

    expect(result).toEqual([]);
  });

  test('defaults nucleus to 13C', () => {
    smilesToHoseCodes.mockReturnValue([]);
    loadDatabase.mockReturnValue({});

    lookupNmrShifts('C');

    expect(smilesToHoseCodes).toHaveBeenCalledWith('C', { nucleus: '13C' });
  });

  test('passes custom nucleus through', () => {
    smilesToHoseCodes.mockReturnValue([]);
    loadDatabase.mockReturnValue({});

    lookupNmrShifts('C', { nucleus: '1H' });

    expect(smilesToHoseCodes).toHaveBeenCalledWith('C', { nucleus: '1H' });
  });
});
