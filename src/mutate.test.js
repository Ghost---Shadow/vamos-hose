import { describe, test, expect } from 'bun:test';
import {
  countCarbons, containsOnlyAtoms, getAtomCount, getMolecularWeight, validateSmiles,
  enumerateSubstitutions, enumerateAdditions, enumerateBondChanges,
  enumerateRemovals, enumerateInlineInsertions, enumerateRingClosures,
  enumerateFragmentAttachments, enumerateAllMutations,
} from './mutate.js';

// ── SMILES analysis utilities ──────────────────────────────────────────────────

describe('countCarbons', () => {
  test('counts uppercase C (not Cl)', () => {
    expect(countCarbons('CCO')).toBe(2);
    expect(countCarbons('ClCCl')).toBe(1);
    expect(countCarbons('CCCCl')).toBe(3);
  });

  test('counts lowercase c (aromatic)', () => {
    expect(countCarbons('c1ccccc1')).toBe(6);
  });

  test('counts mixed C and c', () => {
    expect(countCarbons('c1ccc(CC)cc1')).toBe(8);
  });

  test('returns 0 for no carbons', () => {
    expect(countCarbons('O')).toBe(0);
    expect(countCarbons('[NH3]')).toBe(0);
  });
});

describe('containsOnlyAtoms', () => {
  const allowed = new Set(['C', 'c', 'N', 'n', 'O', 'o']);

  test('passes for allowed atoms', () => {
    expect(containsOnlyAtoms('CCNCO', allowed)).toBe(true);
    expect(containsOnlyAtoms('c1ccncc1', allowed)).toBe(true);
  });

  test('fails for disallowed atoms', () => {
    expect(containsOnlyAtoms('CCCl', allowed)).toBe(false);
    expect(containsOnlyAtoms('c1ccc(Br)cc1', allowed)).toBe(false);
    expect(containsOnlyAtoms('CCS', allowed)).toBe(false);
  });
});

describe('getAtomCount', () => {
  test('counts heavy atoms', () => {
    expect(getAtomCount('CCO')).toBe(3);
    expect(getAtomCount('c1ccccc1')).toBe(6);
    expect(getAtomCount('ClCCl')).toBe(3);
  });

  test('handles bracket atoms', () => {
    expect(getAtomCount('[NH3]C')).toBe(2);
  });
});

describe('getMolecularWeight', () => {
  test('returns 0 for empty/no atoms', () => {
    expect(getMolecularWeight('')).toBe(0);
  });

  test('estimates MW for simple molecules', () => {
    // Methane: CH4 = 12.011 + 4*1.008 = 16.043
    const methaneMW = getMolecularWeight('C');
    expect(methaneMW).toBeCloseTo(16.043, 0);
  });

  test('benzene MW is ~78', () => {
    const mw = getMolecularWeight('c1ccccc1');
    expect(mw).toBeGreaterThan(70);
    expect(mw).toBeLessThan(85);
  });
});

describe('validateSmiles', () => {
  test('valid SMILES', () => {
    const result = validateSmiles('CCO');
    expect(result).toEqual({ valid: true, canonical: expect.any(String) });
  });

  test('invalid SMILES', () => {
    const result = validateSmiles('XYZ!!!');
    expect(result).toEqual({ valid: false, canonical: null });
  });
});

// ── Mutation enumerators ───────────────────────────────────────────────────────

describe('enumerateSubstitutions', () => {
  test('generates atom replacements', () => {
    const results = enumerateSubstitutions('CCO');
    expect(results.length).toBeGreaterThan(0);
    // Should contain N→C, O→N, etc.
    const descriptions = results.map(r => r.description);
    expect(descriptions.some(d => d.includes('sub'))).toBe(true);
  });

  test('every result has smiles and description', () => {
    for (const r of enumerateSubstitutions('CCO')) {
      expect(r.smiles).toBeDefined();
      expect(r.description).toBeDefined();
    }
  });

  test('does not include identity (same SMILES)', () => {
    const canonical = validateSmiles('CCO').canonical;
    const results = enumerateSubstitutions('CCO');
    expect(results.every(r => r.smiles !== canonical)).toBe(true);
  });

  test('returns empty for invalid SMILES', () => {
    expect(enumerateSubstitutions('(((')).toEqual([]);
  });
});

describe('enumerateAdditions', () => {
  test('adds branches at atom positions', () => {
    const results = enumerateAdditions('CC');
    expect(results.length).toBeGreaterThan(0);
    const descriptions = results.map(r => r.description);
    expect(descriptions.some(d => d.includes('add (C)'))).toBe(true);
    expect(descriptions.some(d => d.includes('add (N)'))).toBe(true);
    expect(descriptions.some(d => d.includes('add (O)'))).toBe(true);
  });

  test('results have more atoms than input', () => {
    const inputAtoms = getAtomCount('CC');
    for (const r of enumerateAdditions('CC')) {
      expect(getAtomCount(r.smiles)).toBeGreaterThanOrEqual(inputAtoms);
    }
  });

  test('returns empty for invalid SMILES', () => {
    expect(enumerateAdditions('(((')).toEqual([]);
  });
});

describe('enumerateBondChanges', () => {
  test('modifies bonds in molecules with explicit bonds', () => {
    const results = enumerateBondChanges('CC=CC');
    expect(results.length).toBeGreaterThan(0);
    const descriptions = results.map(r => r.description);
    // Should be able to remove the double bond or change to triple
    expect(descriptions.some(d => d.includes('bond'))).toBe(true);
  });

  test('can insert bonds between single-bonded atoms', () => {
    const results = enumerateBondChanges('CCCC');
    // Adjacent atoms connected by single bond → try inserting = or #
    expect(results.length).toBeGreaterThan(0);
  });

  test('returns empty for single atom', () => {
    expect(enumerateBondChanges('C')).toEqual([]);
  });
});

describe('enumerateRemovals', () => {
  test('removes parenthesized branches', () => {
    const results = enumerateRemovals('CC(O)CC');
    expect(results.length).toBeGreaterThan(0);
    const descriptions = results.map(r => r.description);
    expect(descriptions.some(d => d.includes('remove (O)'))).toBe(true);
  });

  test('removes terminal atom', () => {
    const results = enumerateRemovals('CCCC');
    expect(results.some(r => r.description === 'remove last atom')).toBe(true);
  });

  test('result has fewer atoms than input', () => {
    const inputAtoms = getAtomCount('CC(O)CC');
    for (const r of enumerateRemovals('CC(O)CC')) {
      expect(getAtomCount(r.smiles)).toBeLessThan(inputAtoms);
    }
  });
});

describe('enumerateInlineInsertions', () => {
  test('inserts atoms at boundaries', () => {
    const results = enumerateInlineInsertions('CC(O)C');
    expect(results.length).toBeGreaterThan(0);
    const descriptions = results.map(r => r.description);
    expect(descriptions.some(d => d.includes('inline'))).toBe(true);
  });

  test('inserts after closing paren', () => {
    const results = enumerateInlineInsertions('C(O)C');
    const descriptions = results.map(r => r.description);
    // After ')' should produce inline insertion
    expect(descriptions.some(d => d.includes('inline C'))).toBe(true);
  });

  test('appends at end of SMILES', () => {
    const results = enumerateInlineInsertions('CC');
    // Should have candidates with appended atoms
    expect(results.length).toBeGreaterThan(0);
  });
});

describe('enumerateRingClosures', () => {
  test('generates ring closures for long chains', () => {
    const results = enumerateRingClosures('CCCCCC');
    expect(results.length).toBeGreaterThan(0);
    const descriptions = results.map(r => r.description);
    expect(descriptions.some(d => d.includes('ring close'))).toBe(true);
  });

  test('preserves atom count', () => {
    const inputAtoms = getAtomCount('CCCCCC');
    for (const r of enumerateRingClosures('CCCCCC')) {
      expect(getAtomCount(r.smiles)).toBe(inputAtoms);
    }
  });

  test('returns empty for too-short chains', () => {
    // Need at least 4 atoms to form a ring (3 positions apart)
    expect(enumerateRingClosures('CC')).toEqual([]);
  });
});

describe('enumerateFragmentAttachments', () => {
  test('attaches fragments at atom positions', () => {
    const results = enumerateFragmentAttachments('CC', ['O', 'N']);
    expect(results.length).toBeGreaterThan(0);
    const descriptions = results.map(r => r.description);
    expect(descriptions.some(d => d.includes('attach (O)'))).toBe(true);
    expect(descriptions.some(d => d.includes('attach (N)'))).toBe(true);
  });

  test('handles ring fragments', () => {
    const results = enumerateFragmentAttachments('CC', ['c1ccccc1']);
    expect(results.length).toBeGreaterThan(0);
    const descriptions = results.map(r => r.description);
    expect(descriptions.some(d => d.includes('attach (c1ccccc1)'))).toBe(true);
  });

  test('returns empty with no fragments', () => {
    expect(enumerateFragmentAttachments('CC', [])).toEqual([]);
    expect(enumerateFragmentAttachments('CC', null)).toEqual([]);
  });
});

describe('enumerateAllMutations', () => {
  test('combines all mutation types', () => {
    const results = enumerateAllMutations('CCCCC');
    expect(results.length).toBeGreaterThan(0);

    const descriptions = results.map(r => r.description);
    // Should have mutations from multiple types
    expect(descriptions.some(d => d.startsWith('sub '))).toBe(true);
    expect(descriptions.some(d => d.startsWith('add '))).toBe(true);
    expect(descriptions.some(d => d.startsWith('inline '))).toBe(true);
    expect(descriptions.some(d => d.startsWith('remove '))).toBe(true);
  });

  test('accepts fragment argument', () => {
    const withFrags = enumerateAllMutations('CC', ['c1ccccc1']);
    const withoutFrags = enumerateAllMutations('CC');
    expect(withFrags.length).toBeGreaterThan(withoutFrags.length);
  });

  test('no duplicate SMILES in results', () => {
    const results = enumerateAllMutations('CCCC');
    const smilesList = results.map(r => r.smiles);
    expect(new Set(smilesList).size).toBe(smilesList.length);
  });

  test('returns empty for invalid SMILES', () => {
    expect(enumerateAllMutations('(((')).toEqual([]);
  });
});

// ── Edge cases ─────────────────────────────────────────────────────────────────

describe('edge cases', () => {
  test('aromatic ring mutations', () => {
    const benzene = 'c1ccccc1';
    const subs = enumerateSubstitutions(benzene);
    expect(subs.length).toBeGreaterThan(0);

    const adds = enumerateAdditions(benzene);
    expect(adds.length).toBeGreaterThan(0);
  });

  test('telmisartan-like molecule generates diverse mutations', () => {
    const telmisartan = 'CCCC1=NC2=C(C=C(C=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C(=O)O)C5=NC6=CC=CC=C6N5C)C';
    const adds = enumerateAdditions(telmisartan);
    const subs = enumerateSubstitutions(telmisartan);
    const inlines = enumerateInlineInsertions(telmisartan);
    const frags = enumerateFragmentAttachments(telmisartan, ['C', 'O', 'c1ccccc1']);

    // Each type should generate meaningful candidates
    expect(adds.length).toBeGreaterThan(5);
    expect(subs.length).toBeGreaterThan(5);
    expect(inlines.length).toBeGreaterThan(0);
    expect(frags.length).toBeGreaterThan(5);
  });

  test('single-atom SMILES', () => {
    expect(enumerateSubstitutions('C').length).toBeGreaterThan(0);
    expect(enumerateAdditions('C').length).toBeGreaterThan(0);
    expect(enumerateRemovals('C')).toEqual([]);
    expect(enumerateRingClosures('C')).toEqual([]);
  });
});
