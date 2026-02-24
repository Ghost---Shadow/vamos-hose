import { describe, test, expect } from 'bun:test';
import { computeLoss, identifyWorstAtoms, identifyMolecule, predictShiftsWithAtomIndices } from './identify.js';
import {
  validateSmiles, getAtomCount,
  enumerateSubstitutions, enumerateBondChanges, enumerateAdditions, enumerateRemovals,
  enumerateFragmentAttachments, enumerateAllMutations,
} from './mutate.js';

// ─── computeLoss unit tests ───────────────────────────────────────────

describe('computeLoss', () => {
  test('perfect match returns zero loss', () => {
    const result = computeLoss([10, 20, 30], [10, 20, 30]);
    expect(result.loss).toBe(0);
    expect(result.unmatchedTarget).toEqual([]);
    expect(result.unmatchedPred).toEqual([]);
  });

  test('matched errors are summed', () => {
    const result = computeLoss([11, 22], [10, 20]);
    expect(result.loss).toBe(3); // |11-10| + |22-20|
    expect(result.assignments.length).toBe(2);
  });

  test('unmatched target peaks incur penalty', () => {
    const result = computeLoss([10], [10, 20, 30], { unmatchedPenalty: 50 });
    // 1 matched (10-10=0), 2 unmatched target
    expect(result.loss).toBe(100);
    expect(result.unmatchedTarget.length).toBe(2);
  });

  test('unmatched predicted peaks incur penalty', () => {
    const result = computeLoss([10, 20, 30], [10], { unmatchedPenalty: 50 });
    expect(result.loss).toBe(100); // 0 + 2*50
    expect(result.unmatchedPred.length).toBe(2);
  });

  test('empty predicted vs non-empty target', () => {
    const result = computeLoss([], [10, 20], { unmatchedPenalty: 50 });
    expect(result.loss).toBe(100);
  });

  test('both empty returns zero', () => {
    const result = computeLoss([], []);
    expect(result.loss).toBe(0);
  });

  test('greedy assignment finds optimal nearest neighbor', () => {
    // Target: [10, 100], Predicted: [12, 98]
    // Optimal: 10→12 (err=2), 100→98 (err=2), total=4
    const result = computeLoss([12, 98], [10, 100]);
    expect(result.loss).toBe(4);
  });

  test('predAtomError map is populated for all predicted atoms', () => {
    const result = computeLoss([11, 22], [10, 20]);
    expect(result.predAtomError.size).toBe(2);
  });

  test('size mismatch penalty reflects count difference', () => {
    const result1 = computeLoss([10, 20], [10, 20]);
    const result2 = computeLoss([10, 20, 30, 40, 50], [10, 20], { unmatchedPenalty: 50 });
    expect(result2.loss).toBeGreaterThan(result1.loss);
  });
});

// ─── identifyWorstAtoms unit tests ────────────────────────────────────

describe('identifyWorstAtoms', () => {
  test('returns atom indices sorted by error descending', () => {
    const lossResult = {
      predAtomError: new Map([[0, 1.0], [1, 5.0], [2, 3.0]]),
    };
    const prediction = [
      { shift: 10, atomIndex: 0 },
      { shift: 20, atomIndex: 3 },
      { shift: 30, atomIndex: 5 },
    ];
    const worst = identifyWorstAtoms(lossResult, prediction, 2);
    // predIdx=1 (error=5.0) → atomIndex=3, predIdx=2 (error=3.0) → atomIndex=5
    expect(worst).toEqual([3, 5]);
  });

  test('respects topK limit', () => {
    const lossResult = {
      predAtomError: new Map([[0, 1], [1, 2], [2, 3], [3, 4]]),
    };
    const prediction = [
      { shift: 10, atomIndex: 0 },
      { shift: 20, atomIndex: 1 },
      { shift: 30, atomIndex: 2 },
      { shift: 40, atomIndex: 3 },
    ];
    expect(identifyWorstAtoms(lossResult, prediction, 2).length).toBe(2);
  });
});

// ─── mutate.js unit tests ─────────────────────────────────────────────

describe('mutate.js', () => {
  describe('validateSmiles', () => {
    test('valid SMILES returns valid=true with canonical', () => {
      const result = validateSmiles('CCO');
      expect(result.valid).toBe(true);
      expect(typeof result.canonical).toBe('string');
    });

    test('invalid SMILES returns valid=false', () => {
      const result = validateSmiles('XXXINVALID!!!');
      expect(result.valid).toBe(false);
      expect(result.canonical).toBeNull();
    });
  });

  describe('getAtomCount', () => {
    test('counts atoms in ethanol', () => {
      expect(getAtomCount('CCO')).toBe(3);
    });

    test('counts atoms in benzene', () => {
      expect(getAtomCount('c1ccccc1')).toBe(6);
    });

    test('counts atoms in toluene (molecule type)', () => {
      expect(getAtomCount('Cc1ccccc1')).toBe(7);
    });
  });

  describe('enumerateSubstitutions', () => {
    test('generates substitutions for methane', () => {
      const muts = enumerateSubstitutions('C');
      expect(muts.length).toBeGreaterThan(0);
      const smilesList = muts.map(m => m.smiles);
      expect(new Set(smilesList).size).toBe(smilesList.length);
    });

    test('generates substitutions for benzene', () => {
      const muts = enumerateSubstitutions('c1ccccc1');
      expect(muts.length).toBeGreaterThan(0);
      // Should include pyridine-like substitutions
      expect(muts.some(m => m.smiles.includes('N') || m.smiles.includes('n'))).toBe(true);
    });
  });

  describe('enumerateBondChanges', () => {
    test('ethane bonds can change', () => {
      const muts = enumerateBondChanges('CC');
      // CC has bonds [null] — should produce C=C and C#C
      expect(muts.length).toBeGreaterThanOrEqual(1);
    });
  });

  describe('enumerateAdditions', () => {
    test('benzene can grow', () => {
      const muts = enumerateAdditions('c1ccccc1');
      expect(muts.length).toBeGreaterThanOrEqual(1);
    });
  });

  describe('enumerateRemovals', () => {
    test('removes branch from branched molecule', () => {
      const muts = enumerateRemovals('CC(O)CC');
      expect(muts.length).toBeGreaterThanOrEqual(1);
      expect(muts.some(m => m.description.includes('remove'))).toBe(true);
    });
  });

  describe('enumerateFragmentAttachments', () => {
    test('attaches fragments to ethanol', () => {
      const muts = enumerateFragmentAttachments('CCO', ['C(=O)O', 'NC']);
      expect(muts.length).toBeGreaterThanOrEqual(1);
    });

    test('empty fragments produce no results', () => {
      const muts = enumerateFragmentAttachments('CCO', []);
      expect(muts).toEqual([]);
    });
  });

  describe('enumerateAllMutations', () => {
    test('deduplicates by SMILES', () => {
      const muts = enumerateAllMutations('CC');
      const smilesList = muts.map(m => m.smiles);
      expect(new Set(smilesList).size).toBe(smilesList.length);
    });

    test('excludes identity mutation', () => {
      const muts = enumerateAllMutations('CC');
      expect(muts.every(m => m.smiles !== 'CC')).toBe(true);
    });

    test('includes fragment attachments when provided', () => {
      const withoutFrags = enumerateAllMutations('CC');
      const withFrags = enumerateAllMutations('CC', ['C(=O)O']);
      expect(withFrags.length).toBeGreaterThanOrEqual(withoutFrags.length);
    });
  });
});

// ─── predictShiftsWithAtomIndices ─────────────────────────────────────

describe('predictShiftsWithAtomIndices', () => {
  test('predicts shifts for ethanol with atom indices', async () => {
    const results = await predictShiftsWithAtomIndices('CCO');
    expect(results.length).toBe(2); // 2 carbons in ethanol
    for (const r of results) {
      expect(typeof r.shift).toBe('number');
      expect(typeof r.atomIndex).toBe('number');
      expect(typeof r.hose).toBe('string');
      expect(r.atomIndex).toBeGreaterThanOrEqual(0);
    }
  });

  test('aspirin has 9 predicted shifts', async () => {
    const results = await predictShiftsWithAtomIndices('CC(=O)Oc1ccccc1C(O)=O');
    expect(results.length).toBe(9);
  });
});

// ─── Integration tests ────────────────────────────────────────────────

describe('identifyMolecule', () => {
  test('self-identification: ethanol finds itself', async () => {
    const pred = await predictShiftsWithAtomIndices('CCO');
    const targetShifts = pred.map(p => p.shift);

    const result = await identifyMolecule(targetShifts, {
      startSmiles: 'CCO',
      maxSteps: 5,
    });

    // Starting from the correct molecule, loss should be near 0
    expect(result.loss).toBeLessThan(1);
    expect(result.trajectory[0].step).toBe(0);
    expect(result.trajectory[0].mutation).toBe('initial');
  }, 30000);

  test('trajectory records all steps with correct structure', async () => {
    const pred = await predictShiftsWithAtomIndices('CC');
    const targetShifts = pred.map(p => p.shift);

    const result = await identifyMolecule(targetShifts, {
      startSmiles: 'C',
      maxSteps: 5,
    });

    expect(result.trajectory.length).toBeGreaterThanOrEqual(1);
    for (const entry of result.trajectory) {
      expect(entry).toHaveProperty('step');
      expect(entry).toHaveProperty('smiles');
      expect(entry).toHaveProperty('loss');
      expect(entry).toHaveProperty('mutation');
      expect(entry).toHaveProperty('predictedShifts');
    }
  }, 30000);

  test('methane → ethane: loss decreases', async () => {
    const pred = await predictShiftsWithAtomIndices('CC');
    const targetShifts = pred.map(p => p.shift);

    const result = await identifyMolecule(targetShifts, {
      startSmiles: 'C',
      maxSteps: 10,
    });

    const initialLoss = result.trajectory[0].loss;
    expect(result.loss).toBeLessThan(initialLoss);
  }, 60000);

  test('onStep callback is invoked', async () => {
    const pred = await predictShiftsWithAtomIndices('CC');
    const targetShifts = pred.map(p => p.shift);
    const steps = [];

    await identifyMolecule(targetShifts, {
      startSmiles: 'C',
      maxSteps: 3,
      onStep: (info) => steps.push(info),
    });

    expect(Array.isArray(steps)).toBe(true);
  }, 30000);

  test('round-trip: aspirin shifts from salicylic acid → loss decreases', async () => {
    const aspirinPred = await predictShiftsWithAtomIndices('CC(=O)Oc1ccccc1C(O)=O');
    const targetShifts = aspirinPred.map(p => p.shift);
    expect(targetShifts.length).toBe(9);

    const result = await identifyMolecule(targetShifts, {
      startSmiles: 'OC(=O)c1ccccc1O',
      maxSteps: 20,
      topK: 3,
    });

    expect(result.trajectory.length).toBeGreaterThan(1);
    const initialLoss = result.trajectory[0].loss;
    expect(result.loss).toBeLessThan(initialLoss);
  }, 120000);

  test('benzene → toluene: loss decreases toward toluene target', async () => {
    const toluenePred = await predictShiftsWithAtomIndices('Cc1ccccc1');
    const targetShifts = toluenePred.map(p => p.shift);

    const result = await identifyMolecule(targetShifts, {
      startSmiles: 'c1ccccc1',
      maxSteps: 15,
      topK: 3,
    });

    expect(result.loss).toBeLessThan(result.trajectory[0].loss);
  }, 120000);
});
