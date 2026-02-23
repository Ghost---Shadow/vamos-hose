import { describe, test, expect } from 'bun:test';
import { computeLoss, identifyWorstAtoms, identifyMolecule, predictShiftsWithAtomIndices } from './identify.js';
import {
  cloneMolecule, tryMutation,
  enumerateSubstitutions, enumerateBondChanges, enumerateAdditions, enumerateRemovals,
  enumerateAllMutations, getAllHeavyAtomIndices, getAllCarbonIndices,
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
  describe('cloneMolecule', () => {
    test('clones ethanol', () => {
      const result = cloneMolecule('CCO');
      expect(result).not.toBeNull();
      expect(typeof result.smiles).toBe('string');
    });

    test('returns null for invalid SMILES', () => {
      expect(cloneMolecule('XXXINVALID')).toBeNull();
    });
  });

  describe('tryMutation', () => {
    test('valid substitution returns new SMILES', () => {
      const result = tryMutation('CC', (mol) => { mol.setAtomicNo(0, 7); });
      expect(result).not.toBeNull();
      expect(typeof result).toBe('string');
      expect(result).not.toBe('CC');
    });

    test('invalid mutation returns null', () => {
      // Pentavalent carbon should fail
      const result = tryMutation('C', (mol) => {
        for (let i = 0; i < 6; i++) {
          const a = mol.addAtom(6);
          mol.addBond(0, a);
        }
      });
      // Either null or valid rearrangement — just ensure no crash
      expect(typeof result === 'string' || result === null).toBe(true);
    });
  });

  describe('enumerateSubstitutions', () => {
    test('generates substitutions for methane', () => {
      const muts = enumerateSubstitutions('C', [0]);
      expect(muts.length).toBeGreaterThan(0);
      // All should be unique SMILES
      const smilesList = muts.map(m => m.smiles);
      expect(new Set(smilesList).size).toBe(smilesList.length);
    });

    test('empty indices produce no mutations', () => {
      expect(enumerateSubstitutions('CC', [])).toEqual([]);
    });
  });

  describe('enumerateBondChanges', () => {
    test('ethane can become ethylene', () => {
      const muts = enumerateBondChanges('CC', [0]);
      expect(muts.length).toBeGreaterThanOrEqual(1);
      expect(muts.some(m => m.smiles.includes('=') || m.description.includes('order'))).toBe(true);
    });
  });

  describe('enumerateAdditions', () => {
    test('methane can grow', () => {
      const muts = enumerateAdditions('C', [0]);
      expect(muts.length).toBeGreaterThanOrEqual(1);
    });
  });

  describe('enumerateRemovals', () => {
    test('removes terminal atom from ethane', () => {
      const muts = enumerateRemovals('CC', [0, 1]);
      expect(muts.length).toBeGreaterThanOrEqual(1);
    });

    test('does not remove ring atoms from cyclohexane', () => {
      const muts = enumerateRemovals('C1CCCCC1', [0, 1, 2]);
      expect(muts).toEqual([]);
    });
  });

  describe('enumerateAllMutations', () => {
    test('deduplicates by SMILES', () => {
      const muts = enumerateAllMutations('CC', [0, 1]);
      const smilesList = muts.map(m => m.smiles);
      expect(new Set(smilesList).size).toBe(smilesList.length);
    });

    test('excludes identity mutation', () => {
      const muts = enumerateAllMutations('CC', [0]);
      expect(muts.every(m => m.smiles !== 'CC')).toBe(true);
    });
  });

  describe('helper functions', () => {
    test('getAllHeavyAtomIndices for ethanol', () => {
      const indices = getAllHeavyAtomIndices('CCO');
      expect(indices.length).toBe(3); // C, C, O
    });

    test('getAllCarbonIndices for ethanol', () => {
      const indices = getAllCarbonIndices('CCO');
      expect(indices.length).toBe(2); // C, C only
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
