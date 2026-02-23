# Plan: Molecular Gradient Descent Engine for NMR Spectrum Matching

## Context
User wants to iteratively optimize a molecule to match a target 13C NMR spectrum. Given observed peaks + an initial guess, the engine applies discrete molecular mutations (atom swaps, bond changes, additions, removals), evaluates each via forward HOSE-based prediction, and follows the steepest descent direction. The trajectory of N steps is recorded for visualization.

Inspired by Kwon et al. (2021, Scientific Reports) — "chemical shift similarity" approach using HOSE codes + NMRShiftDB2, which is exactly our existing prediction pipeline.

## Algorithm: Steepest Descent on Molecular Graphs

```
for step = 1..maxSteps:
  1. Predict shifts for current molecule via HOSE codes → queryHose
  2. Compute loss = greedy_match_error(predicted, target) + penalty * |unmatched|
  3. Attribute error to specific atoms (worst-contributing carbons)
  4. Enumerate focused mutations on top-K worst atoms only (~20-30 candidates)
  5. Evaluate each mutation: forward predict → compute loss
  6. Accept mutation with lowest loss (if improved)
  7. Record {step, smiles, loss, mutation} in trajectory
```

### Mutation Operators (discrete "gradient" directions)
1. **substituteAtom** — Change element: C↔N, C↔O, C↔S, C↔F, C↔Cl, C↔Br
2. **changeBondOrder** — Toggle: single↔double, single↔triple
3. **addAtom** — Attach C/N/O to atom with free valence
4. **removeAtom** — Remove terminal (degree-1, non-ring) heavy atoms

### Loss Function
Greedy nearest-neighbor matching of predicted vs target shifts:
- `loss = Σ|matched errors| + 50 * (unmatchedTarget + unmatchedPred)`
- Per-atom error attribution → focus mutations on worst atoms

## New Files (3 new, 1 modified)

### 1. `src/mutate.js` — Mutation operators [Sonnet]
- `cloneMolecule(smiles)` — round-trip via OCL.Molecule.fromSmiles + toSmiles
- `tryMutation(smiles, mutateFn)` — apply mutation, validate, return new SMILES or null
- `enumerateSubstitutions(smiles, atomIndices)` — all element swaps for given atoms
- `enumerateBondChanges(smiles, atomIndices)` — bond order changes for adjacent bonds
- `enumerateAdditions(smiles, atomIndices)` — add C/N/O to atoms with free valence
- `enumerateRemovals(smiles, atomIndices)` — remove terminal atoms
- `enumerateAllMutations(smiles, atomIndices)` — all 4 types, deduplicated by SMILES
- Uses: `import OCL from 'openchemlib'` (v9.20.0, full build, already in deps)

### 2. `src/identify.js` — Gradient descent engine [Opus]
- `computeLoss(predictedShifts, targetShifts, opts)` → `{ loss, assignments, predAtomError }`
- `predictShiftsWithAtomIndices(smiles)` → `[{ shift, atomIndex, hose }]`
  - Calls `smilesToHoseCodes()` for atom indices + `queryHose()` with truncation fallback
  - Duplicates truncation logic from `lookup.js` to get atom index tracking
- `identifyWorstAtoms(lossResult, prediction, topK)` → atom indices with highest error
- `identifyMolecule(targetShifts, options)` → `{ smiles, loss, steps, trajectory, predictedShifts }`
  - Options: `{ startSmiles, maxSteps, topK, unmatchedPenalty, onStep }`
  - Returns full trajectory: `[{ step, smiles, loss, mutation, predictedShifts }]`
- Uses: `smilesToHoseCodes` from `./smiles-to-hose.js`, `queryHose`/`preloadChunks` from `./database.js`

### 3. `src/identify.test.js` — Tests [Sonnet]
**Unit tests:**
- `computeLoss`: perfect match=0, matched errors summed, unmatched penalty, greedy assignment
- `enumerateSubstitutions/BondChanges/Additions/Removals`: valid outputs, dedup, no identity
- `tryMutation`: valid mutation returns SMILES, invalid returns null

**Integration tests (60s timeout):**
- Round-trip aspirin: predict aspirin shifts → target, start from salicylic acid → loss decreases
- Benzene→toluene: target=toluene shifts, initial=benzene → loss decreases
- Self-identification: start from ethanol with ethanol target → loss ≈ 0
- Trajectory structure: step 0 = initial, has smiles/loss/mutation fields

### 4. `index.js` — Add export [Haiku]
```js
export { identifyMolecule, computeLoss } from './src/identify.js';
```

## Key Reused Functions (no changes needed)
- `smilesToHoseCodes(smiles)` → `[{atom, index, hose}]` — `src/smiles-to-hose.js`
- `queryHose(hoseCode)` → `{avgShift, smiles, nucleus}` — `src/database.js`
- `preloadChunks(hoseCodes)` — `src/database.js`
- `computeWeightedAvg(entry)` — `src/database.js`

## Verification
1. `bun test` — all existing 358 tests still pass + new identify tests pass
2. Integration test: aspirin round-trip shows loss decreasing across trajectory
3. Trajectory output has correct structure for future UI visualization
