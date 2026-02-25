# Telmisartan: Gradient Descent Failure Analysis

## The Molecule

Telmisartan: `CCCC1=NC2=C(C=C(C=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C(=O)O)C5=NC6=CC=CC=C6N5C)C`

- 33 carbon environments
- MW ~515 Da
- Key substructures: two benzimidazoles (fused 5+6 rings with N), biphenyl, carboxylic acid, propyl chain, N-methyl, aryl-methyl

Target shifts (sorted):
```
[15.8, 21.4, 21.8, 30.1, 31.6, 53.3, 110.5, 118.9, 120.6, 129.0, 129.0, 129.4,
 130.7, 130.8, 130.8, 131.2, 135.1, 135.5, 135.5, 135.5, 136.5, 137.3, 137.7,
 138.4, 139.3, 141.6, 142.9, 143.4, 144.3, 147.8, 151.7, 155.0, 174.1]
```

5% loss threshold: 193.3 PPM total.

## Current Behavior (from benzene, 60s timeout)

The engine reaches loss=84.1 in 15 steps — **below threshold** — but the 4.5s test timeout kills it at step ~3-4.

### Observed Trajectory (60s run)

```
Step  0: loss=1947.0  C= 6  MW=  78  benzene
Step  1: loss=1410.4  C=10  MW= 134  attach (C(C)CC)     ← propyl-like chain
Step  2: loss=1041.5  C=14  MW= 188  attach (C=C(C)C)
Steps 3-10:                           attach (cc) x8      ← aromatic buildup
Step 10: loss= 264.5  C=30  MW= 397
Step 11: loss= 221.9  C=31  MW= 441  attach (C(=O)O)     ← carboxylic acid
Step 12: loss= 147.8  C=34  MW= 481  attach (CC=C)
Step 13: loss=  93.7  C=33  MW= 546  sub C->Br           ← HALOGEN HACK
Step 15: loss=  84.1  C=33  MW= 624  sub N->Br           ← HALOGEN HACK
Steps 16-18: perturbation, stuck
Step 19: backtrack to step 11
Steps 20-24: same pattern, stuck again at ~84
Step 25: backtrack to step 7, rebuilds, stuck at ~101
```

### Final molecule produced

```
c1c(C(C(C=C(C(cc(BrC=C)(cc(cc))))C))C(c(c(C(=O)O)c(cc))c(cc))C)ccc(cc(Br))c1
```

This is a polyaromatic carbon skeleton with **two bromine atoms** and **zero nitrogen atoms**. It bears no structural resemblance to telmisartan.

## Root Cause Analysis

### Failure 1: Timeout (primary blocker)

The test allows 4.5s. Each step evaluates ~500 candidate mutations, each requiring HOSE code generation + DB lookup. At ~300ms/step, the engine completes ~15 steps in 4.5s but only ~3-4 steps before the bun test timeout kicks in.

The algorithm **does converge** given enough time (84.1 < 193.3 threshold) but can't finish within the test budget.

### Failure 2: No ring-forming mutations

The mutation operators are:
- Substitution: `C→N`, `C→O`, etc.
- Addition: insert `(X)` branch
- Bond change: single↔double↔triple
- Removal: strip `(...)` branches
- Fragment attach: insert `(frag)` branch

**None of these can create ring closures.** SMILES ring closures require paired digits (e.g., `c1ccccc1`), which no mutation operator produces. This means:

- The engine **cannot build benzimidazole** (fused 5+6 ring with N)
- The engine **cannot build any new rings** beyond what the starting molecule has
- Starting from benzene (1 ring), it can only grow tree-shaped branches off that ring

This is the fundamental structural limitation. Telmisartan has 5 rings (two benzimidazoles = 4 rings + biphenyl contributes 2 more). The engine can never build these.

### Failure 3: Halogen count manipulation

The engine discovers that substituting `C→Br` reduces carbon count (and thus the unmatched peak penalty) without changing the aromatic shift profile much. This is a loss function exploit, not real chemistry:

```
Step 12: 34 carbons, loss=147.8
Step 13: sub C->Br → 33 carbons, loss=93.7   (−54 loss from removing 1 unmatched penalty)
```

The unmatched penalty is 50 per peak. Removing one excess carbon saves 50 in penalty while only adding ~4 in match error. The engine learns to "adjust" carbon count with halogens.

### Failure 4: Fragment quality from reverse lookup

The HOSE-guided fragment system (`getGuidedFragments`) queries the inverted index for shifts in the 110-155 PPM range (aromatic region) and extracts parenthesized branches from the resulting SMILES. For telmisartan's aromatic peaks, this produces:

- `(cc)` — generic aromatic pair
- `(C)`, `(N)`, `(O)` — single atoms
- `(C(=O)O)` — carboxylic acid (good!)

But never produces:
- Benzimidazole fragments (`c1nc2ccccc2[nH]1`)
- Biphenyl fragments (`c1ccc(-c2ccccc2)cc1`)
- Any multi-ring fragment

The `extractFragmentsFromSmiles` function only captures parenthesized branches from a single SMILES, which are short linear fragments, never ring systems.

### Failure 5: No nitrogen in the solution

The final molecule has 0 nitrogen atoms despite telmisartan having 4 nitrogens. The engine has no incentive to introduce nitrogen because:

1. Aromatic N shifts (benzimidazole C2, C8 at 143-155 PPM) are close enough to aromatic C shifts (128-140 PPM) that the loss function doesn't strongly prefer N
2. Substituting aromatic `c→n` often produces invalid SMILES (broken ring aromaticity)
3. The engine can match the 143-155 PPM targets by having more aromatic carbons at slightly wrong positions

### Failure 6: Greedy loss assignment allows drift

The greedy nearest-neighbor loss function (`computeLoss`) can assign predicted peaks to wrong targets when there are many similar-valued peaks. With 20+ peaks in the 128-144 PPM range, small errors accumulate without strong gradient toward the correct structure.

## What Needs to Change

### 1. Ring closure mutations (critical)

Add a mutation operator that creates ring closures:
- Pick two atoms in the current molecule
- Add ring closure digits to connect them
- This enables building 5-membered rings (imidazole), fused rings (benzimidazole), etc.

Implementation: insert matching ring digits at two atom positions in the SMILES string.

### 2. Multi-ring fragment library (critical)

Pre-define common heterocyclic fragments that the engine can attach as whole units:
- Benzimidazole: `c1[nH]c2ccccc12`
- Imidazole: `c1cc[nH]n1`
- Pyridine: `c1ccncc1`
- Indole, pyrrole, etc.

These are too complex to discover through single-atom mutations.

### 3. Heteroatom bonus in loss function

Add a term that rewards matching the expected heteroatom count. If the target spectrum has peaks at 150+ PPM (characteristic of C adjacent to N in benzimidazoles), penalize candidates that lack nitrogen.

### 4. Speed optimization

Each step evaluates ~500 candidates × (HOSE generation + DB lookup). Options:
- **Batch HOSE prediction**: predict all candidates' HOSE codes in one pass
- **Early termination**: skip candidates once we find one better than current best by >10
- **Candidate pruning**: only evaluate a random subset of candidates per step
- **MW pre-filter**: reject candidates before HOSE prediction if MW is way off

### 5. Anti-halogen bias

Penalize molecules containing halogens (Br, Cl, F, I) unless the target spectrum shows characteristic halogen-adjacent shifts. For telmisartan (pure C/N/O/H compound), halogen substitution should be discouraged.

## Test Strategy

`test-integration/telmisartan-corruption.test.js` — progressive corruption tests:

| Test | Corruption | Carbons Removed | Expected Difficulty |
|---|---|---|---|
| 0 | Exact telmisartan | 0 | Trivial (loss=0) |
| 1 | Remove propyl (CCC→C) | 2 | Easy |
| 2 | Remove COOH | 2 | Easy |
| 3 | Remove aryl-CH3 | 1 | Easy |
| 4 | Remove N-CH3 | 1 | Easy |
| 5 | Both methyls removed | 2 | Moderate |
| 6 | Remove 2nd benzimidazole | 8 | Hard (requires ring formation) |
| 7 | N→C swap in benzimidazole | 0 (wrong atoms) | Hard (requires heteroatom recognition) |
| 8 | Just biphenyl-COOH | 20 | Very hard (must build from scratch) |

All tests should reach 0% error — the engine must recover from each corruption back to a molecule whose predicted spectrum exactly matches the target.
