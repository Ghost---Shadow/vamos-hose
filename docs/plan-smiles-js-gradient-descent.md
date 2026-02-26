# Gradient Descent Engine: 13C NMR → Molecular Identification

## Goal

Given a set of experimental 13C NMR chemical shifts (peak positions in PPM), find the molecular structure that best explains them. This is an inverse problem: NMR spectra → molecule.

## Algorithm (v2 — Batch-Attach with Formula Filtering)

Peak-driven batch attachment instead of single-atom gradient descent:

```
1. Start with an initial molecule (corrupted/partial structure)
2. Predict its 13C NMR spectrum via HOSE code lookup
3. Identify MISSING peaks (target peaks with no close predicted match)
4. For each missing peak:
   a. Reverse lookup: shift → HOSE codes (via inverted index)
   b. Filter candidates by molecular formula budget (target formula - current formula)
   c. Extract fragment SMILES from matching HOSE environments
5. Attach ALL fragments in one go (one per missing peak)
6. Predict spectrum of the combined molecule
7. Identify SHIFTED peaks — peaks we thought were covered but moved due to
   the combined molecular environment changing their chemical shifts
8. Repeat from step 3 for shifted peaks until converged
```

**Key insight**: The old algorithm added one atom/fragment per step, requiring 20+ steps
for a molecule missing 5 atoms. Batch-attach covers all gaps in a single step, then
fixes secondary shifts caused by the new environment.

**Molecular formula filtering**: Given the target formula (e.g. C33H30N4O2), we know
the exact atom budget remaining. When looking up HOSE codes for a missing peak, we
reject any fragment that would exceed the budget. This dramatically prunes the search
space — instead of trying C, N, O, S, Cl fragments, we only try atoms that are
actually missing from the current structure.

### Fallback: Greedy Gradient Descent (v1)

If batch-attach fails to converge (loss not decreasing), fall back to per-step mutations:

```
1. Generate candidate mutations (substitutions, additions, bond changes, removals, fragments)
2. Evaluate each candidate's predicted spectrum
3. Accept the candidate with lowest loss
4. Repeat until converged or max steps reached
```

Escape mechanisms for local minima:
- **Tabu list**: visited SMILES are never revisited
- **Perturbation**: accept slightly worse candidates when stuck on a plateau (1.2x, 1.5x, 2.0x tolerance)
- **Backtracking**: rewind to earlier trajectory points (60%, 30%, 10% of history)

## HOSE-Guided Fragment Selection

Instead of blind enumeration, the engine uses the nmrshiftdb2 database in reverse to guide mutations:

```
poorly matched target peaks
  → reverse lookup: shift → HOSE codes (via inverted index, 250 PPM files)
  → HOSE → SMILES (via queryHose)
  → extract parenthesized branches as fragment SMILES
  → attach fragments to current molecule
```

This makes mutations targeted — if the spectrum is missing a peak around 170 PPM (carbonyl region), the engine looks up what structural environments produce 170 PPM shifts and tries attaching those fragments.

## String-Level Mutations (current approach)

All mutations work by direct SMILES string manipulation. This is critical — earlier AST-based approaches failed on fused rings (naphthalene) and nested attachments.

**How it works:**
1. Tokenize the SMILES string with regex: `/\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnosp]/g`
2. This finds every atom token regardless of molecule topology — ring atoms, chain atoms, atoms inside branches, fused ring atoms
3. For each atom position, try mutations (substitute, insert branch, attach fragment)
4. Validate each result with `normalize()` — invalid SMILES are silently skipped
5. Deduplicate by canonical SMILES

**Mutation types:**
| Type | Method | Example |
|---|---|---|
| Substitution | Replace atom token | `CCO` → `CCN` (O→N) |
| Addition | Insert `(X)` after atom | `CC` → `C(O)C` |
| Bond change | Replace/insert `=`/`#` | `CC` → `C=C` |
| Removal | Strip `(...)` branches | `CC(O)C` → `CCC` |
| Fragment attach | Insert `(frag)` after atom | `CC` → `C(C(=O)O)C` |

## Pipeline

```
src/identify.js          — gradient descent loop, loss computation, HOSE guidance
src/mutate.js            — SMILES string-level mutation enumeration
src/estimate.js          — reverse lookup: shift → HOSE → SMILES
src/smiles-to-hose.js    — forward prediction: SMILES → HOSE codes
src/hose-generator.js    — CDK two-pass HOSE algorithm
src/database.js          — chunked nmrshiftdb2 lookup (210MB, 1.4M entries)
src/lookup.js            — HOSE code → shift lookup
```

**Forward path** (predict spectrum from molecule):
```
SMILES → OpenChemLib parse → HOSE codes → nmrshiftdb2 lookup → predicted shifts
```

**Reverse path** (guide mutations from target spectrum):
```
target shift → inverted index PPM file → HOSE codes → queryHose → SMILES → extract fragments
```

## Test Files

- `src/identify.test.js` — 34 unit + integration tests (computeLoss, identifyWorstAtoms, mutate.js, predictShiftsWithAtomIndices, identifyMolecule)
- `src/identify-compounds.test.js` — 8 common compounds diagnostic with 5% loss threshold

## Current Status

**7/8 compounds within 5% loss threshold:**

| Compound | Target Carbons | Initial Loss | Final Loss | 5% Threshold | Steps | Status |
|---|---|---|---|---|---|---|
| Aspirin | 9 | 283.3 | 38.7 | 58.0 | 20 | PASS |
| Caffeine | 8 | 579.7 | 22.0 | 40.1 | 25 | PASS |
| Ibuprofen | 13 | 934.9 | 26.3 | 58.2 | 22 | PASS |
| Acetaminophen | 8 | 245.1 | 33.2 | 48.6 | 25 | PASS |
| Naproxen | 14 | 513.5 | 46.6 | 80.7 | 25 | PASS |
| Toluene | 7 | 165.6 | 0.0 | 40.8 | 1 | PASS |
| Ethanol | 2 | 70.7 | 7.0 | 3.9 | 22 | FAIL |
| Acetic acid | 2 | 72.9 | 0.0 | 9.9 | 1 | PASS |

## Remaining Issue: Ethanol

Ethanol (`CCO`, 2 carbons, target shifts [18.6, 58.4]) has a tight 5% threshold of 3.9 PPM. The engine starts from `C` and quickly finds `C(NC)` (loss=21.9) but then gets stuck trying exotic atoms (Cl, S, aromatic c) instead of the simpler `CCO`. The HOSE reverse lookup suggests fragments containing N and Cl for peaks near 18.6 and 58.4 PPM because those shifts appear in many molecular environments.

Possible fixes:
- Prioritize common organic atoms (C, N, O) over halogens in mutation ordering
- Add a "simplicity bias" that penalizes molecules with unusual atom combinations
- Increase the perturbation tolerance to escape the `C(NC)` basin
