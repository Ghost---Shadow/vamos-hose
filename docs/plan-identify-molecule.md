# Plan: Identify Candidate Molecules from 13C NMR Peaks

## Background / Literature

Kwon et al. (2021, Scientific Reports) describe three approaches for NMR→molecule search:
1. **Chemical shift similarity** — predict shifts for candidates, compare to observed peaks. Fast.
2. **Spectral similarity** — compare raw spectra. Requires raw data we don't have.
3. **Molecule-to-spectrum estimation** — L-BFGS optimization. Too heavy for client-side JS.

We implement approach #1 using our existing HOSE-code forward prediction pipeline
(`lookupNmrShifts`), which uses the same NMRShiftDB2 database as Kwon's HOSE baseline.

## Algorithm: Two-Phase Candidate Identification

### Phase 1 — Candidate Collection (inverted index scan)

For each observed peak (within tolerance), query the inverted index to get top-K
HOSE code matches. For each HOSE code, resolve the source SMILES from the chunk
database via `queryHose()`. Build a frequency map:

```
candidateMap: Map<SMILES, {
  peaksMatched: Set<number>,   // which observed peaks this molecule explains
  totalError: number,          // sum of |dbShift - observedPeak|
  hoses: string[],             // HOSE codes that matched
}>
```

A molecule that appears as the source for HOSE codes matching MANY different peaks
is a strong candidate. Filter to candidates matching ≥ 30% of observed peaks.

**Parameters**: K = 20 matches per peak (wider net than the current top-1).

### Phase 2 — Forward Scoring (precise re-ranking)

For each candidate SMILES (up to top-M by frequency from Phase 1):

1. Run `lookupNmrShifts(smiles)` → array of predicted shifts for ALL carbons
2. Score by greedy bipartite matching of predicted shifts vs observed peaks:
   - Sort predicted shifts ascending
   - Sort observed peaks ascending
   - Greedy assign: for each predicted shift, find nearest unassigned observed
     peak within tolerance. Record matched pairs.
3. Compute composite score:
   ```
   coverage  = matchedPeaks / observedPeakCount          (0–1)
   accuracy  = 1 / (1 + meanAbsError)                    (0–1)
   sizePenalty = 1 - |predictedCount - observedCount| / max(predictedCount, observedCount)
   score = coverage * accuracy * sizePenalty
   ```
4. Carbon count filter: skip candidates where |carbonCount - peakCount| > peakCount
   (molecule has more than 2x or less than 0.5x the expected carbons)

### Phase 3 — Return top 5

Sort by score descending. Return:
```
[{
  rank: 1,
  smiles: "CC(=O)Oc1ccccc1C(=O)O",
  score: 0.85,
  predictedShifts: [21.1, 123.4, 124.8, ...],
  matchedPeaks: 9,
  totalPeaks: 9,
  meanError: 1.2,
  carbonCount: 9,
}]
```

## API Design

```js
// src/identify.js
export async function identifyCandidates(options) {
  // options: { peaks: number[], tolerance?: number, maxCandidates?: number,
  //            topKPerPeak?: number, onProgress?: (phase, pct) => void }
  // returns: Promise<Array<CandidateResult>>
}
```

## File Changes

| File | Change | Delegate to |
|------|--------|-------------|
| `src/identify.js` | New file — core algorithm | **Opus** (algorithmic) |
| `src/identify.test.js` | New file — unit + integration tests | **Sonnet** |
| `index.js` | Add export for `identifyCandidates` | **Haiku** |

## Cost Analysis

For aspirin (9 peaks, tolerance ±2):
- Phase 1: loads ~45 PPM files (~9 MB), resolves ~180 HOSE codes → ~30-60 unique SMILES
- Phase 2: runs `lookupNmrShifts` on top ~30 candidates → loads ~50-100 chunks (~40-80 MB)
- Total: ~50-90 MB network, ~10-20s on fast connection
- Returns 5 ranked candidate molecules with scores

## Edge Cases

- **Fewer SMILES than expected**: If inverted index yields < 5 unique SMILES, return what we have
- **No forward prediction hits**: If `lookupNmrShifts` returns empty for a candidate, score = 0
- **Duplicate SMILES**: Canonical SMILES deduplication via string equality
- **Very large molecules** (30+ carbons): Carbon count filter prevents wasting time on poor matches
