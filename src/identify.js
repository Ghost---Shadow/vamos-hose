import { smilesToHoseCodes } from './smiles-to-hose.js';
import { queryHose, preloadChunks } from './database.js';
import {
  enumerateAllMutations, enumerateAdditions, enumerateInlineInsertions,
  enumerateFragmentAttachments, enumerateSubstitutions, enumerateBondChanges,
  enumerateRemovals, enumerateRingClosures,
  getAtomCount, getMolecularWeight, containsOnlyAtoms, countCarbons,
} from './mutate.js';
import { estimateFromSpectra, resolveHoseSmiles } from './estimate.js';
import OCL from 'openchemlib';

function normalize(smiles) {
  return OCL.Molecule.fromSmiles(smiles).toSmiles();
}

// Fragment SMILES strings (pre-computed, no smiles-js dependency needed)
const COMMON_FRAGMENTS = [
  'C', 'CC', 'CCC', 'O', 'N', 'C(=O)O', 'NC',
];

const HETEROCYCLIC_FRAGMENTS = [
  'c1nc[nH]c1',   // imidazole
  'c1cnccc1',     // pyridine
  'c1c[nH]cc1',   // pyrrole
  'c1ccccc1',     // phenyl
  'C(=O)O',       // carboxyl
  'N',            // amino
  'C#N',          // cyano
];

/**
 * Extract useful sub-fragments from a full molecule SMILES.
 * Pulls out parenthesized branches which represent functional groups.
 */
function extractFragmentsFromSmiles(smiles) {
  const fragments = new Set();
  // Simple branches: (O), (N), (CC), etc.
  const simpleRegex = /\(([^()]+)\)/g;
  let match;
  while ((match = simpleRegex.exec(smiles)) !== null) {
    const inner = match[1];
    const cleaned = inner.replace(/^[=#:]/, '');
    if (cleaned.length > 0) fragments.add(cleaned);
  }
  // Nested branches: (C(=O)O), (NC(=O)), etc.
  const nestedRegex = /\(([^()]*\([^()]+\)[^()]*)\)/g;
  while ((match = nestedRegex.exec(smiles)) !== null) {
    const inner = match[1];
    const cleaned = inner.replace(/^[=#:]/, '');
    if (cleaned.length > 0) fragments.add(cleaned);
  }
  return [...fragments];
}

/**
 * Get HOSE-guided fragments for poorly-matched peaks.
 * Uses reverse lookup: target shift → HOSE codes → SMILES → extract fragments.
 */
async function getGuidedFragments(targetShifts, currentLoss, currentPrediction) {
  const peaksToLookup = [];

  // Unmatched target peaks
  for (const idx of currentLoss.unmatchedTarget) {
    peaksToLookup.push(targetShifts[idx]);
  }

  // Worst-matched peaks (error > 2 ppm)
  const sortedAssignments = [...currentLoss.assignments].sort((a, b) => b.error - a.error);
  for (const a of sortedAssignments.slice(0, 5)) {
    if (a.error > 2.0) {
      peaksToLookup.push(targetShifts[a.targetIdx]);
    }
  }

  if (peaksToLookup.length === 0) return [...COMMON_FRAGMENTS, ...HETEROCYCLIC_FRAGMENTS];

  // Dedupe peaks (round to nearest integer for PPM bucket lookup)
  const uniquePeaks = [...new Set(peaksToLookup.map(p => Math.round(p)))];

  try {
    const hoseResults = await estimateFromSpectra({
      peaks: uniquePeaks,
      tolerance: 5.0,
      maxResults: 3,
    });

    const enriched = await resolveHoseSmiles(hoseResults);

    const fragments = new Set([...COMMON_FRAGMENTS, ...HETEROCYCLIC_FRAGMENTS]);
    for (const result of Object.values(enriched)) {
      if (!result || !result.smiles) continue;
      for (const frag of extractFragmentsFromSmiles(result.smiles)) {
        fragments.add(frag);
      }
    }
    return [...fragments];
  } catch {
    return [...COMMON_FRAGMENTS];
  }
}

/**
 * Safe wrapper for predictShiftsWithAtomIndices that returns null on parse errors.
 */
async function safePredictShifts(smiles) {
  try {
    return await predictShiftsWithAtomIndices(smiles);
  } catch {
    return null;
  }
}

// Module-level prediction cache: shared across all identifyMolecule calls
const _globalPredictionCache = new Map();
async function globalCachedPredict(smiles) {
  if (_globalPredictionCache.has(smiles)) return _globalPredictionCache.get(smiles);
  const pred = await safePredictShifts(smiles);
  if (pred) _globalPredictionCache.set(smiles, pred);
  return pred;
}

// Static fragment list: avoids expensive reverse HOSE lookups in getGuidedFragments
const ALL_FRAGMENTS = [...COMMON_FRAGMENTS, ...HETEROCYCLIC_FRAGMENTS];

/**
 * Predict 13C shifts with atom index tracking.
 * Returns atom indices for error attribution.
 */
export async function predictShiftsWithAtomIndices(smiles) {
  const hoseCodes = smilesToHoseCodes(smiles, { nucleus: '13C' });
  await preloadChunks(hoseCodes.map(e => e.hose));

  const results = [];
  for (const entry of hoseCodes) {
    let hit = null;
    let hoseToUse = entry.hose;

    hit = await queryHose(hoseToUse);

    if (!hit) {
      let truncated = hoseToUse;
      for (let attempt = 0; attempt < 8 && !hit; attempt++) {
        const lastDelimIdx = Math.max(
          truncated.lastIndexOf('/'),
          truncated.lastIndexOf(')'),
          truncated.lastIndexOf('('),
        );
        if (lastDelimIdx <= 0) break;
        const beforeDelim = truncated.substring(0, lastDelimIdx);
        const delim = truncated[lastDelimIdx];
        truncated = beforeDelim + delim;
        hit = await queryHose(truncated);
        if (hit) { hoseToUse = truncated; break; }
        truncated = beforeDelim;
        hit = await queryHose(truncated);
        if (hit) { hoseToUse = truncated; break; }
      }
    }

    if (!hit && hoseToUse.match(/^H+/)) {
      const withoutH = hoseToUse.replace(/^H+/, '');
      hit = await queryHose(withoutH);
      if (hit) hoseToUse = withoutH;
    }

    if (hit) {
      results.push({
        shift: hit.avgShift,
        atomIndex: entry.index,
        hose: hoseToUse,
      });
    }
  }
  return results;
}

/**
 * Compute loss between predicted and target 13C NMR shifts.
 * Uses greedy nearest-neighbor assignment.
 */
export function computeLoss(predictedShifts, targetShifts, options = {}) {
  const { unmatchedPenalty = 50 } = options;

  if (predictedShifts.length === 0 && targetShifts.length === 0) {
    return { loss: 0, assignments: [], unmatchedTarget: [], unmatchedPred: [], predAtomError: new Map() };
  }

  if (predictedShifts.length === 0) {
    return {
      loss: targetShifts.length * unmatchedPenalty,
      assignments: [],
      unmatchedTarget: targetShifts.map((_, i) => i),
      unmatchedPred: [],
      predAtomError: new Map(),
    };
  }

  const predSorted = predictedShifts.map((s, i) => ({ val: s, origIdx: i }))
    .sort((a, b) => a.val - b.val);
  const targetSorted = targetShifts.map((s, i) => ({ val: s, origIdx: i }))
    .sort((a, b) => a.val - b.val);

  const usedPred = new Set();
  const assignments = [];

  for (const target of targetSorted) {
    let bestDist = Infinity;
    let bestPred = null;

    for (const pred of predSorted) {
      if (usedPred.has(pred.origIdx)) continue;
      const dist = Math.abs(pred.val - target.val);
      if (dist < bestDist) {
        bestDist = dist;
        bestPred = pred;
      }
    }

    if (bestPred !== null) {
      usedPred.add(bestPred.origIdx);
      assignments.push({
        targetIdx: target.origIdx,
        predIdx: bestPred.origIdx,
        error: bestDist,
      });
    }
  }

  const matchedTarget = new Set(assignments.map(a => a.targetIdx));
  const unmatchedTarget = targetSorted
    .filter(t => !matchedTarget.has(t.origIdx))
    .map(t => t.origIdx);
  const unmatchedPred = predSorted
    .filter(p => !usedPred.has(p.origIdx))
    .map(p => p.origIdx);

  const matchedError = assignments.reduce((sum, a) => sum + a.error, 0);
  const penalty = (unmatchedTarget.length + unmatchedPred.length) * unmatchedPenalty;
  const loss = matchedError + penalty;

  const predAtomError = new Map();
  for (const a of assignments) {
    predAtomError.set(a.predIdx, a.error);
  }
  for (const idx of unmatchedPred) {
    predAtomError.set(idx, unmatchedPenalty);
  }

  return { loss, assignments, unmatchedTarget, unmatchedPred, predAtomError };
}

/**
 * Identify the top-K predicted atoms with highest error contribution.
 */
export function identifyWorstAtoms(lossResult, prediction, topK = 5) {
  const atomErrors = [];

  for (const [predIdx, error] of lossResult.predAtomError) {
    if (predIdx < prediction.length) {
      atomErrors.push({ atomIndex: prediction[predIdx].atomIndex, error });
    }
  }

  atomErrors.sort((a, b) => b.error - a.error);
  return atomErrors.slice(0, topK).map(e => e.atomIndex);
}

/**
 * Identify a molecule whose predicted 13C NMR spectrum matches the target shifts.
 * Uses steepest descent with HOSE-guided fragment mutations.
 */
export async function identifyMolecule(targetShifts, options = {}) {
  const {
    startSmiles = 'C',
    maxSteps = 50,
    topK = 5,
    unmatchedPenalty = 50,
    convergenceThreshold = 0.01,
    onStep = null,
    maxBacktracks = 5,
    timeoutMs = 0,
    targetMW = 0,
    mwTolerance = 0.3,
    allowedAtoms = null,
  } = options;

  const allowedAtomSet = allowedAtoms ? new Set(allowedAtoms) : null;

  const lossOpts = { unmatchedPenalty };
  const startTime = timeoutMs > 0 ? Date.now() : 0;

  let currentSmiles = startSmiles;
  let currentPrediction = await globalCachedPredict(currentSmiles);
  let currentShifts = currentPrediction.map(p => p.shift);
  let currentLoss = computeLoss(currentShifts, targetShifts, lossOpts);

  const trajectory = [{
    step: 0,
    smiles: currentSmiles,
    loss: currentLoss.loss,
    mutation: 'initial',
    predictedShifts: [...currentShifts],
  }];

  // Early return if already solved
  if (currentLoss.loss <= convergenceThreshold) {
    return {
      smiles: currentSmiles,
      predictedShifts: currentShifts,
      loss: currentLoss.loss,
      steps: 0,
      trajectory,
    };
  }

  let plateauCount = 0;
  let backtracksUsed = 0;
  const visitedSmiles = new Set([currentSmiles]);

  for (let step = 1; step <= maxSteps; step++) {
    // Time budget check
    if (timeoutMs > 0 && (Date.now() - startTime) >= timeoutMs) break;

    // Use static fragment list (avoids expensive reverse HOSE lookups)
    const guidedFragments = ALL_FRAGMENTS;

    // Staged enumeration: additive mutations first (fast path for structural recovery)
    const filterCandidates = (raw) => {
      let cands = raw.filter(c => !visitedSmiles.has(c.smiles));
      if (allowedAtomSet) cands = cands.filter(c => containsOnlyAtoms(c.smiles, allowedAtomSet));
      if (targetMW > 0) {
        const maxMW = targetMW * (1 + mwTolerance);
        cands = cands.filter(c => { const mw = getMolecularWeight(c.smiles); return mw > 0 && mw <= maxMW; });
      }
      return cands;
    };

    // Enumerate mutations in two stages: additive first (fast path), then structural tweaks
    let candidates = filterCandidates([
      ...enumerateAdditions(currentSmiles),
      ...enumerateInlineInsertions(currentSmiles),
      ...enumerateFragmentAttachments(currentSmiles, guidedFragments),
    ]);
    // Flag: if stage 1 finds loss=0, skip stage 2 enumeration entirely
    let needStage2 = true;

    if (candidates.length === 0) {
      // No candidates after tabu + MW + atom filter — jump to earliest unexhausted point
      if (backtracksUsed < maxBacktracks && trajectory.length >= 2) {
        // When candidates are exhausted, go directly to the start for maximum diversity
        const rewindPoint = trajectory[0];
        currentSmiles = rewindPoint.smiles;
        currentPrediction = await globalCachedPredict(currentSmiles);
        currentShifts = currentPrediction.map(p => p.shift);
        currentLoss = computeLoss(currentShifts, targetShifts, lossOpts);
        plateauCount = 0;
        backtracksUsed++;
        trajectory.push({
          step, smiles: currentSmiles, loss: currentLoss.loss,
          mutation: `[backtrack to step ${rewindPoint.step}]`,
          predictedShifts: [...currentShifts],
        });
        if (onStep) onStep(trajectory[trajectory.length - 1]);
        continue;
      }
      break;
    }

    // Sort by carbon count closeness, pre-filter to top 128, then normalize-dedup
    const targetCarbons = targetShifts.length;
    const currentCarbons = countCarbons(currentSmiles);
    candidates.sort((a, b) => Math.abs(countCarbons(a.smiles) - targetCarbons) - Math.abs(countCarbons(b.smiles) - targetCarbons));
    if (candidates.length > 128) candidates.length = 128;
    // Normalize-dedup: eliminates non-canonical duplicates from fastValidate
    const seenNorm = new Set();
    candidates = candidates.filter(c => {
      try {
        const norm = normalize(c.smiles);
        if (seenNorm.has(norm)) return false;
        seenNorm.add(norm);
        c.smiles = norm;
        return true;
      } catch { return false; }
    });
    if (candidates.length > 64) candidates.length = 64;

    // Evaluate candidates in parallel batches, keeping top-K for lookahead
    let bestCandidate = null;
    let bestLoss = currentLoss.loss;
    let leastWorseCand = null;
    let leastWorseLoss = Infinity;
    const topK_evaluated = []; // best candidate per mutation type for lookahead

    const BATCH_SIZE = 64;
    let earlyExit = false;

    // Helper to process evaluation results
    const processResults = (results) => {
      for (const result of results) {
        if (result.status !== 'fulfilled' || !result.value) continue;
        const { candidate, pred, shifts, loss } = result.value;
        const evaluated = {
          smiles: candidate.smiles,
          description: candidate.description,
          prediction: pred,
          shifts,
          lossResult: loss,
        };
        if (loss.loss < bestLoss) {
          bestLoss = loss.loss;
          bestCandidate = evaluated;
        }
        if (loss.loss <= convergenceThreshold) { earlyExit = true; }
        if (loss.loss < leastWorseLoss && loss.loss >= currentLoss.loss) {
          leastWorseLoss = loss.loss;
          leastWorseCand = evaluated;
        }
        // Track best candidate per mutation type for diverse lookahead
        const descType = candidate.description.split(' ')[0];
        const existing = topK_evaluated.find(e => e.description.split(' ')[0] === descType);
        if (!existing || loss.loss < existing.lossResult.loss) {
          if (existing) topK_evaluated.splice(topK_evaluated.indexOf(existing), 1);
          topK_evaluated.push(evaluated);
        }
      }
    };

    // Evaluate all candidates (capped at 100)
    for (let bStart = 0; bStart < candidates.length; bStart += BATCH_SIZE) {
      if (earlyExit) break;
      if (timeoutMs > 0 && (Date.now() - startTime) >= timeoutMs) break;
      const batch = candidates.slice(bStart, bStart + BATCH_SIZE);
      const results = await Promise.allSettled(
        batch.map(async (candidate) => {
          const pred = await globalCachedPredict(candidate.smiles);
          if (!pred) return null;
          const shifts = pred.map(p => p.shift);
          const loss = computeLoss(shifts, targetShifts, lossOpts);
          return { candidate, pred, shifts, loss };
        })
      );
      processResults(results);
    }

    // Stage 2: structural tweaks (only when carbon count already correct AND stage 1 didn't find loss=0)
    if (!earlyExit && needStage2 && currentCarbons === targetCarbons) {
      let stage2 = filterCandidates([
        ...enumerateSubstitutions(currentSmiles),
        ...enumerateBondChanges(currentSmiles),
        ...enumerateRemovals(currentSmiles),
        ...enumerateRingClosures(currentSmiles),
      ]);
      stage2.sort((a, b) => Math.abs(countCarbons(a.smiles) - targetCarbons) - Math.abs(countCarbons(b.smiles) - targetCarbons));
      if (stage2.length > 64) stage2.length = 64;
      for (let bStart = 0; bStart < stage2.length; bStart += BATCH_SIZE) {
        if (earlyExit) break;
        if (timeoutMs > 0 && (Date.now() - startTime) >= timeoutMs) break;
        const batch = stage2.slice(bStart, bStart + BATCH_SIZE);
        const results = await Promise.allSettled(
          batch.map(async (candidate) => {
            const pred = await globalCachedPredict(candidate.smiles);
            if (!pred) return null;
            const shifts = pred.map(p => p.shift);
            const loss = computeLoss(shifts, targetShifts, lossOpts);
            return { candidate, pred, shifts, loss };
          })
        );
        processResults(results);
      }
    }

    // Sort by loss for lookahead priority
    topK_evaluated.sort((a, b) => a.lossResult.loss - b.lossResult.loss);

    // 2-step lookahead: check if any top candidate leads to a better 2-step path
    let lookaheadOverride = false;
    const elapsed = timeoutMs > 0 ? Date.now() - startTime : 0;
    if (bestLoss > convergenceThreshold && topK_evaluated.length > 1 &&
        !(timeoutMs > 0 && elapsed >= timeoutMs * 0.9)) {
      let bestStep2Loss = bestLoss;
      let bestIntermediate = null;
      const lookaheadBudget = timeoutMs > 0 ? (timeoutMs - elapsed) * 0.4 : Infinity;
      const lookaheadEnd = Date.now() + lookaheadBudget;

      // Check non-best candidates in REVERSE loss order (most different from best = most potential)
      // Use only additive mutations for step 2 (fast enumeration)
      const lookaheadCands = topK_evaluated.slice(1).reverse();
      for (const cand of lookaheadCands) {
        if (timeoutMs > 0 && Date.now() >= lookaheadEnd) break;
        // Fast: only additive mutations for step 2
        let step2Cands = filterCandidates([
          ...enumerateAdditions(cand.smiles),
          ...enumerateInlineInsertions(cand.smiles),
          ...enumerateFragmentAttachments(cand.smiles, guidedFragments),
        ]);
        step2Cands = step2Cands.filter(c => c.smiles !== cand.smiles);
        // Only evaluate exact carbon count matches
        step2Cands = step2Cands.filter(c => countCarbons(c.smiles) === targetCarbons);
        const subset = step2Cands.slice(0, 32);
        if (subset.length === 0) continue;
        const step2Results = await Promise.allSettled(
          subset.map(async (c2) => {
            const pred = await globalCachedPredict(c2.smiles);
            if (!pred) return null;
            const shifts = pred.map(p => p.shift);
            const loss = computeLoss(shifts, targetShifts, lossOpts);
            return { loss: loss.loss };
          })
        );
        for (const r of step2Results) {
          if (r.status !== 'fulfilled' || !r.value) continue;
          if (r.value.loss < bestStep2Loss) {
            bestStep2Loss = r.value.loss;
            bestIntermediate = cand;
          }
        }
      }

      // Only accept lookahead if 2-step path is significantly better (>30% lower loss)
      if (bestIntermediate && bestStep2Loss < bestLoss * 0.7) {
        bestCandidate = bestIntermediate;
        bestLoss = bestIntermediate.lossResult.loss;
        bestCandidate.description = `[lookahead] ${bestCandidate.description}`;
        lookaheadOverride = true;
        plateauCount = 0;
      }
    }

    // Check convergence (skip if lookahead override forced an intermediate step)
    if (!lookaheadOverride && (bestCandidate === null || (currentLoss.loss - bestLoss) < convergenceThreshold)) {
      plateauCount++;
      const tolerance = [1.2, 1.5, 2.0][Math.min(plateauCount - 1, 2)];
      if (plateauCount <= 3 && leastWorseCand && leastWorseLoss <= currentLoss.loss * tolerance) {
        bestCandidate = leastWorseCand;
        bestLoss = leastWorseLoss;
        bestCandidate.description = `[perturb] ${bestCandidate.description}`;
      } else if (backtracksUsed < maxBacktracks && trajectory.length >= 3) {
        const rewindFractions = [0.6, 0.3, 0.1, 0.0];
        const fraction = rewindFractions[Math.min(backtracksUsed, rewindFractions.length - 1)];
        const rewindIdx = Math.max(0, Math.floor(trajectory.length * fraction));
        const rewindPoint = trajectory[rewindIdx];
        currentSmiles = rewindPoint.smiles;
        currentPrediction = await globalCachedPredict(currentSmiles);
        currentShifts = currentPrediction.map(p => p.shift);
        currentLoss = computeLoss(currentShifts, targetShifts, lossOpts);
        plateauCount = 0;
        backtracksUsed++;
        trajectory.push({
          step, smiles: currentSmiles, loss: currentLoss.loss,
          mutation: `[backtrack to step ${rewindPoint.step}]`,
          predictedShifts: [...currentShifts],
        });
        if (onStep) onStep(trajectory[trajectory.length - 1]);
        continue;
      } else {
        break;
      }
    } else {
      plateauCount = 0;
    }

    // Accept best mutation
    currentSmiles = bestCandidate.smiles;
    visitedSmiles.add(currentSmiles);
    currentPrediction = bestCandidate.prediction;
    currentShifts = bestCandidate.shifts;
    currentLoss = bestCandidate.lossResult;

    const stepInfo = {
      step,
      smiles: currentSmiles,
      loss: currentLoss.loss,
      mutation: bestCandidate.description,
      predictedShifts: [...currentShifts],
    };
    trajectory.push(stepInfo);

    if (onStep) onStep(stepInfo);

    if (currentLoss.loss < convergenceThreshold) break;
  }

  // Return best point across trajectory
  let bestPoint = trajectory[0];
  for (const point of trajectory) {
    if (point.loss < bestPoint.loss) bestPoint = point;
  }

  let bestPrediction = currentPrediction;
  if (bestPoint.smiles !== currentSmiles) {
    bestPrediction = await globalCachedPredict(bestPoint.smiles);
  }

  return {
    smiles: bestPoint.smiles,
    loss: bestPoint.loss,
    steps: trajectory.length - 1,
    trajectory,
    predictedShifts: bestPrediction,
  };
}
