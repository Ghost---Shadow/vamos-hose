import { smilesToHoseCodes } from './smiles-to-hose.js';
import { queryHose, preloadChunks } from './database.js';
import { enumerateAllMutations, getAtomCount } from './mutate.js';
import { estimateFromSpectra, resolveHoseSmiles } from './estimate.js';

const COMMON_FRAGMENTS = ['C', 'CC', 'O', 'N', 'C(=O)O', 'NC', 'S'];

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

  if (peaksToLookup.length === 0) return [...COMMON_FRAGMENTS];

  // Dedupe peaks (round to nearest integer for PPM bucket lookup)
  const uniquePeaks = [...new Set(peaksToLookup.map(p => Math.round(p)))];

  try {
    const hoseResults = await estimateFromSpectra({
      peaks: uniquePeaks,
      tolerance: 5.0,
      maxResults: 3,
    });

    const enriched = await resolveHoseSmiles(hoseResults);

    const fragments = new Set(COMMON_FRAGMENTS);
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
 * Uses steepest descent with HOSE-guided fragment mutations via smiles-js.
 */
export async function identifyMolecule(targetShifts, options = {}) {
  const {
    startSmiles = 'C',
    maxSteps = 50,
    topK = 5,
    unmatchedPenalty = 50,
    convergenceThreshold = 0.01,
    onStep = null,
    maxBacktracks = 3,
    timeoutMs = 0,
  } = options;

  const lossOpts = { unmatchedPenalty };
  const startTime = timeoutMs > 0 ? Date.now() : 0;

  let currentSmiles = startSmiles;
  let currentPrediction = await predictShiftsWithAtomIndices(currentSmiles);
  let currentShifts = currentPrediction.map(p => p.shift);
  let currentLoss = computeLoss(currentShifts, targetShifts, lossOpts);

  const trajectory = [{
    step: 0,
    smiles: currentSmiles,
    loss: currentLoss.loss,
    mutation: 'initial',
    predictedShifts: [...currentShifts],
  }];

  let plateauCount = 0;
  let backtracksUsed = 0;
  const visitedSmiles = new Set([currentSmiles]);

  for (let step = 1; step <= maxSteps; step++) {
    // Time budget check
    if (timeoutMs > 0 && (Date.now() - startTime) >= timeoutMs) break;

    // Get HOSE-guided fragments based on poorly-matched peaks
    const guidedFragments = await getGuidedFragments(targetShifts, currentLoss, currentPrediction);

    // Enumerate mutations with HOSE-guided fragments
    let candidates = enumerateAllMutations(currentSmiles, guidedFragments);

    // Filter tabu list
    candidates = candidates.filter(c => !visitedSmiles.has(c.smiles));

    if (candidates.length === 0) break;

    // Evaluate candidates in parallel batches
    let bestCandidate = null;
    let bestLoss = currentLoss.loss;
    let leastWorseCand = null;
    let leastWorseLoss = Infinity;

    const BATCH_SIZE = 32;
    for (let bStart = 0; bStart < candidates.length; bStart += BATCH_SIZE) {
      if (timeoutMs > 0 && (Date.now() - startTime) >= timeoutMs) break;
      const batch = candidates.slice(bStart, bStart + BATCH_SIZE);
      const results = await Promise.allSettled(
        batch.map(async (candidate) => {
          const pred = await predictShiftsWithAtomIndices(candidate.smiles);
          const shifts = pred.map(p => p.shift);
          const loss = computeLoss(shifts, targetShifts, lossOpts);
          return { candidate, pred, shifts, loss };
        })
      );

      for (const result of results) {
        if (result.status !== 'fulfilled') continue;
        const { candidate, pred, shifts, loss } = result.value;
        if (loss.loss < bestLoss) {
          bestLoss = loss.loss;
          bestCandidate = {
            smiles: candidate.smiles,
            description: candidate.description,
            prediction: pred,
            shifts,
            lossResult: loss,
          };
        }
        if (loss.loss < leastWorseLoss && loss.loss >= currentLoss.loss) {
          leastWorseLoss = loss.loss;
          leastWorseCand = { smiles: candidate.smiles, description: candidate.description, prediction: pred, shifts, lossResult: loss };
        }
      }
    }

    // Check convergence
    if (bestCandidate === null || (currentLoss.loss - bestLoss) < convergenceThreshold) {
      plateauCount++;
      const tolerance = [1.2, 1.5, 2.0][Math.min(plateauCount - 1, 2)];
      if (plateauCount <= 3 && leastWorseCand && leastWorseLoss <= currentLoss.loss * tolerance) {
        bestCandidate = leastWorseCand;
        bestLoss = leastWorseLoss;
        bestCandidate.description = `[perturb] ${bestCandidate.description}`;
      } else if (backtracksUsed < maxBacktracks && trajectory.length >= 3) {
        const rewindFractions = [0.6, 0.3, 0.1];
        const fraction = rewindFractions[Math.min(backtracksUsed, rewindFractions.length - 1)];
        const rewindIdx = Math.max(1, Math.floor(trajectory.length * fraction));
        const rewindPoint = trajectory[rewindIdx];
        currentSmiles = rewindPoint.smiles;
        currentPrediction = await predictShiftsWithAtomIndices(currentSmiles);
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
    bestPrediction = await predictShiftsWithAtomIndices(bestPoint.smiles);
  }

  return {
    smiles: bestPoint.smiles,
    loss: bestPoint.loss,
    steps: trajectory.length - 1,
    trajectory,
    predictedShifts: bestPrediction,
  };
}
