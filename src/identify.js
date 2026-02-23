import { smilesToHoseCodes } from './smiles-to-hose.js';
import { queryHose, preloadChunks } from './database.js';
import { enumerateAllMutations, getAllHeavyAtomIndices, getAllCarbonIndices } from './mutate.js';

/**
 * Predict 13C shifts with atom index tracking.
 * Unlike lookupNmrShifts, this returns atom indices for error attribution.
 *
 * @param {string} smiles
 * @returns {Promise<Array<{ shift: number, atomIndex: number, hose: string }>>}
 */
export async function predictShiftsWithAtomIndices(smiles) {
  const hoseCodes = smilesToHoseCodes(smiles, { nucleus: '13C' });
  await preloadChunks(hoseCodes.map(e => e.hose));

  const results = [];
  for (const entry of hoseCodes) {
    let hit = null;
    let hoseToUse = entry.hose;

    // Exact match
    hit = await queryHose(hoseToUse);

    // Truncation fallback (same logic as lookup.js)
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

    // Try without leading H's
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
 *
 * @param {number[]} predictedShifts - predicted shift values
 * @param {number[]} targetShifts - target shift values
 * @param {object} [options]
 * @param {number} [options.unmatchedPenalty=50] - penalty per unmatched peak (ppm)
 * @returns {{ loss: number, assignments: Array, unmatchedTarget: number[], unmatchedPred: number[], predAtomError: Map<number, number> }}
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

  // For each target peak, greedily assign nearest unassigned predicted peak
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

  // Unmatched indices
  const matchedTarget = new Set(assignments.map(a => a.targetIdx));
  const unmatchedTarget = targetSorted
    .filter(t => !matchedTarget.has(t.origIdx))
    .map(t => t.origIdx);
  const unmatchedPred = predSorted
    .filter(p => !usedPred.has(p.origIdx))
    .map(p => p.origIdx);

  // Total loss
  const matchedError = assignments.reduce((sum, a) => sum + a.error, 0);
  const penalty = (unmatchedTarget.length + unmatchedPred.length) * unmatchedPenalty;
  const loss = matchedError + penalty;

  // Per-predicted-atom error map
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
 * Returns their molecule atom indices for focused mutation.
 *
 * @param {{ predAtomError: Map<number, number> }} lossResult
 * @param {Array<{ shift: number, atomIndex: number }>} prediction
 * @param {number} [topK=5]
 * @returns {number[]} molecule atom indices with highest error
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
 * Uses steepest descent with discrete molecular mutations.
 *
 * @param {number[]} targetShifts - target 13C NMR chemical shifts (ppm)
 * @param {object} [options]
 * @param {string} [options.startSmiles='C'] - starting molecule SMILES
 * @param {number} [options.maxSteps=50] - maximum optimization steps
 * @param {number} [options.topK=5] - number of worst atoms to focus mutations on
 * @param {number} [options.unmatchedPenalty=50] - penalty per unmatched peak
 * @param {number} [options.convergenceThreshold=0.01] - stop if loss improvement < this
 * @param {boolean} [options.widenOnPlateau=true] - try all-atom mutations if focused fails
 * @param {function} [options.onStep] - callback(stepInfo) for progress tracking
 * @returns {Promise<{ smiles: string, loss: number, steps: number, trajectory: Array, predictedShifts: Array }>}
 */
export async function identifyMolecule(targetShifts, options = {}) {
  const {
    startSmiles = 'C',
    maxSteps = 50,
    topK = 5,
    unmatchedPenalty = 50,
    convergenceThreshold = 0.01,
    widenOnPlateau = true,
    onStep = null,
  } = options;

  const lossOpts = { unmatchedPenalty };

  // Initial prediction
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

  for (let step = 1; step <= maxSteps; step++) {
    // Error attribution: identify worst atoms
    const worstAtomIndices = identifyWorstAtoms(currentLoss, currentPrediction, topK);

    let atomIndicesToMutate = worstAtomIndices.length > 0
      ? worstAtomIndices
      : getAllCarbonIndices(currentSmiles);

    // Enumerate focused mutations
    let candidates = enumerateAllMutations(currentSmiles, atomIndicesToMutate);

    // Widen if no focused candidates found
    if (candidates.length === 0 && widenOnPlateau) {
      atomIndicesToMutate = getAllHeavyAtomIndices(currentSmiles);
      candidates = enumerateAllMutations(currentSmiles, atomIndicesToMutate);
    }

    if (candidates.length === 0) break;

    // Evaluate all candidates
    let bestCandidate = null;
    let bestLoss = currentLoss.loss;

    for (const candidate of candidates) {
      try {
        const pred = await predictShiftsWithAtomIndices(candidate.smiles);
        const shifts = pred.map(p => p.shift);
        const loss = computeLoss(shifts, targetShifts, lossOpts);

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
      } catch {
        // Skip candidates that fail prediction
      }
    }

    // Check convergence
    if (bestCandidate === null || (currentLoss.loss - bestLoss) < convergenceThreshold) {
      break;
    }

    // Accept best mutation
    currentSmiles = bestCandidate.smiles;
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

    // Early termination: perfect match
    if (currentLoss.loss < convergenceThreshold) break;
  }

  return {
    smiles: currentSmiles,
    loss: currentLoss.loss,
    steps: trajectory.length - 1,
    trajectory,
    predictedShifts: currentPrediction,
  };
}
