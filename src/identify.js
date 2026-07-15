import { smilesToHoseCodes } from './smiles-to-hose.js';
import { queryHose, preloadChunks } from './database.js';
import {
  enumerateAllMutations, enumerateAdditions, enumerateInlineInsertions,
  enumerateFragmentAttachments, enumerateSubstitutions, enumerateBondChanges,
  enumerateRemovals, enumerateRingClosures, enumerateAtomSwaps,
  getAtomCount, getMolecularWeight, containsOnlyAtoms, countCarbons,
  getMolecularFormula, formulaBudget, fragmentFitsBudget, renumberRingDigits, atomToElement,
} from './mutate.js';
import { estimateFromSpectra, resolveHoseSmiles } from './estimate.js';

import { buildSMILES, Ring, FusedRing, normalize } from 'smiles-js';
import {
  methyl, ethyl, propyl, hydroxyl, carboxyl, amino, cyano,
  benzene, imidazole, pyridine, pyrrole, phenyl,
} from 'smiles-js/common';

// Build fragment SMILES strings from smiles-js constructors (guaranteed valid)
const COMMON_FRAGMENTS = [
  'C', 'CC', 'CCC', 'O', 'N', 'C(=O)O', 'NC',
];

// Heterocyclic ring fragments built from smiles-js Ring constructors
const HETEROCYCLIC_FRAGMENTS = [
  buildSMILES(imidazole),    // c1nc[nH]c1
  buildSMILES(pyridine),     // c1ccncc1
  buildSMILES(pyrrole),      // c1cc[nH]c1
  buildSMILES(phenyl),       // c1ccccc1
  buildSMILES(carboxyl),     // C(=O)O
  buildSMILES(amino),        // N
  buildSMILES(cyano),        // C#N
  'c1nc2ccccc2[nH]1',        // benzimidazole aromatic (7C, 2N)
  'c1nc2ccccc2n1C',           // N-methylbenzimidazole aromatic (8C, 2N)
  'c1ccc2[nH]cnc2c1',        // benzimidazole (alt notation)
  'c1ccc2nc[nH]c2c1',        // benzimidazole (alt notation 2)
  'C1=NC2=CC=CC=C2N1',       // benzimidazole Kekulé (7C, 2N)
  'C1=NC2=CC=CC=C2N1C',      // N-methylbenzimidazole Kekulé (8C, 2N)
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
 * Batch-attach: collect all budget-fitting fragments, try each at every atom position,
 * evaluate OVERALL loss (not per-peak). This handles large fragments like benzimidazole
 * that cover multiple missing peaks at once.
 *
 * @param {boolean} includeShifted - if true, also include poorly-matched peaks (for round 2+)
 * Returns array of candidate SMILES with descriptions.
 */
async function batchAttachForMissingPeaks(currentSmiles, currentLoss, targetShifts, targetFormula, allowedAtomSet, includeShifted = false) {
  const currentFormula = getMolecularFormula(currentSmiles);
  const budget = formulaBudget(targetFormula, currentFormula);

  const hasPositiveBudget = Object.values(budget).some(v => v > 0);

  // Strategy 0: Pairwise atom swaps — when budget has both negative and positive entries
  // e.g., {C:-2, N:+2} means swap 2 C's to N's
  const swapCandidates = [];
  const negElements = Object.entries(budget).filter(([, v]) => v < 0);
  const posElements = Object.entries(budget).filter(([, v]) => v > 0);
  if (negElements.length > 0 && posElements.length > 0) {
    const positions = findAtomPositionsForAttach(currentSmiles);
    // Handle simple case: swap element A→B at K positions
    for (const [fromEl, deficit] of negElements) {
      for (const [toEl, surplus] of posElements) {
        if (Math.abs(deficit) !== surplus) continue; // only exact swaps
        const swapCount = surplus;
        const fromLower = fromEl.toLowerCase();
        const fromPositions = positions.filter(p => atomToElement(p.atom) === fromEl);
        if (fromPositions.length === 0 || swapCount > fromPositions.length) continue;

        if (swapCount === 1) {
          // Single swap: try each position
          for (const pos of fromPositions) {
            const toAtom = pos.atom === pos.atom.toLowerCase() ? toEl.toLowerCase() : toEl;
            const s = currentSmiles.slice(0, pos.start) + toAtom + currentSmiles.slice(pos.end);
            const validated = fastValidate(s);
            if (!validated) continue;
            try {
              const norm = normalize(validated);
              swapCandidates.push({ smiles: norm, description: `swap ${fromEl}→${toEl} @${positions.indexOf(pos)}` });
            } catch { /* invalid */ }
          }
        } else if (swapCount === 2) {
          // Pairwise: try all pairs
          for (let i = 0; i < fromPositions.length; i++) {
            for (let j = i + 1; j < fromPositions.length; j++) {
              const pi = fromPositions[i], pj = fromPositions[j];
              const toI = pi.atom === pi.atom.toLowerCase() ? toEl.toLowerCase() : toEl;
              const toJ = pj.atom === pj.atom.toLowerCase() ? toEl.toLowerCase() : toEl;
              let s;
              if (pi.start < pj.start) {
                s = currentSmiles.slice(0, pi.start) + toI + currentSmiles.slice(pi.end, pj.start) + toJ + currentSmiles.slice(pj.end);
              } else {
                s = currentSmiles.slice(0, pj.start) + toJ + currentSmiles.slice(pj.end, pi.start) + toI + currentSmiles.slice(pi.end);
              }
              const validated = fastValidate(s);
              if (!validated) continue;
              try {
                const norm = normalize(validated);
                swapCandidates.push({ smiles: norm, description: `swap ${fromEl}→${toEl} x2` });
              } catch { /* invalid */ }
            }
          }
        } else if (swapCount <= 4) {
          // K-swap: enumerate K-tuples (cap at reasonable count)
          const indices = fromPositions.map((_, idx) => idx);
          const combos = [];
          const combo = (start, chosen) => {
            if (chosen.length === swapCount) { combos.push([...chosen]); return; }
            if (combos.length >= 500) return;
            for (let k = start; k < indices.length; k++) combo(k + 1, [...chosen, k]);
          };
          combo(0, []);
          for (const c of combos) {
            const posArr = c.map(k => fromPositions[k]).sort((a, b) => b.start - a.start);
            let s = currentSmiles;
            for (const pos of posArr) {
              const toAtom = pos.atom === pos.atom.toLowerCase() ? toEl.toLowerCase() : toEl;
              s = s.slice(0, pos.start) + toAtom + s.slice(pos.end);
            }
            const validated = fastValidate(s);
            if (!validated) continue;
            try {
              const norm = normalize(validated);
              swapCandidates.push({ smiles: norm, description: `swap ${fromEl}→${toEl} x${swapCount}` });
            } catch { /* invalid */ }
          }
        }
      }
    }
  }

  if (!hasPositiveBudget && !includeShifted && swapCandidates.length === 0) return [];

  // Collect budget-fitting fragments for branch attachment
  const allFrags = new Set();
  // Separate set for inline-replacement fragments (can exceed branch budget by 1 atom)
  const replaceFrags = new Set();

  // Relaxed budget: for inline replacement, replacing 1 atom gives back 1 unit
  const maxReplaceBudget = {};
  for (const el of Object.keys(budget)) maxReplaceBudget[el] = (budget[el] || 0) + 1;

  function addFrag(frag) {
    if (allowedAtomSet && !containsOnlyAtoms(frag, allowedAtomSet)) return;
    const ff = getMolecularFormula(frag);
    if (fragmentFitsBudget(ff, budget)) {
      allFrags.add(frag);
    } else if (fragmentFitsBudget(ff, maxReplaceBudget)) {
      replaceFrags.add(frag);
    }
  }

  // Simple atoms and chains
  for (const atom of ['C', 'N', 'O']) {
    if (!allowedAtomSet || allowedAtomSet.has(atom)) addFrag(atom);
  }
  if ((budget['C'] || 0) >= 2) allFrags.add('CC');
  if ((budget['C'] || 0) >= 3) allFrags.add('CCC');

  // ALL_FRAGMENTS (heterocyclics, common groups)
  for (const frag of ALL_FRAGMENTS) addFrag(frag);

  // Reverse HOSE lookup for unmatched peaks
  const missingShifts = currentLoss.unmatchedTarget.map(idx => Math.round(targetShifts[idx]));
  if (includeShifted) {
    const sorted = [...currentLoss.assignments].sort((a, b) => b.error - a.error);
    for (const a of sorted.slice(0, 3)) {
      if (a.error > 5.0) missingShifts.push(Math.round(targetShifts[a.targetIdx]));
    }
  }
  const uniqueShifts = [...new Set(missingShifts)];

  for (const shift of uniqueShifts.slice(0, 5)) {
    try {
      const hoseResults = await estimateFromSpectra({ peaks: [shift], tolerance: 5.0, maxResults: 5 });
      const enriched = await resolveHoseSmiles(hoseResults);
      for (const result of Object.values(enriched)) {
        if (!result || !result.smiles) continue;
        for (const frag of extractFragmentsFromSmiles(result.smiles)) {
          addFrag(frag);
        }
      }
    } catch { /* skip */ }
  }

  if (allFrags.size === 0 && replaceFrags.size === 0) return [];

  // Sort fragments: largest first (fragments covering more of the budget are more valuable)
  const sortBySize = (a, b) => {
    const fa = Object.values(getMolecularFormula(a)).reduce((s, v) => s + v, 0);
    const fb = Object.values(getMolecularFormula(b)).reduce((s, v) => s + v, 0);
    return fb - fa;
  };
  const fragList = [...allFrags].sort(sortBySize);
  const replaceFragList = [...allFrags, ...replaceFrags].sort(sortBySize);

  // Try each fragment at every atom position, generate candidates
  const candidates = [];
  const positions = findAtomPositionsForAttach(currentSmiles);

  // Strategy 1: Branch attachment — (fragment) after atom
  for (const frag of fragList) {
    const rn = renumberRingDigits(frag, currentSmiles);
    if (!rn) continue;
    for (let posIdx = 0; posIdx < positions.length; posIdx++) {
      const pos = positions[posIdx];
      const testSmiles = currentSmiles.slice(0, pos.end) + `(${rn})` + currentSmiles.slice(pos.end);
      const validated = fastValidate(testSmiles);
      if (!validated) continue;
      try {
        const norm = normalize(validated);
        candidates.push({ smiles: norm, description: `batch-attach (${frag}) @${posIdx}` });
      } catch { /* invalid */ }
    }
  }

  // Strategy 2: Inline replacement — replace atom with fragment
  // Net formula change = fragFormula - {replacedAtomElement: 1}
  // Uses replaceFragList which includes fragments 1 atom larger than branch budget
  for (const frag of replaceFragList) {
    const rn = renumberRingDigits(frag, currentSmiles);
    if (!rn) continue;
    const fragFormula = getMolecularFormula(frag);
    for (let posIdx = 0; posIdx < positions.length; posIdx++) {
      const pos = positions[posIdx];
      const atomEl = atomToElement(pos.atom);
      // Check net formula change fits budget (directional: can't over-add or wrong-direction)
      let netFits = true;
      const checkedEls = new Set();
      for (const [atom, count] of Object.entries(fragFormula)) {
        checkedEls.add(atom);
        const net = atom === atomEl ? count - 1 : count;
        const b = budget[atom] || 0;
        if (net < Math.min(0, b) || net > Math.max(0, b)) { netFits = false; break; }
      }
      // Check replaced element if not in fragment (net = -1 for that element)
      if (netFits && !checkedEls.has(atomEl)) {
        const b = budget[atomEl] || 0;
        if (-1 < Math.min(0, b) || -1 > Math.max(0, b)) netFits = false;
      }
      if (!netFits) continue;
      const testSmiles = currentSmiles.slice(0, pos.start) + rn + currentSmiles.slice(pos.end);
      const validated = fastValidate(testSmiles);
      if (!validated) continue;
      try {
        const norm = normalize(validated);
        candidates.push({ smiles: norm, description: `batch-replace ${pos.atom}->(${frag}) @${posIdx}` });
      } catch { /* invalid */ }
    }
  }

  // Merge swap candidates with attachment candidates
  const allCandidates = [...swapCandidates, ...candidates];

  // Dedupe
  const seen = new Set();
  return allCandidates.filter(c => {
    if (seen.has(c.smiles)) return false;
    seen.add(c.smiles);
    return true;
  });
}

/** Attach a fragment at the end of a SMILES string, return fastValidate'd result or null. Handles ring digit renumbering. */
function fastValidateAttach(smiles, frag) {
  const rn = renumberRingDigits(frag, smiles) || frag;
  const atomRe = new RegExp(/\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnosp]/.source, 'g');
  let lastMatch = null;
  let m;
  while ((m = atomRe.exec(smiles)) !== null) lastMatch = m;
  if (!lastMatch) return null;
  const insertPos = lastMatch.index + lastMatch[0].length;
  const result = smiles.slice(0, insertPos) + `(${rn})` + smiles.slice(insertPos);
  return fastValidate(result);
}

/** Find atom positions for branch attachment. */
function findAtomPositionsForAttach(smiles) {
  const re = new RegExp(/\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnosp]/.source, 'g');
  const positions = [];
  let m;
  while ((m = re.exec(smiles)) !== null) {
    positions.push({ start: m.index, end: m.index + m[0].length, atom: m[0] });
  }
  return positions;
}

/** Fast validation: balanced parens + non-empty. */
function fastValidate(smiles) {
  if (!smiles || smiles.length === 0) return null;
  let depth = 0;
  for (let i = 0; i < smiles.length; i++) {
    if (smiles[i] === '(') depth++;
    else if (smiles[i] === ')') { depth--; if (depth < 0) return null; }
  }
  return depth === 0 ? smiles : null;
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
    maxBacktracks = 5,
    timeoutMs = 0,
    targetMW = 0,
    mwTolerance = 0.3,
    allowedAtoms = null,
    targetFormula = null,
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

  // ── Phase 1: Beam-search batch-attach (when targetFormula is provided) ───
  if (targetFormula && currentLoss.loss > convergenceThreshold) {
    const BEAM_K = 3;
    const maxBatchRounds = Math.min(maxSteps, 5);

    let beams = [{
      smiles: currentSmiles,
      prediction: currentPrediction,
      shifts: currentShifts,
      lossResult: currentLoss,
      description: 'initial',
    }];

    for (let batchRound = 0; batchRound < maxBatchRounds; batchRound++) {
      if (timeoutMs > 0 && (Date.now() - startTime) >= timeoutMs) break;
      if (beams.some(b => b.lossResult.loss <= convergenceThreshold)) break;

      const includeShifted = batchRound > 0;
      const nextBeams = [];
      const seenSmiles = new Set(beams.map(b => b.smiles));

      for (const beam of beams) {
        if (timeoutMs > 0 && (Date.now() - startTime) >= timeoutMs) break;

        const batchCandidates = await batchAttachForMissingPeaks(
          beam.smiles, beam.lossResult, targetShifts, targetFormula, allowedAtomSet, includeShifted
        );
        if (batchCandidates.length === 0) continue;

        const evaluated = [];
        for (const cand of batchCandidates) {
          const pred = await globalCachedPredict(cand.smiles);
          if (!pred) continue;
          const shifts = pred.map(p => p.shift);
          const loss = computeLoss(shifts, targetShifts, lossOpts);
          if (loss.loss < beam.lossResult.loss) {
            evaluated.push({ smiles: cand.smiles, description: cand.description, prediction: pred, shifts, lossResult: loss });
          }
        }

        evaluated.sort((a, b) => a.lossResult.loss - b.lossResult.loss);
        let taken = 0;
        for (const e of evaluated) {
          if (taken >= BEAM_K) break;
          if (seenSmiles.has(e.smiles)) continue;
          seenSmiles.add(e.smiles);
          nextBeams.push(e);
          taken++;
        }
      }

      if (nextBeams.length === 0) break;

      nextBeams.sort((a, b) => a.lossResult.loss - b.lossResult.loss);
      beams = nextBeams.slice(0, BEAM_K);

      // Log best beam's step to trajectory
      const bestNow = beams[0];
      const stepInfo = {
        step: trajectory.length,
        smiles: bestNow.smiles,
        loss: bestNow.lossResult.loss,
        mutation: bestNow.description,
        predictedShifts: [...bestNow.shifts],
      };
      trajectory.push(stepInfo);
      if (onStep) onStep(stepInfo);
    }

    // Pick globally best beam
    const bestBeam = beams.reduce((best, b) => b.lossResult.loss < best.lossResult.loss ? b : best, beams[0]);
    currentSmiles = bestBeam.smiles;
    currentPrediction = bestBeam.prediction;
    currentShifts = bestBeam.shifts;
    currentLoss = bestBeam.lossResult;

    if (currentLoss.loss <= convergenceThreshold) {
      return {
        smiles: currentSmiles,
        predictedShifts: currentPrediction,
        loss: currentLoss.loss,
        steps: trajectory.length - 1,
        trajectory,
      };
    }
  }

  // ── Phase 2: Greedy gradient descent (fallback / refinement) ──────────────
  let plateauCount = 0;
  let backtracksUsed = 0;
  const visitedSmiles = new Set([currentSmiles]);

  const remainingSteps = maxSteps - (trajectory.length - 1);
  for (let step = 1; step <= remainingSteps; step++) {
    // Time budget check
    if (timeoutMs > 0 && (Date.now() - startTime) >= timeoutMs) break;

    // Use static fragment list (avoids expensive reverse HOSE lookups)
    const guidedFragments = ALL_FRAGMENTS;

    // Staged enumeration: additive mutations first (fast path for structural recovery)
    const filterCandidates = (raw) => {
      let cands = raw.filter(c => !visitedSmiles.has(c.smiles));
      if (allowedAtomSet) cands = cands.filter(c => containsOnlyAtoms(c.smiles, allowedAtomSet));
      if (targetFormula) {
        cands = cands.filter(c => {
          const cf = getMolecularFormula(c.smiles);
          for (const [atom, count] of Object.entries(cf)) {
            if (count > (targetFormula[atom] || 0)) return false;
          }
          return true;
        });
      } else if (targetMW > 0) {
        const maxMW = targetMW * (1 + mwTolerance);
        cands = cands.filter(c => { const mw = getMolecularWeight(c.smiles); return mw > 0 && mw <= maxMW; });
      }
      return cands;
    };

    // Enumerate mutations in two stages: additive first (fast path), then structural tweaks
    const targetCarbons = targetShifts.length;
    const currentCarbons = countCarbons(currentSmiles);

    let candidates = filterCandidates([
      ...enumerateAdditions(currentSmiles),
      ...enumerateInlineInsertions(currentSmiles),
      ...enumerateFragmentAttachments(currentSmiles, guidedFragments),
    ]);
    // Flag: if stage 1 finds loss=0, skip stage 2 enumeration entirely
    let needStage2 = true;

    // If stage 1 is empty and carbon count matches, try stage 2 immediately (swaps, subs)
    if (candidates.length === 0 && currentCarbons === targetCarbons) {
      candidates = filterCandidates([
        ...enumerateSubstitutions(currentSmiles),
        ...enumerateBondChanges(currentSmiles),
        ...enumerateRemovals(currentSmiles),
        ...enumerateRingClosures(currentSmiles),
        ...enumerateAtomSwaps(currentSmiles),
      ]);
      needStage2 = false; // already ran stage 2
    }

    if (candidates.length === 0) {
      // No candidates after tabu + MW + atom filter — jump to earliest unexhausted point
      if (backtracksUsed < maxBacktracks && trajectory.length >= 2) {
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
        ...enumerateAtomSwaps(currentSmiles),
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
