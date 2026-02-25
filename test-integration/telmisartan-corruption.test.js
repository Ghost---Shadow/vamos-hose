import { describe, test, expect, beforeAll } from 'bun:test';
import { identifyMolecule, predictShiftsWithAtomIndices, computeLoss } from '../src/identify.js';
import { getMolecularWeight } from '../src/mutate.js';
import { normalize } from 'smiles-js';
// DB loads lazily on demand — no preload needed

const telmisartan = 'CCCC1=NC2=C(C=C(C=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C(=O)O)C5=NC6=CC=CC=C6N5C)C';

let targetShifts;
let targetMagnitude;
let fivePctThreshold;
let targetMW;

beforeAll(async () => {
  const pred = await predictShiftsWithAtomIndices(telmisartan);
  targetShifts = pred.map(p => p.shift);
  targetMagnitude = targetShifts.reduce((s, v) => s + Math.abs(v), 0);
  fivePctThreshold = targetMagnitude * 0.05;
  targetMW = getMolecularWeight(telmisartan);
  console.log(`Target: ${targetShifts.length} carbons, MW=${targetMW.toFixed(0)}, magnitude=${targetMagnitude.toFixed(0)}, 5%=${fivePctThreshold.toFixed(1)}`);
  console.log(`Target shifts: [${targetShifts.sort((a, b) => a - b).map(s => s.toFixed(1)).join(', ')}]`);
});

// Telmisartan contains only C, N, O (plus H implicit)
const TELMISARTAN_ATOMS = ['C', 'c', 'N', 'n', 'O', 'o'];

async function runCorruption(name, startSmiles, maxSteps = 25, timeoutMs = 30000) {
  let validStart;
  try {
    validStart = normalize(startSmiles);
  } catch {
    console.log(`  [${name}] INVALID SMILES: ${startSmiles}`);
    return null;
  }

  const startPred = await predictShiftsWithAtomIndices(validStart);
  const startShifts = startPred.map(p => p.shift);
  const startLoss = computeLoss(startShifts, targetShifts, { unmatchedPenalty: 50 });
  const startMW = getMolecularWeight(validStart);

  console.log(`\n=== ${name} ===`);
  console.log(`  Start SMILES: ${validStart}`);
  console.log(`  Start MW: ${startMW.toFixed(0)} (target: ${targetMW.toFixed(0)})`);
  console.log(`  Start carbons: ${startShifts.length} (target: ${targetShifts.length})`);
  console.log(`  Start loss: ${startLoss.loss.toFixed(1)}`);
  console.log(`  Start shifts: [${startShifts.sort((a, b) => a - b).map(s => s.toFixed(1)).join(', ')}]`);

  const result = await identifyMolecule(targetShifts, {
    startSmiles: validStart,
    maxSteps,
    topK: 5,
    unmatchedPenalty: 50,
    timeoutMs,
    targetMW,
    allowedAtoms: TELMISARTAN_ATOMS,
  });

  const finalMW = getMolecularWeight(result.smiles);
  console.log(`  Final SMILES: ${result.smiles}`);
  console.log(`  Final MW: ${finalMW.toFixed(0)}, carbons: ${result.predictedShifts.length}`);
  console.log(`  Final loss: ${result.loss.toFixed(1)} (threshold: ${fivePctThreshold.toFixed(1)})`);
  console.log(`  Steps: ${result.steps}`);
  console.log(`  ${result.loss <= fivePctThreshold ? 'PASS' : 'FAIL'}`);

  // Show trajectory
  for (const step of result.trajectory) {
    const mw = getMolecularWeight(step.smiles);
    console.log(`    Step ${String(step.step).padStart(2)}: loss=${step.loss.toFixed(1).padStart(7)} C=${String((step.predictedShifts || []).length).padStart(2)} MW=${mw.toFixed(0).padStart(4)} "${step.mutation}"`);
  }

  // Show worst matches if failed
  if (result.loss > fivePctThreshold) {
    const finalPred = await predictShiftsWithAtomIndices(result.smiles);
    const fs = finalPred.map(p => p.shift);
    const fl = computeLoss(fs, targetShifts, { unmatchedPenalty: 50 });
    console.log(`  Unmatched target: ${fl.unmatchedTarget.length}, unmatched pred: ${fl.unmatchedPred.length}`);
    const sorted = [...fl.assignments].sort((a, b) => b.error - a.error);
    for (const a of sorted.slice(0, 5)) {
      console.log(`    pred=${fs[a.predIdx]?.toFixed(1)} ↔ target=${targetShifts[a.targetIdx]?.toFixed(1)} err=${a.error.toFixed(1)}`);
    }
  }

  return result;
}

describe('Telmisartan Corruption Tests', () => {
  test('0: Exact telmisartan (sanity check - should be ~0 loss immediately)', async () => {
    const result = await runCorruption('Exact telmisartan', telmisartan, 5, 5000);
    expect(result).not.toBeNull();
    expect(result.loss).toBeLessThan(5);
    expect(result.steps).toBeLessThanOrEqual(1);
  }, 10000);

  test('1: Missing propyl chain (CCC→C, -2 carbons)', async () => {
    const corrupted = telmisartan.replace('CCCC1', 'CC1');
    const result = await runCorruption('Missing propyl (CCC→C)', corrupted, 10, 10000);
    expect(result).not.toBeNull();
    expect(result.loss).toBeLessThan(fivePctThreshold);
  }, 15000);

  test('2: Missing carboxylic acid (C(=O)O removed)', async () => {
    const corrupted = telmisartan.replace('C(=O)O', '');
    const result = await runCorruption('Missing COOH', corrupted, 10, 10000);
    expect(result).not.toBeNull();
    expect(result.loss).toBeLessThan(fivePctThreshold);
  }, 15000);

  test('3: Missing aryl methyl (terminal CH3)', async () => {
    const corrupted = telmisartan.replace(/\)C$/, ')');
    const result = await runCorruption('Missing aryl-CH3', corrupted, 10, 10000);
    expect(result).not.toBeNull();
    expect(result.loss).toBeLessThan(fivePctThreshold);
  }, 15000);

  test('4: Missing N-methyl on benzimidazole', async () => {
    const corrupted = telmisartan.replace('N5C)', 'N5)');
    const result = await runCorruption('Missing N-CH3', corrupted, 10, 10000);
    expect(result).not.toBeNull();
    expect(result.loss).toBeLessThan(fivePctThreshold);
  }, 15000);

  test('5: Both methyls removed (-2 carbons)', async () => {
    let corrupted = telmisartan.replace('N5C)', 'N5)');
    corrupted = corrupted.replace(/\)C$/, ')');
    const result = await runCorruption('Both methyls removed', corrupted, 10, 15000);
    expect(result).not.toBeNull();
    expect(result.loss).toBeLessThan(fivePctThreshold);
  }, 20000);

  test('6: Missing second benzimidazole (replaced with H)', async () => {
    const corrupted = telmisartan.replace('C5=NC6=CC=CC=C6N5C', 'C');
    const result = await runCorruption('Missing 2nd benzimidazole', corrupted, 10, 15000);
    expect(result).not.toBeNull();
    console.log(`  Loss ratio: ${(result.loss / fivePctThreshold * 100).toFixed(0)}% of threshold`);
  }, 20000);

  test('7: N→C swap in first benzimidazole (wrong heterocycle)', async () => {
    let corrupted = telmisartan.replace('=NC2=', '=CC2=');
    corrupted = corrupted.replace('N1CC3', 'C1CC3');
    const result = await runCorruption('N→C in 1st benzimidazole', corrupted, 10, 15000);
    expect(result).not.toBeNull();
    console.log(`  Loss ratio: ${(result.loss / fivePctThreshold * 100).toFixed(0)}% of threshold`);
  }, 20000);

  test('8: Just the biphenyl-COOH fragment', async () => {
    const corrupted = 'c1ccc(cc1)c2ccccc2C(=O)O';
    const result = await runCorruption('Biphenyl-COOH only', corrupted, 10, 15000);
    expect(result).not.toBeNull();
    console.log(`  Loss ratio: ${(result.loss / fivePctThreshold * 100).toFixed(0)}% of threshold`);
  }, 20000);
});
