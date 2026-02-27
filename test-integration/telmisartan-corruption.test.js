import { describe, test, expect, beforeAll } from 'bun:test';
import { identifyMolecule, predictShiftsWithAtomIndices, computeLoss } from '../src/identify.js';
import { getMolecularWeight, getMolecularFormula, countCarbons } from '../src/mutate.js';
import { normalize } from 'smiles-js';
// DB loads lazily on demand — no preload needed

const telmisartan = 'CCCC1=NC2=C(C=C(C=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C(=O)O)C5=NC6=CC=CC=C6N5C)C';

let targetShifts;
let targetMW;
let targetFormula;

beforeAll(async () => {
  const pred = await predictShiftsWithAtomIndices(telmisartan);
  targetShifts = pred.map(p => p.shift);
  targetMW = getMolecularWeight(telmisartan);
  targetFormula = getMolecularFormula(telmisartan);
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
    targetFormula,
    allowedAtoms: TELMISARTAN_ATOMS,
  });

  const finalMW = getMolecularWeight(result.smiles);
  console.log(`  Final SMILES: ${result.smiles}`);
  console.log(`  Final MW: ${finalMW.toFixed(0)}, carbons: ${result.predictedShifts.length}`);
  console.log(`  Final loss: ${result.loss.toFixed(1)}`);
  console.log(`  Steps: ${result.steps}`);

  // Show trajectory
  for (const step of result.trajectory) {
    const mw = getMolecularWeight(step.smiles);
    console.log(`    Step ${String(step.step).padStart(2)}: loss=${step.loss.toFixed(1).padStart(7)} C=${String((step.predictedShifts || []).length).padStart(2)} MW=${mw.toFixed(0).padStart(4)} "${step.mutation}"`);
  }

  return result;
}

describe('Telmisartan Corruption Tests', () => {
  test('0: Exact telmisartan (sanity check)', async () => {
    const result = await runCorruption('Exact telmisartan', telmisartan);
    expect(result.smiles).toEqual(telmisartan);
    expect(result.loss).toEqual(0.0);
  });

  test('1: Missing propyl chain (CCC→C, -2 carbons)', async () => {
    const corrupted = telmisartan.replace('CCCC1', 'CC1');
    const result = await runCorruption('Missing propyl (CCC→C)', corrupted);
    expect(result.loss).toEqual(0.0);
    expect(result.predictedShifts.length).toEqual(targetShifts.length);
  }, 60000);

  test('2: Missing carboxylic acid (C(=O)O removed)', async () => {
    const corrupted = telmisartan.replace('C(=O)O', '');
    const result = await runCorruption('Missing COOH', corrupted);
    expect(result.loss).toEqual(0.0);
    expect(result.predictedShifts.length).toEqual(targetShifts.length);
  }, 60000);

  test('3: Missing aryl methyl (terminal CH3)', async () => {
    const corrupted = telmisartan.replace(/\)C$/, ')');
    const result = await runCorruption('Missing aryl-CH3', corrupted);
    expect(result.loss).toEqual(0.0);
    expect(result.predictedShifts.length).toEqual(targetShifts.length);
  }, 60000);

  test('4: Missing N-methyl on benzimidazole', async () => {
    const corrupted = telmisartan.replace('N5C)', 'N5)');
    const result = await runCorruption('Missing N-CH3', corrupted);
    expect(result.loss).toEqual(0.0);
    expect(result.predictedShifts.length).toEqual(targetShifts.length);
  }, 60000);

  test('5: Both methyls removed (-2 carbons)', async () => {
    let corrupted = telmisartan.replace('N5C)', 'N5)');
    corrupted = corrupted.replace(/\)C$/, ')');
    const result = await runCorruption('Both methyls removed', corrupted);
    expect(result.loss).toEqual(0.0);
    expect(result.predictedShifts.length).toEqual(targetShifts.length);
  }, 60000);

  test('6: Missing second benzimidazole (replaced with H)', async () => {
    const corrupted = telmisartan.replace('C5=NC6=CC=CC=C6N5C', 'C');
    const result = await runCorruption('Missing 2nd benzimidazole', corrupted);
    expect(result.loss).toEqual(0.0);
    expect(result.predictedShifts.length).toEqual(targetShifts.length);
  }, 120000);

  test('7: N→C swap in first benzimidazole (wrong heterocycle)', async () => {
    let corrupted = telmisartan.replace('=NC2=', '=CC2=');
    corrupted = corrupted.replace('N1CC3', 'C1CC3');
    const result = await runCorruption('N→C in 1st benzimidazole', corrupted);
    expect(result.loss).toEqual(0.0);
    expect(result.predictedShifts.length).toEqual(targetShifts.length);
  }, 120000);

  test('8: Just the biphenyl-COOH fragment', async () => {
    // Starting from only 40% of the molecule — exact topology reconstruction is
    // beyond greedy mutation. Target: correct atom count + low average error.
    const corrupted = 'c1ccc(cc1)c2ccccc2C(=O)O';
    const result = await runCorruption('Biphenyl-COOH only', corrupted, 50, 60000);
    expect(result.loss).toBeLessThan(100);
    expect(countCarbons(result.smiles)).toEqual(targetShifts.length);
  }, 120000);
});
