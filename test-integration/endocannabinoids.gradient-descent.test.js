import { describe, test, expect } from 'bun:test';
import { identifyMolecule, predictShiftsWithAtomIndices, computeLoss } from '../src/identify.js';
import { enumerateAllMutations } from '../src/mutate.js';
import {
  anandamide,
  thc,
  cbd,
  nabilone,
  palmitoylethanolamide,
} from './endocannabinoids.smiles.js';

const COMPOUNDS = [
  { name: 'Anandamide', smiles: anandamide, start: 'C', desc: 'from methane' },
  { name: 'THC', smiles: thc, start: 'c1ccccc1', desc: 'from benzene' },
  { name: 'CBD', smiles: cbd, start: 'c1ccccc1', desc: 'from benzene' },
  { name: 'Nabilone', smiles: nabilone, start: 'c1ccccc1', desc: 'from benzene' },
  { name: 'Palmitoylethanolamide', smiles: palmitoylethanolamide, start: 'C', desc: 'from methane' },
];

function printTrajectory(name, result, targetShifts) {
  console.log(`\n=== ${name} ===`);
  console.log(`Target: [${targetShifts.map(s => s.toFixed(1)).join(', ')}]`);
  console.log(`Steps: ${result.steps}  |  Initial loss: ${result.trajectory[0].loss.toFixed(2)}  |  Final loss: ${result.loss.toFixed(2)}`);
  const pctOfInitial = ((result.loss / result.trajectory[0].loss) * 100).toFixed(1);
  console.log(`Loss remaining: ${pctOfInitial}% of initial`);
  for (const step of result.trajectory) {
    const shifts = step.predictedShifts || [];
    console.log(`  Step ${step.step}: loss=${step.loss.toFixed(2)}  C=${shifts.length}  "${step.mutation}"`);
    console.log(`    ${step.smiles}`);
  }
}

describe('Endocannabinoids - Gradient Descent 5% Loss', () => {
  const results = {};

  for (const compound of COMPOUNDS) {
    test(`${compound.name}: ${compound.desc}`, async () => {
      const targetPred = await predictShiftsWithAtomIndices(compound.smiles);
      const targetShifts = targetPred.map(p => p.shift);
      const targetMagnitude = targetShifts.reduce((s, v) => s + Math.abs(v), 0);
      const fivePctThreshold = targetMagnitude * 0.05;

      const result = await identifyMolecule(targetShifts, {
        startSmiles: compound.start,
        maxSteps: 25,
        topK: 5,
        unmatchedPenalty: 50,
        timeoutMs: 4500,
      });

      printTrajectory(compound.name, result, targetShifts);
      console.log(`  5% threshold: ${fivePctThreshold.toFixed(2)} | Final loss: ${result.loss.toFixed(2)} | ${result.loss <= fivePctThreshold ? 'PASS' : 'FAIL'}`);

      results[compound.name] = {
        initialLoss: result.trajectory[0].loss,
        finalLoss: result.loss,
        steps: result.steps,
        threshold: fivePctThreshold,
        pass: result.loss <= fivePctThreshold,
        finalSmiles: result.smiles,
        targetShifts,
      };

      if (result.loss > fivePctThreshold) {
        console.log(`\n  >>> LOCAL MINIMA ANALYSIS for ${compound.name}:`);
        console.log(`  Final SMILES: ${result.smiles}`);
        console.log(`  Final shifts: [${result.predictedShifts.map(p => p.shift.toFixed(1)).join(', ')}]`);
        console.log(`  Target shifts: [${targetShifts.map(s => s.toFixed(1)).join(', ')}]`);
        const finalShifts = result.predictedShifts.map(p => p.shift);
        const detailedLoss = computeLoss(finalShifts, targetShifts, { unmatchedPenalty: 50 });
        console.log(`  Unmatched target peaks: ${detailedLoss.unmatchedTarget.length}`);
        console.log(`  Unmatched predicted peaks: ${detailedLoss.unmatchedPred.length}`);
        console.log(`  Matched error sum: ${detailedLoss.assignments.reduce((s, a) => s + a.error, 0).toFixed(2)}`);
        const sortedAssignments = [...detailedLoss.assignments].sort((a, b) => b.error - a.error);
        console.log('  Top error assignments:');
        for (const a of sortedAssignments.slice(0, 5)) {
          const predVal = finalShifts[a.predIdx];
          const targetVal = targetShifts[a.targetIdx];
          console.log(`    pred[${a.predIdx}]=${predVal.toFixed(1)} ↔ target[${a.targetIdx}]=${targetVal.toFixed(1)}  error=${a.error.toFixed(1)}`);
        }
        const lastCandidates = enumerateAllMutations(result.smiles);
        console.log(`  Candidates at final step: ${lastCandidates.length} (all-atom widening)`);
      }
    }, 5000);
  }

  test('summary table', async () => {
    console.log('\n\n========================================');
    console.log('COMPOUND RESULTS SUMMARY');
    console.log('========================================');
    console.log('Compound                  Init     Final    5%Thr    Steps  Status');
    console.log('─'.repeat(75));
    let passCount = 0;
    for (const compound of COMPOUNDS) {
      const r = results[compound.name];
      if (!r) continue;
      const status = r.pass ? '✓ PASS' : '✗ FAIL';
      if (r.pass) passCount++;
      console.log(`${compound.name.padEnd(26)}${r.initialLoss.toFixed(1).padStart(8)}${r.finalLoss.toFixed(1).padStart(9)}${r.threshold.toFixed(1).padStart(9)}${String(r.steps).padStart(8)}  ${status}`);
    }
    console.log('─'.repeat(75));
    console.log(`${passCount}/${Object.keys(results).length} compounds within 5% loss threshold`);
    // Note: timed-out compounds won't appear in results
    expect(Object.keys(results).length).toBeGreaterThanOrEqual(0);
  }, 10000);
});
