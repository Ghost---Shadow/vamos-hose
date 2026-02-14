// Endocannabinoids and Cannabinoid System Analgesics
// Endogenous and exogenous cannabinoid receptor agonists for pain relief

import { Fragment } from 'smiles-js';

// === ANANDAMIDE ===
// N-arachidonoylethanolamine (AEA) - endogenous CB1/CB2 agonist
// SMILES: CCCCCC=CCC=CCC=CCC=CCCCC(=O)NCCO
// "Bliss molecule" - partial agonist at CB1 receptors

export const anandamide = Fragment('CCCCCC=CCC=CCC=CCC=CCCCC(=O)NCCO');

// === 2-ARACHIDONOYLGLYCEROL ===
// 2-AG - most abundant endocannabinoid
// SMILES: CCCCCC=CCC=CCC=CCC=CCCCC(=O)OC(CO)CO
// Full agonist at both CB1 and CB2 receptors

export const arachidonoylglycerol2 = Fragment('CCCCCC=CCC=CCC=CCC=CCCCC(=O)OC(CO)CO');

// === THC (Δ9-TETRAHYDROCANNABINOL) ===
// Primary psychoactive component of cannabis
// SMILES: CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O
// Partial agonist at CB1/CB2 receptors

export const thc = Fragment('CCCCCC1=CC(=C2C3C=C(CCC3C(OC2=C1)(C)C)C)O');

// === CBD (CANNABIDIOL) ===
// Non-psychoactive cannabinoid with analgesic properties
// SMILES: CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O
// Negative allosteric modulator of CB1, indirect effects on endocannabinoid system

export const cbd = Fragment('CCCCCC1=CC(=C(C(=C1)O)C2C=C(CCC2C(=C)C)C)O');

// === NABILONE ===
// Synthetic cannabinoid analog of THC
// SMILES: CCCCCCC(C)(C)C1=CC(=C2C3CC(=O)CCC3C(OC2=C1)(C)C)O
// Used for chemotherapy-induced nausea and neuropathic pain

export const nabilone = Fragment('CCCCCCC(C)(C)C1=CC(=C2C3CC(=O)CCC3C(OC2=C1)(C)C)O');

// === PALMITOYLETHANOLAMIDE ===
// PEA - endogenous fatty acid amide with anti-inflammatory properties
// SMILES: CCCCCCCCCCCCCCCC(=O)NCCO
// Indirect cannabinoid activity via entourage effect

export const palmitoylethanolamide = Fragment('CCCCCCCCCCCCCCCC(=O)NCCO');
