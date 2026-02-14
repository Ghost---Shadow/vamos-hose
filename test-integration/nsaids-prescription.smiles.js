// Prescription NSAIDs and COX-2 Inhibitors
// Stronger anti-inflammatory drugs requiring prescription

import { Fragment } from 'smiles-js';

// === CELECOXIB ===
// COX-2 selective inhibitor
// SMILES: CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F

export const celecoxib = Fragment('CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F');

// === MELOXICAM ===
// Oxicam class NSAID with thiazole
// SMILES: CC1=C(N=C(S1)NC(=O)C2=C(C3=CC=CC=C3S(=O)(=O)N2C)O)C

export const meloxicam = Fragment('CC1=C(N=C(S1)NC(=O)C2=C(C3=CC=CC=C3S(=O)(=O)N2C)O)C');

// === PIROXICAM ===
// Oxicam class NSAID
// SMILES: CN1C(=C(C2=CC=CC=C2S1(=O)=O)O)C(=O)NC3=CC=CC=N3

export const piroxicam = Fragment('CN1C(=C(C2=CC=CC=C2S1(=O)=O)O)C(=O)NC3=CC=CC=N3');

// === ETODOLAC ===
// Indole acetic acid derivative
// SMILES: CCC1=CC2=C(C=C1CC(=O)O)NC3=C2CCOC3(CC)CC

export const etodolac = Fragment('CCC1=CC2=C(C=C1CC(=O)O)NC3=C2CCOC3(CC)CC');

// === KETOROLAC ===
// Pyrrolizine carboxylic acid with benzoyl
// SMILES: OC(=O)C1CCN2C1=CC=C2C(=O)C3=CC=CC=C3

export const ketorolac = Fragment('OC(=O)C1CCN2C1=CC=C2C(=O)C3=CC=CC=C3');

// === ROFECOXIB (Withdrawn 2004) ===
// COX-2 selective (withdrawn due to cardiovascular risk)
// SMILES: CS(=O)(=O)C1=CC=C(C=C1)C2=C(C(=O)OC2)C3=CC=CC=C3

export const rofecoxib = Fragment('CS(=O)(=O)C1=CC=C(C=C1)C2=C(C(=O)OC2)C3=CC=CC=C3');

// === ETORICOXIB ===
// COX-2 selective inhibitor
// SMILES: CC1=NC=C(C=C1)C2=CC=C(C=C2)S(=O)(=O)C3=CC=CC=C3

export const etoricoxib = Fragment('CC1=NC=C(C=C1)C2=CC=C(C=C2)S(=O)(=O)C3=CC=CC=C3');

// === NABUMETONE ===
// Prodrug (converted to 6-methoxy-2-naphthylacetic acid)
// SMILES: COC1=CC2=CC(=CC=C2C=C1)CCC(=O)C

export const nabumetone = Fragment('COC1=CC2=CC(=CC=C2C=C1)CCC(=O)C');

// === OXAPROZIN ===
// Long half-life NSAID with oxazole ring
// SMILES: OC(=O)CCC1=NC(=C(O1)C2=CC=CC=C2)C3=CC=CC=C3

export const oxaprozin = Fragment('OC(=O)CCC1=NC(=C(O1)C2=CC=CC=C2)C3=CC=CC=C3');
