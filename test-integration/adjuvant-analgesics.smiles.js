// Adjuvant Analgesics and Neuropathic Pain Medications
// Drugs not designed as painkillers but help with certain pain types

import { Fragment } from 'smiles-js';

// === GABAPENTIN ===
// Structural analog of GABA (gamma-aminobutyric acid)
// SMILES: C1CCC(CC1)(CC(=O)O)CN
// Used for neuropathic pain and seizures

export const gabapentin = Fragment('C1CCC(CC1)(CC(=O)O)CN');

// === PREGABALIN ===
// (S)-3-(aminomethyl)-5-methylhexanoic acid
// SMILES: CC(C)CC(CC(=O)O)CN
// More potent than gabapentin, binds alpha-2-delta subunit

export const pregabalin = Fragment('CC(C)CC(CC(=O)O)CN');

// === AMITRIPTYLINE ===
// Tricyclic antidepressant used for neuropathic pain
// SMILES: CN(C)CCC=C1C2=CC=CC=C2CCC3=CC=CC=C31

export const amitriptyline = Fragment('CN(C)CCC=C1C2=CC=CC=C2CCC3=CC=CC=C31');

// === DULOXETINE ===
// SNRI (serotonin-norepinephrine reuptake inhibitor)
// SMILES: CNCCC(C1=CC=CS1)OC2=CC=CC3=CC=CC=C32
// Used for diabetic neuropathy, fibromyalgia

export const duloxetine = Fragment('CNCCC(C1=CC=CS1)OC2=CC=CC3=CC=CC=C32');

// === CARBAMAZEPINE ===
// Anticonvulsant used for trigeminal neuralgia
// SMILES: C1=CC=C2C(=C1)C=CC3=CC=CC=C3N2C(=O)N

export const carbamazepine = Fragment('C1=CC=C2C(=C1)C=CC3=CC=CC=C3N2C(=O)N');

// === VALPROIC ACID ===
// 2-propylpentanoic acid (branched fatty acid)
// SMILES: CCCC(CCC)C(=O)O
// Anticonvulsant, used for migraine prophylaxis

export const valproicAcid = Fragment('CCCC(CCC)C(=O)O');
