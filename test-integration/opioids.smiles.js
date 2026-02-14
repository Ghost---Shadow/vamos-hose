// Opioid Analgesics
// Strong prescription painkillers acting on opioid receptors

import { Fragment } from 'smiles-js';

// === MORPHINE ===
// Natural opioid alkaloid - base structure for many opioids
// SMILES: CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O
// Has phenolic OH at C3 and allylic OH at C6

export const morphine = Fragment('CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O');

// === CODEINE ===
// Morphine with C3 phenolic OH methylated to methoxy
// SMILES: CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O
// Less potent than morphine, used for cough suppression

export const codeine = Fragment('CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O');

// === OXYCODONE ===
// Semi-synthetic: C6 ketone, C14 hydroxyl, C3 methoxy
// SMILES: CN1CCC23C4C(=O)CCC2(C1CC5=C3C(=C(C=C5)OC)O4)O
// More potent than morphine

export const oxycodone = Fragment('CN1CCC23C4C(=O)CCC2(C1CC5=C3C(=C(C=C5)OC)O4)O');

// === HYDROCODONE ===
// Semi-synthetic: C6 ketone, C3 methoxy (differs from oxycodone by lacking C14 hydroxyl)
// SMILES: CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4
// CHEMISTRY FIX: Correct structure differs from oxycodone

export const hydrocodone = Fragment('CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4');

// === FENTANYL ===
// Synthetic phenylpiperidine (NOT morphinan structure)
// SMILES: CCC(=O)N(C1CCN(CC1)CCC2=CC=CC=C2)C3=CC=CC=C3
// 50-100x more potent than morphine

export const fentanyl = Fragment('CCC(=O)N(C1CCN(CC1)CCC2=CC=CC=C2)C3=CC=CC=C3');

// === TRAMADOL ===
// Synthetic centrally acting analgesic (weak opioid + SNRI activity)
// SMILES: CN(C)CC1CCCCC1(C2=CC(=CC=C2)OC)O
// Cyclohexanol with dimethylaminomethyl and m-methoxyphenyl substituents

export const tramadol = Fragment('CN(C)CC1CCCCC1(C2=CC(=CC=C2)OC)O');

// === METHADONE ===
// Synthetic acyclic opioid (lacks morphinan structure)
// SMILES: CCC(=O)C(CC(C)N(C)C)(C1=CC=CC=C1)C2=CC=CC=C2
// Used for opioid addiction treatment

export const methadone = Fragment('CCC(=O)C(CC(C)N(C)C)(C1=CC=CC=C1)C2=CC=CC=C2');
