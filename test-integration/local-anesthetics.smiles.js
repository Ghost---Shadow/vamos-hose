// Local Anesthetics and Topical Analgesics
// Pain relief through nerve block or topical application

import { Fragment } from 'smiles-js';

// === LIDOCAINE ===
// Amide-type local anesthetic (first of this class)
// SMILES: CCN(CC)CC(=O)NC1=C(C)C=CC=C1C

export const lidocaine = Fragment('CCN(CC)CC(=O)NC1=C(C)C=CC=C1C');

// === BUPIVACAINE ===
// Long-acting amide anesthetic
// SMILES: CCCCN1CCCCC1C(=O)NC2=C(C)C=CC=C2C

export const bupivacaine = Fragment('CCCCN1CCCCC1C(=O)NC2=C(C)C=CC=C2C');

// === ROPIVACAINE ===
// Similar to bupivacaine but with n-propyl instead of n-butyl
// SMILES: CCCN1CCCCC1C(=O)NC2=C(C)C=CC=C2C
// Less cardiotoxic than bupivacaine

export const ropivacaine = Fragment('CCCN1CCCCC1C(=O)NC2=C(C)C=CC=C2C');

// === MEPIVACAINE ===
// Intermediate-duration amide anesthetic
// SMILES: CN1CCCCC1C(=O)NC2=C(C)C=CC=C2C

export const mepivacaine = Fragment('CN1CCCCC1C(=O)NC2=C(C)C=CC=C2C');

// === PRILOCAINE ===
// Amide anesthetic with lower toxicity
// SMILES: CCCNC(C)C(=O)NC1=CC=CC=C1C

export const prilocaine = Fragment('CCCNC(C)C(=O)NC1=CC=CC=C1C');

// === BENZOCAINE ===
// Ester-type topical anesthetic (no systemic use due to poor water solubility)
// SMILES: CCOC(=O)C1=CC=C(C=C1)N

export const benzocaine = Fragment('CCOC(=O)C1=CC=C(C=C1)N');

// === TETRACAINE ===
// Long-acting ester anesthetic
// SMILES: CCCCNC1=CC=C(C=C1)C(=O)OCCN(C)C

export const tetracaine = Fragment('CCCCNC1=CC=C(C=C1)C(=O)OCCN(C)C');

// === PROCAINE ===
// Classic ester anesthetic (Novocain)
// SMILES: CCN(CC)CCOC(=O)C1=CC=C(C=C1)N

export const procaine = Fragment('CCN(CC)CCOC(=O)C1=CC=C(C=C1)N');
