// Over-the-Counter NSAIDs
// Common painkillers available without prescription

import { Fragment } from 'smiles-js';

// === ASPIRIN ===
// Acetylsalicylic acid
// SMILES: CC(=O)Oc1ccccc1C(=O)O

export const aspirin = Fragment('CC(=O)Oc1ccccc1C(=O)O');

// === IBUPROFEN ===
// 2-(4-isobutylphenyl)propionic acid
// SMILES: CC(C)Cc1ccc(cc1)C(C)C(=O)O

export const ibuprofen = Fragment('CC(C)Cc1ccc(cc1)C(C)C(=O)O');

// === NAPROXEN ===
// (S)-2-(6-methoxynaphthalen-2-yl)propionic acid
// SMILES: COc1ccc2cc(ccc2c1)C(C)C(=O)O

export const naproxen = Fragment('COc1ccc2cc(ccc2c1)C(C)C(=O)O');

// === KETOPROFEN ===
// 2-(3-benzoylphenyl)propionic acid
// SMILES: CC(c1cccc(c1)C(=O)c2ccccc2)C(=O)O

export const ketoprofen = Fragment('CC(c1cccc(c1)C(=O)c2ccccc2)C(=O)O');

// === DICLOFENAC ===
// 2-[2-(2,6-dichloroanilino)phenyl]acetic acid
// SMILES: OC(=O)Cc1ccccc1Nc2c(Cl)cccc2Cl

export const diclofenac = Fragment('OC(=O)Cc1ccccc1Nc2c(Cl)cccc2Cl');
