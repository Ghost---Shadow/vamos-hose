// Acetaminophen and Related Compounds
// Non-NSAID analgesics and antipyretics

import {
  Fragment, Ring, Linear, Molecule,
} from '../src/index.js';

// === FRAGMENTS ===

// Acetamido: CH3-C(=O)-NH- -> CC(=O)N
export const acetamido = Fragment('CC(=O)N');

// Benzene ring
export const benzene = Ring({ atoms: 'c', size: 6 });

// Hydroxyl group
export const hydroxyl = Linear(['O']);

// Para-hydroxybenzene (phenol at position 4)
export const paraHydroxybenzene = benzene.attach(4, hydroxyl);

// === ACETAMINOPHEN ===
// SMILES: CC(=O)NC1=CC=C(O)C=C1

export const acetaminophen = Fragment('CC(=O)Nc1ccc(O)cc1');
export const acetaminophenBuilt = Molecule([acetamido, paraHydroxybenzene]);
export const paracetamol = Fragment('CC(=O)Nc1ccc(O)cc1');

// === PHENACETIN ===
// SMILES: CC(=O)NC1=CC=C(OCC)C=C1

// Ethoxy group: -O-CH2-CH3
export const ethoxy = Linear(['O', 'C', 'C']);

// Para-ethoxybenzene (ethoxy at position 4)
export const paraEthoxybenzene = Ring({ atoms: 'c', size: 6 }).attach(4, ethoxy);

export const phenacetin = Fragment('CC(=O)Nc1ccc(OCC)cc1');
export const phenacetinBuilt = Molecule([acetamido, paraEthoxybenzene]);
