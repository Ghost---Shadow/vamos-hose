// Endogenous Opioid Peptides
// Natural peptides produced by the body that bind opioid receptors

import { Fragment } from 'smiles-js';

// === MET-ENKEPHALIN ===
// Tyr-Gly-Gly-Phe-Met - first endogenous opioid discovered (1975)
// SMILES: CSCCC(C(=O)O)NC(=O)C(CC1=CC=CC=C1)NC(=O)CNC(=O)CNC(=O)C(CC2=CC=C(C=C2)O)N
// Binds delta opioid receptors, involved in pain modulation

export const metEnkephalin = Fragment('CSCCC(C(=O)O)NC(=O)C(CC1=CC=CC=C1)NC(=O)CNC(=O)CNC(=O)C(CC2=CC=C(C=C2)O)N');

// === LEU-ENKEPHALIN ===
// Tyr-Gly-Gly-Phe-Leu - differs from met-enkephalin only at C-terminal residue
// SMILES: CC(C)CC(C(=O)O)NC(=O)C(CC1=CC=CC=C1)NC(=O)CNC(=O)CNC(=O)C(CC2=CC=C(C=C2)O)N
// Also binds delta opioid receptors

export const leuEnkephalin = Fragment('CC(C)CC(C(=O)O)NC(=O)C(CC1=CC=CC=C1)NC(=O)CNC(=O)CNC(=O)C(CC2=CC=C(C=C2)O)N');

// === ENDOMORPHIN-1 ===
// Tyr-Pro-Trp-Phe-NH2 - highest affinity endogenous mu-opioid receptor ligand
// SMILES: C1CC(N(C1)C(=O)C(CC2=CC=C(C=C2)O)N)C(=O)NC(CC3=CNC4=CC=CC=C43)C(=O)NC(CC5=CC=CC=C5)C(=O)N
// Selective for mu receptors (same target as morphine)

export const endomorphin1 = Fragment('C1CC(N(C1)C(=O)C(CC2=CC=C(C=C2)O)N)C(=O)NC(CC3=CNC4=CC=CC=C43)C(=O)NC(CC5=CC=CC=C5)C(=O)N');

// === ENDOMORPHIN-2 ===
// Tyr-Pro-Phe-Phe-NH2 - differs from endomorphin-1 at position 3 (Phe instead of Trp)
// SMILES: C1CC(N(C1)C(=O)C(CC2=CC=C(C=C2)O)N)C(=O)NC(CC3=CC=CC=C3)C(=O)NC(CC4=CC=CC=C4)C(=O)N
// Also selective for mu receptors

export const endomorphin2 = Fragment('C1CC(N(C1)C(=O)C(CC2=CC=C(C=C2)O)N)C(=O)NC(CC3=CC=CC=C3)C(=O)NC(CC4=CC=CC=C4)C(=O)N');
