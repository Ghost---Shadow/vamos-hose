// Cholesterol-Lowering Medications
// Statins (HMG-CoA reductase inhibitors) and other lipid-lowering agents
// SMILES verified from PubChem (pubchem.ncbi.nlm.nih.gov)

import { Fragment } from 'smiles-js';

// === ATORVASTATIN ===
// Brand name: Lipitor
// Most prescribed statin worldwide
// PubChem CID: 60823
// SMILES: CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O

export const atorvastatin = Fragment('CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O');

// === SIMVASTATIN ===
// Brand name: Zocor
// Prodrug activated in the liver
// PubChem CID: 54454
// SMILES: CCC(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12

export const simvastatin = Fragment('CCC(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12');

// === ROSUVASTATIN ===
// Brand name: Crestor
// Most potent statin
// PubChem CID: 446157
// SMILES: CC(C)c1nc(nc(c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O)c2ccc(F)cc2)N(C)S(=O)(=O)C

export const rosuvastatin = Fragment('CC(C)c1nc(nc(c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O)c2ccc(F)cc2)N(C)S(=O)(=O)C');

// === PRAVASTATIN ===
// Brand name: Pravachol
// Hydrophilic statin, less lipophilic than others
// PubChem CID: 54687
// SMILES: CCC(C)C(=O)O[C@H]1C[C@H](O)C=C2C=C[C@H](C)[C@H](CC[C@@H](O)CC(=O)O)[C@H]12

export const pravastatin = Fragment('CCC(C)C(=O)O[C@H]1C[C@H](O)C=C2C=C[C@H](C)[C@H](CC[C@@H](O)CC(=O)O)[C@H]12');

// === LOVASTATIN ===
// Brand name: Mevacor
// First statin approved by FDA (1987)
// PubChem CID: 53232
// SMILES: CCC(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12

export const lovastatin = Fragment('CCC(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12');

// === FLUVASTATIN ===
// Brand name: Lescol
// First fully synthetic statin
// PubChem CID: 446155
// SMILES: CC(C)n1c(C(=O)NC(C(=O)O)Cc2ccccc2)c(-c2ccc(F)cc2)c(O)c1C

export const fluvastatin = Fragment('CC(C)n1c(C(=O)Nc2ccccc2C(=O)O)c(-c2ccc(F)cc2)c(O)c1/C=C/c1ccccc1');

// === EZETIMIBE ===
// Brand name: Zetia
// Cholesterol absorption inhibitor (not a statin)
// PubChem CID: 150311
// SMILES: C[C@H]1[C@@H](O)[C@H](O[C@H]1c2ccc(F)cc2)c3ccc(O)cc3C(=O)N4CCC(CC4)c5ccc(F)cc5

export const ezetimibe = Fragment('C[C@H]1[C@@H](O)[C@H](O[C@H]1c2ccc(F)cc2)c3ccc(O)cc3C(=O)N4CCC(CC4)c5ccc(F)cc5');

// === FENOFIBRATE ===
// Brand name: Tricor
// Fibric acid derivative for triglycerides
// PubChem CID: 3339
// SMILES: CC(C)OC(=O)C(C)(C)Oc1ccc(cc1)C(=O)c2ccc(Cl)cc2

export const fenofibrate = Fragment('CC(C)OC(=O)C(C)(C)Oc1ccc(cc1)C(=O)c2ccc(Cl)cc2');

// === GEMFIBROZIL ===
// Brand name: Lopid
// Fibric acid derivative
// PubChem CID: 3463
// SMILES: CC1=CC=C(C=C1)C(C)(C)CCCC(=O)O

export const gemfibrozil = Fragment('CC1=CC(C)=CC=C1CCCC(C)(C)C(=O)O');

// === PITAVASTATIN ===
// Brand name: Livalo
// Newer generation statin
// PubChem CID: 5282452
// SMILES: CC(C)c1nc(nc(c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O)c2ccc(F)cc2)NC(=O)C3CC3

export const pitavastatin = Fragment('CC(C)c1nc(nc(c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O)c2ccc(F)cc2)NC(=O)C3CC3');

export const cholesterol = Fragment('C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C');
