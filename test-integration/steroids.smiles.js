// Corticosteroids (Anti-Inflammatory Steroids)
// NOT to be confused with anabolic steroids
// Reduce inflammation by suppressing immune response

import { Fragment } from 'smiles-js';

// Note: Steroids have complex ring structures that are better represented as
// complete SMILES strings. The modular fragment approach from SELFIES doesn't
// translate well to JavaScript composition. These are provided as complete molecules

// === CORTISONE ===
// 11-keto glucocorticoid
// Trade names: Cortone Acetate
// Uses: Inflammation, adrenal insufficiency (prodrug)
// SMILES: CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C(C(=O)CO)CCC12

export const cortisone = Fragment('CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C(C(=O)CO)CCC12');

// === HYDROCORTISONE (CORTISOL) ===
// 11-hydroxy active form
// Trade names: Cortef, Hydrocortone
// Uses: Inflammation, topical anti-inflammatory
// SMILES: CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C(C(=O)CO)CCC12O

export const hydrocortisone = Fragment('CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C(C(=O)CO)CCC12O');
export const cortisol = hydrocortisone;

// === PREDNISONE ===
// 11-keto with delta-1
// Trade names: Deltasone, Rayos
// Uses: Inflammation, autoimmune conditions
// SMILES: CC12CC(=O)C=CC1=CC(O)C1C2CCC2(C)C(C(=O)CO)CCC12

export const prednisone = Fragment('CC12CC(=O)C=CC1=CC(O)C1C2CCC2(C)C(C(=O)CO)CCC12');

// === PREDNISOLONE ===
// 11-hydroxy with delta-1 (active form)
// Trade names: Prelone, Orapred
// Uses: Same as prednisone (hepatic activation not needed)
// SMILES: CC12CC(=O)C=CC1=CC(O)C1C2C(O)CC2(C)C(C(=O)CO)CCC12

export const prednisolone = Fragment('CC12CC(=O)C=CC1=CC(O)C1C2C(O)CC2(C)C(C(=O)CO)CCC12');

// === METHYLPREDNISOLONE ===
// Trade names: Medrol, Solu-Medrol
// C6-methyl substitution increases potency
// Uses: Severe inflammation, MS relapses
// SMILES: CC12CC(=O)C(C)=CC1=CC(O)C1C2C(O)CC2(C)C(C(=O)CO)CCC12

export const methylprednisolone = Fragment('CC12CC(=O)C(C)=CC1=CC(O)C1C2C(O)CC2(C)C(C(=O)CO)CCC12');

// === DEXAMETHASONE ===
// Trade names: Decadron, Dexasone
// 9-fluoro, 16-methyl substitutions (25x more potent than cortisol)
// Uses: Severe inflammation, COVID-19, cerebral edema
// SMILES: CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO

export const dexamethasone = Fragment('CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO');

// === BETAMETHASONE ===
// Trade names: Celestone
// Similar to dexamethasone (16-beta-methyl vs 16-alpha-methyl)
// Uses: Fetal lung maturation, inflammation
// SMILES: CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO

export const betamethasone = Fragment('CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO');

// === TRIAMCINOLONE ===
// Trade names: Kenalog, Aristospan
// 9-fluoro, 16-hydroxyl substitutions
// Uses: Joint injections, topical inflammation
// SMILES: CC12CC(O)C3C(CCC4=CC(=O)C=CC34C)C1(F)C(O)CC2(C)C(=O)CO

export const triamcinolone = Fragment('CC12CC(O)C3C(CCC4=CC(=O)C=CC34C)C1(F)C(O)CC2(C)C(=O)CO');

// === BUDESONIDE ===
// Trade names: Pulmicort, Entocort
// Topical with high first-pass metabolism (less systemic effects)
// Uses: Asthma, IBD, nasal allergies
// SMILES: CCCC1OC2CC3C4CCC5=CC(=O)C=CC5(C)C4(F)C(O)CC3(C)C2(O1)C(=O)CO

export const budesonide = Fragment('CCCC1OC2CC3C4CCC5=CC(=O)C=CC5(C)C4(F)C(O)CC3(C)C2(O1)C(=O)CO');

// === FLUTICASONE ===
// Trade names: Flovent, Flonase, Advair
// High topical potency, minimal systemic absorption
// Uses: Asthma, allergic rhinitis
// SMILES: CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=S)OCF

export const fluticasone = Fragment('CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=S)OCF');

// === BECLOMETHASONE ===
// Trade names: Qvar, Beconase
// Inhaled corticosteroid for asthma
// Uses: Asthma, allergic rhinitis
// SMILES: CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(Cl)C(O)CC2(C)C1(O)C(=O)CO

export const beclomethasone = Fragment('CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(Cl)C(O)CC2(C)C1(O)C(=O)CO');

// === FLUDROCORTISONE ===
// Trade names: Florinef
// 9-fluoro increases mineralocorticoid activity
// Uses: Adrenal insufficiency, orthostatic hypotension
// SMILES: CC12CCC(=O)C=C1CCC1C2(F)C(O)CC2(C)C(C(=O)CO)CCC12

export const fludrocortisone = Fragment('CC12CCC(=O)C=C1CCC1C2(F)C(O)CC2(C)C(C(=O)CO)CCC12');

// === MOMETASONE ===
// Trade names: Nasonex, Elocon
// Topical with chloro substitutions
// Uses: Nasal allergies, skin inflammation
// SMILES: CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(Cl)C(O)CC2(C)C1(O)C(=O)CCl

export const mometasone = Fragment('CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(Cl)C(O)CC2(C)C1(O)C(=O)CCl');
