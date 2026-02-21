// Opioid Analgesics
// Strong prescription painkillers acting on opioid receptors

// === MORPHINE ===
// Natural opioid alkaloid - base structure for many opioids
export const morphine = 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O';

// === CODEINE ===
// Morphine with C3 phenolic OH methylated to methoxy
export const codeine = 'CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O';

// === OXYCODONE ===
// Semi-synthetic: C6 ketone, C14 hydroxyl, C3 methoxy
export const oxycodone = 'CN1CCC23C4C(=O)CCC2(C1CC5=C3C(=C(C=C5)OC)O4)O';

// === HYDROCODONE ===
// Semi-synthetic: C6 ketone, C3 methoxy
export const hydrocodone = 'CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4';

// === FENTANYL ===
// Synthetic phenylpiperidine (NOT morphinan structure)
export const fentanyl = 'CCC(=O)N(C1CCN(CC1)CCC2=CC=CC=C2)C3=CC=CC=C3';

// === TRAMADOL ===
// Synthetic centrally acting analgesic (weak opioid + SNRI activity)
export const tramadol = 'CN(C)CC1CCCCC1(C2=CC(=CC=C2)OC)O';

// === METHADONE ===
// Synthetic acyclic opioid
export const methadone = 'CCC(=O)C(CC(C)N(C)C)(C1=CC=CC=C1)C2=CC=CC=C2';
