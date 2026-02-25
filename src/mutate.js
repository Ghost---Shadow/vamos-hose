import OCL from 'openchemlib';

function normalize(smiles) {
  return OCL.Molecule.fromSmiles(smiles).toSmiles();
}

// ── Constants ──────────────────────────────────────────────────────────────────

const SUBSTITUTION_ATOMS = ['C', 'c', 'N', 'n', 'O', 'o', 'S', 's', 'F', 'Cl', 'Br'];
const ADDABLE_GROUPS = ['C', 'N', 'O'];

const HEAVY_ATOM_WEIGHT = {
  'B': 10.811, 'b': 10.811, 'C': 12.011, 'c': 12.011,
  'N': 14.007, 'n': 14.007, 'O': 15.999, 'o': 15.999,
  'F': 18.998, 'P': 30.974, 'p': 30.974, 'S': 32.065,
  's': 32.065, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904,
};

const MAX_VALENCE = {
  'B': 3, 'b': 2, 'C': 4, 'c': 3, 'N': 3, 'n': 2,
  'O': 2, 'o': 1, 'F': 1, 'P': 3, 'p': 2, 'S': 2,
  's': 1, 'Cl': 1, 'Br': 1, 'I': 1,
};

// Match atom tokens: bracket atoms, two-letter organics, single-letter organics
const ATOM_RE = /\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnosp]/;

// Full SMILES tokenizer: atoms, bonds, ring digits, parens, etc.
const TOKEN_RE = /\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnosp]|[=#]|%\d{2}|\d|[().\/\\:\-]/g;

// ── Helpers ────────────────────────────────────────────────────────────────────

function tryNormalize(smiles) {
  try { return normalize(smiles); } catch { return null; }
}

/** Fast validation: balanced parens + non-empty. Invalid SMILES caught later by safePredictShifts. */
function fastValidate(smiles) {
  if (!smiles || smiles.length === 0) return null;
  let depth = 0;
  for (let i = 0; i < smiles.length; i++) {
    if (smiles[i] === '(') depth++;
    else if (smiles[i] === ')') { depth--; if (depth < 0) return null; }
  }
  return depth === 0 ? smiles : null;
}

/** Find all atom positions in a SMILES string. */
function findAtomPositions(smiles) {
  const re = new RegExp(ATOM_RE.source, 'g');
  const positions = [];
  let m;
  while ((m = re.exec(smiles)) !== null) {
    positions.push({ start: m.index, end: m.index + m[0].length, atom: m[0] });
  }
  return positions;
}

/** Tokenize a SMILES string into atoms, bonds, parens, digits. */
function tokenize(smiles) {
  const re = new RegExp(TOKEN_RE.source, 'g');
  const tokens = [];
  let m;
  while ((m = re.exec(smiles)) !== null) {
    tokens.push({ text: m[0], index: m.index, end: m.index + m[0].length });
  }
  return tokens;
}

/** Remove duplicate candidates (by SMILES string). */
function dedup(results, canonical) {
  const seen = new Set([canonical]);
  return results.filter(r => {
    if (seen.has(r.smiles)) return false;
    seen.add(r.smiles);
    return true;
  });
}

// ── Factory: normalize → generate → dedup ──────────────────────────────────────

/**
 * Wraps an enumeration generator with normalize-on-input and dedup-on-output.
 * generator(canonical) returns raw candidate array.
 */
function withNormDedup(smiles, generator) {
  const canonical = tryNormalize(smiles);
  if (!canonical) return [];
  return dedup(generator(canonical), canonical);
}

/**
 * Create a validated candidate if it passes fastValidate and differs from canonical.
 * Returns null if invalid or identical.
 */
function candidate(smiles, description, canonical) {
  const norm = fastValidate(smiles);
  if (norm && norm !== canonical) return { smiles: norm, description };
  return null;
}

// ── SMILES analysis utilities ──────────────────────────────────────────────────

/** Extract the element symbol from an atom token, normalizing aromatic to uppercase. */
function elementFromToken(token) {
  if (token.startsWith('[')) {
    const m = token.match(/\[([A-Za-z][a-z]?)/);
    if (!m) return null;
    const e = m[1];
    return e[0].toUpperCase() + e.slice(1);
  }
  return token[0].toUpperCase() + token.slice(1);
}

/** Count atoms by element in a SMILES string. Aromatic atoms normalized to uppercase. */
export function countAtomsByElement(smiles) {
  const counts = {};
  const re = new RegExp(ATOM_RE.source, 'g');
  let m;
  while ((m = re.exec(smiles)) !== null) {
    const elem = elementFromToken(m[0]);
    if (elem) counts[elem] = (counts[elem] || 0) + 1;
  }
  return counts;
}

/** Fast carbon count from SMILES (no parsing needed). */
export function countCarbons(smiles) {
  let count = 0;
  for (let i = 0; i < smiles.length; i++) {
    const ch = smiles[i];
    if (ch === 'C' && smiles[i + 1] !== 'l') count++;
    else if (ch === 'c') count++;
  }
  return count;
}

export function containsOnlyAtoms(smiles, allowedSet) {
  const re = new RegExp(ATOM_RE.source, 'g');
  let m;
  while ((m = re.exec(smiles)) !== null) {
    if (!allowedSet.has(m[0])) return false;
  }
  return true;
}

export function validateSmiles(smiles) {
  try {
    const canonical = normalize(smiles);
    return { valid: true, canonical };
  } catch {
    return { valid: false, canonical: null };
  }
}

export function getAtomCount(smiles) {
  return findAtomPositions(smiles).length;
}

/**
 * Fast molecular weight estimate from SMILES string.
 * Sums heavy atom weights + estimates implicit hydrogens from valence rules.
 */
export function getMolecularWeight(smiles) {
  const atoms = findAtomPositions(smiles);
  if (atoms.length === 0) return 0;

  let heavyWeight = 0;
  let totalMaxValence = 0;
  for (const { atom } of atoms) {
    heavyWeight += HEAVY_ATOM_WEIGHT[atom] || 12.011;
    totalMaxValence += MAX_VALENCE[atom] || 4;
  }

  const ringClosures = (smiles.match(/[0-9%]/g) || []).length / 2;
  const explicitBonds = (atoms.length - 1) + ringClosures;
  const doubleBonds = (smiles.match(/=/g) || []).length;
  const tripleBonds = (smiles.match(/#/g) || []).length;
  const usedValence = (explicitBonds + doubleBonds + tripleBonds * 2) * 2;
  const implicitH = Math.max(0, totalMaxValence - usedValence);

  return heavyWeight + implicitH * 1.008;
}

// ── Mutation enumerators ───────────────────────────────────────────────────────

/** Replace each atom with every substitution atom. */
export function enumerateSubstitutions(smiles) {
  return withNormDedup(smiles, canonical => {
    const results = [];
    for (const pos of findAtomPositions(canonical)) {
      for (const newAtom of SUBSTITUTION_ATOMS) {
        if (newAtom === pos.atom) continue;
        const c = candidate(
          canonical.slice(0, pos.start) + newAtom + canonical.slice(pos.end),
          `sub ${pos.atom}->${newAtom}`, canonical
        );
        if (c) results.push(c);
      }
    }
    return results;
  });
}

/** Add a branch (atom) after each atom position. */
export function enumerateAdditions(smiles) {
  return withNormDedup(smiles, canonical => {
    const results = [];
    for (const pos of findAtomPositions(canonical)) {
      for (const atom of ADDABLE_GROUPS) {
        const c = candidate(
          canonical.slice(0, pos.end) + `(${atom})` + canonical.slice(pos.end),
          `add (${atom})`, canonical
        );
        if (c) results.push(c);
      }
    }
    return results;
  });
}

/** Modify existing bonds or insert new ones between adjacent atoms. */
export function enumerateBondChanges(smiles) {
  return withNormDedup(smiles, canonical => {
    const results = [];
    const tokens = tokenize(canonical);

    for (let i = 0; i < tokens.length; i++) {
      const t = tokens[i];

      // Change or remove existing explicit bond
      if (t.text === '=' || t.text === '#') {
        for (const alt of ['=', '#', ''].filter(b => b !== t.text)) {
          const c = candidate(
            canonical.slice(0, t.index) + alt + canonical.slice(t.end),
            `bond ${t.text}->${alt || 'single'}`, canonical
          );
          if (c) results.push(c);
        }
      }

      // Insert bond between adjacent atom tokens
      if (ATOM_RE.test(t.text) && i + 1 < tokens.length && ATOM_RE.test(tokens[i + 1].text)) {
        for (const bond of ['=', '#']) {
          const c = candidate(
            canonical.slice(0, t.end) + bond + canonical.slice(t.end),
            `bond single->${bond}`, canonical
          );
          if (c) results.push(c);
        }
      }
    }
    return results;
  });
}

/** Remove parenthesized branches and terminal atoms. */
export function enumerateRemovals(smiles) {
  return withNormDedup(smiles, canonical => {
    const results = [];

    // Remove parenthesized branches
    for (const match of canonical.matchAll(/\([^()]+\)/g)) {
      const c = candidate(
        canonical.slice(0, match.index) + canonical.slice(match.index + match[0].length),
        `remove ${match[0]}`, canonical
      );
      if (c) results.push(c);
    }

    // Remove terminal atom
    const positions = findAtomPositions(canonical);
    if (positions.length > 1) {
      const last = positions[positions.length - 1];
      if (last.end === canonical.length) {
        const c = candidate(canonical.slice(0, last.start), 'remove last atom', canonical);
        if (c) results.push(c);
      }
    }

    return results;
  });
}

/**
 * Insert atoms inline at key token boundaries.
 * Covers: after ring-close digits, between atom and ')', after ')', and at end.
 */
export function enumerateInlineInsertions(smiles) {
  return withNormDedup(smiles, canonical => {
    const results = [];
    const tokens = tokenize(canonical);

    for (let i = 0; i < tokens.length; i++) {
      const t = tokens[i];
      const next = tokens[i + 1];
      if (!next) continue;

      const shouldInsert =
        // After ring closure digit, before ')'
        (/^\d$|^%\d{2}$/.test(t.text) && next.text === ')') ||
        // Between atom and ')'
        (ATOM_RE.test(t.text) && next.text === ')') ||
        // After ')'
        (t.text === ')');

      if (shouldInsert) {
        for (const atom of ADDABLE_GROUPS) {
          const c = candidate(
            canonical.slice(0, t.end) + atom + canonical.slice(t.end),
            `inline ${atom}`, canonical
          );
          if (c) results.push(c);
        }
      }
    }

    // Append atom at end of SMILES
    for (const atom of ADDABLE_GROUPS) {
      const c = candidate(canonical + atom, `inline ${atom}`, canonical);
      if (c) results.push(c);
    }

    return results;
  });
}

/** Form new rings by inserting ring closure digits between atom pairs 3-7 apart. */
export function enumerateRingClosures(smiles) {
  return withNormDedup(smiles, canonical => {
    const positions = findAtomPositions(canonical);
    const nextDigit = findNextRingDigit(canonical);
    if (!nextDigit) return [];

    const origAtomCount = positions.length;
    const results = [];

    for (let i = 0; i < positions.length; i++) {
      for (let j = i + 3; j < Math.min(i + 8, positions.length); j++) {
        const posA = positions[i];
        const posB = positions[j];
        const s = canonical.slice(0, posA.end) + nextDigit +
                  canonical.slice(posA.end, posB.end) + nextDigit +
                  canonical.slice(posB.end);
        const norm = fastValidate(s);
        if (!norm || norm === canonical) continue;
        if (findAtomPositions(norm).length !== origAtomCount) continue;
        results.push({ smiles: norm, description: `ring close ${nextDigit}` });
      }
    }
    return results;
  });
}

/** Attach fragments as branches at each atom position. */
export function enumerateFragmentAttachments(smiles, fragments) {
  if (!fragments || fragments.length === 0) return [];
  return withNormDedup(smiles, canonical => {
    const results = [];
    const positions = findAtomPositions(canonical);

    for (const frag of fragments) {
      for (const pos of positions) {
        const c = candidate(
          canonical.slice(0, pos.end) + `(${frag})` + canonical.slice(pos.end),
          `attach (${frag})`, canonical
        );
        if (c) results.push(c);
      }
    }
    return results;
  });
}

/** Enumerate all mutation types at once. */
export function enumerateAllMutations(smiles, hoseFragments = []) {
  const canonical = tryNormalize(smiles);
  if (!canonical) return [];

  const all = [
    ...enumerateSubstitutions(canonical),
    ...enumerateAdditions(canonical),
    ...enumerateInlineInsertions(canonical),
    ...enumerateBondChanges(canonical),
    ...enumerateRemovals(canonical),
    ...enumerateRingClosures(canonical),
    ...enumerateFragmentAttachments(canonical, hoseFragments),
  ];

  return dedup(all, canonical);
}

// ── Internal helpers ───────────────────────────────────────────────────────────

function findNextRingDigit(smiles) {
  const usedDigits = new Set();
  for (const ch of smiles) {
    if (ch >= '0' && ch <= '9') usedDigits.add(ch);
  }
  for (let d = 1; d <= 9; d++) {
    if (!usedDigits.has(String(d))) return String(d);
  }
  return null;
}
