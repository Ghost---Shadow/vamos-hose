import { normalize } from 'smiles-js';

const SUBSTITUTION_ATOMS = ['C', 'c', 'N', 'n', 'O', 'o', 'S', 's', 'F', 'Cl', 'Br'];

const HEAVY_ATOM_WEIGHT = {
  'B': 10.811, 'b': 10.811, 'C': 12.011, 'c': 12.011,
  'N': 14.007, 'n': 14.007, 'O': 15.999, 'o': 15.999,
  'F': 18.998, 'P': 30.974, 'p': 30.974, 'S': 32.065,
  's': 32.065, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904,
};

// Max valence for implicit H estimation
const MAX_VALENCE = {
  'B': 3, 'b': 2, 'C': 4, 'c': 3, 'N': 3, 'n': 2,
  'O': 2, 'o': 1, 'F': 1, 'P': 3, 'p': 2, 'S': 2,
  's': 1, 'Cl': 1, 'Br': 1, 'I': 1,
};
const ADDABLE_GROUPS = ['C', 'N', 'O'];

// Match atom tokens in SMILES: bracket atoms, two-letter organics, single-letter organics
const ATOM_RE = /\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnosp]/;

function tryNormalize(smiles) {
  try { return normalize(smiles); } catch { return null; }
}

function findAtomPositions(smiles) {
  const re = new RegExp(ATOM_RE.source, 'g');
  const positions = [];
  let m;
  while ((m = re.exec(smiles)) !== null) {
    positions.push({ start: m.index, end: m.index + m[0].length, atom: m[0] });
  }
  return positions;
}

function dedup(results, canonical) {
  const seen = new Set([canonical]);
  return results.filter(r => {
    if (seen.has(r.smiles)) return false;
    seen.add(r.smiles);
    return true;
  });
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
 * Accurate enough for MW-based candidate filtering (~±5 Da for typical organics).
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

  // Count explicit bonds: each bond in SMILES connects two atoms
  // Bonds = (number of atoms - 1) for a tree, plus ring closures
  const ringClosures = (smiles.match(/[0-9%]/g) || []).length / 2;
  const explicitBonds = (atoms.length - 1) + ringClosures;

  // Double/triple bonds consume extra valence
  const doubleBonds = (smiles.match(/=/g) || []).length;
  const tripleBonds = (smiles.match(/#/g) || []).length;
  const extraBondOrders = doubleBonds + tripleBonds * 2;

  // Each bond uses 1 valence from each endpoint = 2 valence units total
  // Implicit H fills remaining valence
  const usedValence = (explicitBonds + extraBondOrders) * 2;
  const implicitH = Math.max(0, totalMaxValence - usedValence);

  return heavyWeight + implicitH * 1.008;
}

export function enumerateSubstitutions(smiles) {
  const canonical = tryNormalize(smiles);
  if (!canonical) return [];
  const results = [];
  const positions = findAtomPositions(canonical);

  for (const pos of positions) {
    for (const newAtom of SUBSTITUTION_ATOMS) {
      if (newAtom === pos.atom) continue;
      const s = canonical.slice(0, pos.start) + newAtom + canonical.slice(pos.end);
      const norm = tryNormalize(s);
      if (norm && norm !== canonical) {
        results.push({ smiles: norm, description: `sub ${pos.atom}->${newAtom}` });
      }
    }
  }

  return dedup(results, canonical);
}

export function enumerateAdditions(smiles) {
  const canonical = tryNormalize(smiles);
  if (!canonical) return [];
  const results = [];
  const positions = findAtomPositions(canonical);

  for (const pos of positions) {
    for (const atom of ADDABLE_GROUPS) {
      const s = canonical.slice(0, pos.end) + `(${atom})` + canonical.slice(pos.end);
      const norm = tryNormalize(s);
      if (norm && norm !== canonical) {
        results.push({ smiles: norm, description: `add (${atom})` });
      }
    }
  }

  return dedup(results, canonical);
}

export function enumerateBondChanges(smiles) {
  const canonical = tryNormalize(smiles);
  if (!canonical) return [];
  const results = [];

  // Tokenize to find bond positions
  const tokenRe = /\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnosp]|[=#]|%\d{2}|\d|[().\/\\:\-]/g;
  const tokens = [];
  let m;
  while ((m = tokenRe.exec(canonical)) !== null) {
    tokens.push({ text: m[0], index: m.index, end: m.index + m[0].length });
  }

  for (let i = 0; i < tokens.length; i++) {
    const t = tokens[i];

    // Change or remove existing explicit bond
    if (t.text === '=' || t.text === '#') {
      for (const alt of ['=', '#', ''].filter(b => b !== t.text)) {
        const s = canonical.slice(0, t.index) + alt + canonical.slice(t.end);
        const norm = tryNormalize(s);
        if (norm && norm !== canonical) {
          results.push({ smiles: norm, description: `bond ${t.text}->${alt || 'single'}` });
        }
      }
    }

    // Insert bond between adjacent atom tokens
    if (ATOM_RE.test(t.text) && i + 1 < tokens.length && ATOM_RE.test(tokens[i + 1].text)) {
      for (const bond of ['=', '#']) {
        const s = canonical.slice(0, t.end) + bond + canonical.slice(t.end);
        const norm = tryNormalize(s);
        if (norm && norm !== canonical) {
          results.push({ smiles: norm, description: `bond single->${bond}` });
        }
      }
    }
  }

  return dedup(results, canonical);
}

export function enumerateRemovals(smiles) {
  const canonical = tryNormalize(smiles);
  if (!canonical) return [];
  const results = [];

  // Remove parenthesized branches
  for (const match of canonical.matchAll(/\([^()]+\)/g)) {
    const s = canonical.slice(0, match.index) + canonical.slice(match.index + match[0].length);
    const norm = tryNormalize(s);
    if (norm && norm !== canonical) {
      results.push({ smiles: norm, description: `remove ${match[0]}` });
    }
  }

  // Remove terminal atom (last atom if at end of string)
  const positions = findAtomPositions(canonical);
  if (positions.length > 1) {
    const last = positions[positions.length - 1];
    if (last.end === canonical.length) {
      const s = canonical.slice(0, last.start);
      const norm = tryNormalize(s);
      if (norm && norm !== canonical) {
        results.push({ smiles: norm, description: 'remove last atom' });
      }
    }
  }

  return dedup(results, canonical);
}

export function enumerateFragmentAttachments(smiles, fragments) {
  const canonical = tryNormalize(smiles);
  if (!canonical || !fragments || fragments.length === 0) return [];
  const results = [];
  const positions = findAtomPositions(canonical);

  for (const frag of fragments) {
    for (const pos of positions) {
      const s = canonical.slice(0, pos.end) + `(${frag})` + canonical.slice(pos.end);
      const norm = tryNormalize(s);
      if (norm && norm !== canonical) {
        results.push({ smiles: norm, description: `attach (${frag})` });
      }
    }
  }

  return dedup(results, canonical);
}

export function enumerateAllMutations(smiles, hoseFragments = []) {
  const canonical = tryNormalize(smiles);
  if (!canonical) return [];

  const all = [
    ...enumerateSubstitutions(canonical),
    ...enumerateAdditions(canonical),
    ...enumerateBondChanges(canonical),
    ...enumerateRemovals(canonical),
    ...enumerateFragmentAttachments(canonical, hoseFragments),
  ];

  return dedup(all, canonical);
}
