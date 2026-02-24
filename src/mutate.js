import { normalize } from 'smiles-js';

const SUBSTITUTION_ATOMS = ['C', 'c', 'N', 'n', 'O', 'o', 'S', 's', 'F', 'Cl', 'Br'];
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
