# Plan: HOSE → SMILES + Molecule Rendering for Reverse Tab

## Context
Reverse lookup returns opaque HOSE code strings. User wants each peak to show a **molecule diagram** so they can visually identify the chemical environment. HOSE codes → SMILES via existing `queryHose()` → `entry.s`. SmilesDrawer already loaded. Develop on localhost, push gh-pages only at the end.

## 3 file changes, no new files

### 1. `src/estimate.js` — Add `resolveHoseSmiles()` [Sonnet]
- Add `queryHose` to existing import from `./database.js`
- New exported function:
  ```
  resolveHoseSmiles(results) → { [peak]: { hose, shift, error, smiles: string|null } }
  ```
- Takes output of `estimateFromSpectra`, extracts top-1 per peak
- Calls `queryHose(hose)` in parallel → gets `{ smiles }` from chunks
- `smiles: null` if chunk lookup fails (graceful fallback)

### 2. `index.js` — Export new function [Haiku]
- Add `resolveHoseSmiles` to line 3 re-export

### 3. `index.html` (gh-pages) — Rewrite reverse tab UI [Opus]
- Import `resolveHoseSmiles` from `./index.js`
- Remove "Results per peak" input (always 1 now)
- Add `.peak-card` CSS (card per peak with molecule SVG)
- Rewrite `predictReverse()`:
  1. `estimateFromSpectra({ peaks, tolerance, maxResults: 1 })`
  2. `resolveHoseSmiles(raw)` — loads ~7MB of chunks for 9 peaks
  3. Render card-per-peak: peak label + molecule SVG + HOSE + SMILES
  4. `drawMolecule(smiles, svgId)` for each (reuse existing fn)
- Two-phase spinner: "Searching…" → "Loading structures…"

### Key existing functions (no changes)
- `queryHose(hoseCode)` → `{ smiles, avgShift, ... }` (database.js, already exported)
- `drawMolecule(smiles, svgId)` (index.html, already renders SVG via SmilesDrawer)

## Verification
1. `bun test` on main (estimate.js changes)
2. Copy changed files to gh-pages, test on localhost:8765
3. Push gh-pages only after Chrome verification
