# vamos-hose

> **Pre-release / Work in Progress** -- This project is under active development. APIs, data formats, and file structures will change without notice.

## Goal

Given a SMILES string, return predicted 13C NMR chemical shifts by looking up HOSE codes in a preprocessed version of the NMRShiftDB2 database (1.7M entries, 30 elements).

The final product is a **pure JavaScript library** that runs entirely in the browser with no server-side dependencies. The database will be code-split into chunks so it can be lazy-loaded over the network.

### Target API

#### `lookupNmrShifts` -- SMILES to predicted shifts

```js
import { lookupNmrShifts } from 'vamos-hose';

const shifts = await lookupNmrShifts('CCCCC1=NC(Cl)=C(CO)N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4', {
  nucleus: '13C',
});
// Returns: [{ shift: 13.9, atom: 'C', hose: '...', smiles: '...' }, ...]
```

#### `plotSpectra` -- Render a synthetic NMR spectrum

```js
import { plotSpectra } from 'vamos-hose';

const canvas = document.getElementById('nmr-canvas');
plotSpectra(canvas, shifts, {
  range: [0, 220],    // ppm range
  width: 800,
  height: 300,
  lineWidth: 0.5,     // Lorentzian half-width in ppm
});
```

#### `estimateFromSpectra` -- Reverse lookup: peaks to candidate structures

```js
import { estimateFromSpectra } from 'vamos-hose';

const candidates = await estimateFromSpectra({
  nucleus: '13C',
  peaks: [14.1, 22.7, 32.0, 127.5, 128.3, 130.6, 137.0, 174.1],
  tolerance: 2.0,     // ppm tolerance per peak
  minMatches: 5,      // minimum peaks that must match
});
// Returns: [{ smiles: '...', hose: '...', matchedPeaks: 7, score: 0.87 }, ...]
```

## Current Status

### Done

- [x] HOSE-to-SMILES converter (Python) -- 100% conversion rate on 1.7M entries
- [x] Database preprocessor -- HOSE-keyed JSON with SMILES as secondary field
- [x] Full database build -- 1,421,845 unique HOSE keys, 210 MB JSON
- [x] Integration test suite for hypertension medications (losartan, valsartan, irbesartan, telmisartan)

- [x] `lookupNmrShifts` -- JavaScript lookup engine (`src/lookup.js`)
- [x] `plotSpectra` -- Canvas-based spectrum renderer (`src/plot.js`)
- [x] `estimateFromSpectra` -- Reverse peak-to-structure search (`src/estimate.js`)
- [x] SMILES-to-HOSE conversion in JS (openchemlib-js)
- [x] Database code-splitting for lazy loading (256 hash-based chunks, ~0.8 MB each)
- [x] Browser-compatible module (import.meta.url chunk resolution, npm package fields)

- [x] Integration tests for all 10 drug categories (276 tests across 35+ molecules)

### To Do

- [ ] npm publish

## Development

```bash
npm install
npm test
```

## Acknowledgments

This project builds on the work of several open-source projects and academic publications:

- **[nmrshiftdb2](https://nmrshiftdb.nmr.uni-koeln.de/)** — NMR chemical shift database (Stefan Kuhn, Christoph Steinbeck et al.). The HOSE code generator in `src/hose-generator.js` is a JavaScript port of the [ExtendedHOSECodeGenerator](https://sourceforge.net/p/nmrshiftdb2/code/HEAD/tree/trunk/nmrshiftdb2/src/java/org/openscience/nmrshiftdb/util/ExtendedHOSECodeGenerator.java) (AGPL v3).
- **[CDK (Chemistry Development Kit)](https://cdk.github.io/)** — The canonical labeling algorithm is ported from CDK's [CanonicalLabeler](https://github.com/cdk/cdk/blob/main/base/standard/src/main/java/org/openscience/cdk/graph/invariant/CanonicalLabeler.java) (Oliver Horlacher, LGPL 2.1), which implements the Weininger algorithm (D. Weininger et al., *J. Chem. Inf. Comput. Sci.*, 1989, 29, 97-101).
- **[openchemlib-js](https://github.com/cheminfo/openchemlib-js)** — SMILES parsing and molecule representation (BSD-3-Clause).
- **W. Bremser** — Original HOSE code concept (*Analytica Chimica Acta*, 1978, 103, 355-365).

## License

AGPL-3.0-or-later

The HOSE code generator (`src/hose-generator.js`) is a JavaScript port of nmrshiftdb2's [ExtendedHOSECodeGenerator](https://sourceforge.net/p/nmrshiftdb2/code/HEAD/tree/trunk/nmrshiftdb2/src/java/org/openscience/nmrshiftdb/util/ExtendedHOSECodeGenerator.java), which is licensed under AGPL v3. It also ports CDK's [CanonicalLabeler](https://github.com/cdk/cdk/blob/main/base/standard/src/main/java/org/openscience/cdk/graph/invariant/CanonicalLabeler.java), licensed under LGPL 2.1. Since the algorithms were directly ported (not linked as a library), the derived code is a derivative work and the AGPL v3 copyleft applies to the entire project. See `LICENSE` and the `src/hose-generator.js` file header for details.

Note: The AGPL applies to the **software** only. Chemical shifts, NMR predictions, drug candidates, or any other scientific results produced by running this library are not covered by the AGPL and remain the property of the user.
