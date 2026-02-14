# vamos-hose

> **Pre-release / Work in Progress** -- This project is under active development. APIs, data formats, and file structures will change without notice.

## Goal

Given a SMILES string, return predicted 13C NMR chemical shifts by looking up HOSE codes in a preprocessed version of the NMRShiftDB2 database (1.7M entries, 30 elements).

The final product is a **pure JavaScript library** that runs entirely in the browser with no server-side dependencies. The database will be code-split into chunks so it can be lazy-loaded over the network.

### Target API

#### `lookupNmrShifts` -- SMILES to predicted shifts

```js
import { lookupNmrShifts } from 'vamos-hose';

const shifts = lookupNmrShifts('CCCCC1=NC(Cl)=C(CO)N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4', {
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

const candidates = estimateFromSpectra({
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

### To Do

- [ ] `lookupNmrShifts` -- JavaScript lookup engine (`src/lookup.js`)
- [ ] `plotSpectra` -- Canvas-based spectrum renderer (`src/plot.js`)
- [ ] `estimateFromSpectra` -- Reverse peak-to-structure search (`src/estimate.js`)
- [ ] SMILES-to-HOSE conversion in JS (openchemlib-js)
- [ ] Database code-splitting for lazy loading
- [ ] Integration tests for remaining drug categories
- [ ] Browser bundle / npm package

## Development

```bash
npm install
npm test
```

## License

Apache-2.0
