import { describe, test, expect } from 'bun:test';
import { buildIntensityCurve, ppmToX, niceStep, plotSpectra } from './plot.js';

describe('niceStep', () => {
  test('returns 1 for raw ~1', () => {
    expect(niceStep(1.2)).toBe(1);
  });

  test('returns 2 for raw ~3', () => {
    expect(niceStep(2.8)).toBe(2);
  });

  test('returns 5 for raw ~6', () => {
    expect(niceStep(5.5)).toBe(5);
  });

  test('returns 10 for raw ~9', () => {
    expect(niceStep(8.5)).toBe(10);
  });

  test('returns 20 for raw ~25', () => {
    expect(niceStep(25)).toBe(20);
  });

  test('returns 50 for raw ~55', () => {
    expect(niceStep(55)).toBe(50);
  });
});

describe('ppmToX', () => {
  test('ppmMax maps to left edge', () => {
    expect(ppmToX(220, 0, 220, 800, 20)).toBe(20);
  });

  test('ppmMin maps to right edge', () => {
    expect(ppmToX(0, 0, 220, 800, 20)).toBe(820);
  });

  test('midpoint maps to center', () => {
    expect(ppmToX(110, 0, 220, 800, 20)).toBe(420);
  });
});

describe('buildIntensityCurve', () => {
  test('returns array of correct length', () => {
    const shifts = [{ shift: 100 }];
    const curve = buildIntensityCurve(shifts, 500, 0, 220, 0.5);
    expect(curve.length).toBe(500);
  });

  test('all values are non-negative', () => {
    const shifts = [{ shift: 50 }, { shift: 150 }];
    const curve = buildIntensityCurve(shifts, 200, 0, 220, 0.5);
    for (const v of curve) {
      expect(v).toBeGreaterThanOrEqual(0);
    }
  });

  test('peak is near the shift value', () => {
    const shifts = [{ shift: 110 }];
    const curve = buildIntensityCurve(shifts, 1000, 0, 220, 0.5);
    const maxIdx = curve.indexOf(Math.max(...curve));
    // Pixel 0 = 220 ppm, pixel 999 = 0 ppm. shift=110 → middle → ~idx 500
    expect(maxIdx).toBeGreaterThan(450);
    expect(maxIdx).toBeLessThan(550);
  });

  test('peak maximum equals 1.0 for a single Lorentzian at exact center', () => {
    // At exactly the shift position, Lorentzian = gamma^2 / (0 + gamma^2) = 1
    const shifts = [{ shift: 110 }];
    const curve = buildIntensityCurve(shifts, 1001, 0, 220, 0.5);
    const maxVal = Math.max(...curve);
    expect(maxVal).toBeCloseTo(1.0, 1);
  });

  test('two peaks produce two maxima', () => {
    const shifts = [{ shift: 50 }, { shift: 170 }];
    const curve = buildIntensityCurve(shifts, 1000, 0, 220, 0.5);
    // Find local maxima (value greater than both neighbors)
    const maxima = [];
    for (let i = 1; i < curve.length - 1; i++) {
      if (curve[i] > curve[i - 1] && curve[i] > curve[i + 1]) {
        maxima.push(i);
      }
    }
    expect(maxima.length).toBe(2);
  });

  test('empty shifts gives all zeros', () => {
    const curve = buildIntensityCurve([], 100, 0, 220, 0.5);
    for (const v of curve) {
      expect(v).toBe(0);
    }
  });

  test('wider lineWidth produces broader peaks', () => {
    const shifts = [{ shift: 110 }];
    const narrow = buildIntensityCurve(shifts, 1000, 0, 220, 0.1);
    const wide = buildIntensityCurve(shifts, 1000, 0, 220, 5.0);
    // Count pixels above 50% of max
    const countAboveHalf = (curve) => {
      const half = Math.max(...curve) / 2;
      return curve.filter((v) => v >= half).length;
    };
    expect(countAboveHalf(wide)).toBeGreaterThan(countAboveHalf(narrow));
  });
});

describe('plotSpectra', () => {
  test('works with mock canvas', () => {
    const calls = [];
    const mockCtx = {
      clearRect: (...args) => calls.push(['clearRect', ...args]),
      beginPath: () => calls.push(['beginPath']),
      moveTo: (...args) => calls.push(['moveTo', ...args]),
      lineTo: (...args) => calls.push(['lineTo', ...args]),
      closePath: () => calls.push(['closePath']),
      fill: () => calls.push(['fill']),
      stroke: () => calls.push(['stroke']),
      fillText: (...args) => calls.push(['fillText', ...args]),
      set fillStyle(v) {},
      set strokeStyle(v) {},
      set lineWidth(v) {},
      set font(v) {},
      set textAlign(v) {},
    };
    const mockCanvas = {
      width: 0,
      height: 0,
      getContext: () => mockCtx,
    };

    const shifts = [{ shift: 50 }, { shift: 130 }];
    plotSpectra(mockCanvas, shifts, { width: 400, height: 200 });

    expect(mockCanvas.width).toBe(400);
    expect(mockCanvas.height).toBe(200);
    expect(calls.some((c) => c[0] === 'clearRect')).toBe(true);
    expect(calls.some((c) => c[0] === 'fill')).toBe(true);
    expect(calls.some((c) => c[0] === 'stroke')).toBe(true);
    expect(calls.some((c) => c[0] === 'fillText')).toBe(true);
  });

  test('uses custom range', () => {
    const calls = [];
    const mockCtx = {
      clearRect: () => {},
      beginPath: () => {},
      moveTo: () => {},
      lineTo: () => {},
      closePath: () => {},
      fill: () => {},
      stroke: () => {},
      fillText: (...args) => calls.push(args),
      set fillStyle(v) {},
      set strokeStyle(v) {},
      set lineWidth(v) {},
      set font(v) {},
      set textAlign(v) {},
    };
    const mockCanvas = {
      width: 0,
      height: 0,
      getContext: () => mockCtx,
    };

    plotSpectra(mockCanvas, [{ shift: 5 }], { range: [0, 10], width: 400, height: 200 });

    // Tick labels should be within [0, 10]
    const labelValues = calls.map((c) => parseFloat(c[0])).filter((v) => !isNaN(v));
    for (const v of labelValues) {
      expect(v).toBeGreaterThanOrEqual(0);
      expect(v).toBeLessThanOrEqual(10);
    }
  });
});
