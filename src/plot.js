/**
 * plotSpectra — Render a synthetic NMR spectrum on an HTML Canvas.
 *
 * Draws Lorentzian peaks at each predicted chemical shift, with
 * x-axis inverted (high ppm left, low ppm right) per NMR convention.
 *
 * @param {HTMLCanvasElement} canvas - target canvas element
 * @param {Array<{shift: number}>} shifts - array from lookupNmrShifts
 * @param {object} [options]
 * @param {[number, number]} [options.range=[0, 220]] - ppm range [min, max]
 * @param {number} [options.width=800]  - canvas width in px
 * @param {number} [options.height=300] - canvas height in px
 * @param {number} [options.lineWidth=0.5] - Lorentzian half-width at half-max (ppm)
 */
export function plotSpectra(canvas, shifts, options = {}) {
  const {
    range = [0, 220],
    width = 800,
    height = 300,
    lineWidth = 0.5,
  } = options;

  canvas.width = width;
  canvas.height = height;

  const ctx = canvas.getContext('2d');

  const margin = { top: 10, right: 20, bottom: 40, left: 20 };
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;

  const [ppmMin, ppmMax] = range;

  // Build the intensity curve: one value per pixel column
  const intensities = buildIntensityCurve(shifts, plotW, ppmMin, ppmMax, lineWidth);

  // Find max intensity for normalisation
  const maxIntensity = Math.max(...intensities, 1e-12);

  // --- draw ---
  ctx.clearRect(0, 0, width, height);

  // Spectrum fill
  ctx.beginPath();
  ctx.moveTo(margin.left, margin.top + plotH);
  for (let i = 0; i < plotW; i++) {
    const y = margin.top + plotH - (intensities[i] / maxIntensity) * plotH;
    ctx.lineTo(margin.left + i, y);
  }
  ctx.lineTo(margin.left + plotW, margin.top + plotH);
  ctx.closePath();
  ctx.fillStyle = 'rgba(0, 100, 200, 0.15)';
  ctx.fill();

  // Spectrum line
  ctx.beginPath();
  for (let i = 0; i < plotW; i++) {
    const y = margin.top + plotH - (intensities[i] / maxIntensity) * plotH;
    if (i === 0) ctx.moveTo(margin.left + i, y);
    else ctx.lineTo(margin.left + i, y);
  }
  ctx.strokeStyle = 'rgba(0, 100, 200, 1)';
  ctx.lineWidth = 1;
  ctx.stroke();

  // Axes
  drawAxes(ctx, margin, plotW, plotH, ppmMin, ppmMax);
}

/**
 * Build an array of summed Lorentzian intensities, one per pixel column.
 *
 * Exported for unit-testing.
 *
 * @param {Array<{shift: number}>} shifts
 * @param {number} numPoints - number of pixel columns
 * @param {number} ppmMin
 * @param {number} ppmMax
 * @param {number} gamma - half-width at half-max (ppm)
 * @returns {number[]} intensity array
 */
export function buildIntensityCurve(shifts, numPoints, ppmMin, ppmMax, gamma) {
  const intensities = new Array(numPoints).fill(0);
  const g2 = gamma * gamma;

  for (let i = 0; i < numPoints; i++) {
    // NMR convention: leftmost pixel = ppmMax, rightmost = ppmMin
    const ppm = ppmMax - (i / (numPoints - 1)) * (ppmMax - ppmMin);

    for (const s of shifts) {
      const diff = ppm - s.shift;
      intensities[i] += g2 / (diff * diff + g2);
    }
  }
  return intensities;
}

/**
 * Convert a ppm value to an x pixel coordinate (NMR convention: inverted).
 *
 * Exported for unit-testing.
 */
export function ppmToX(ppm, ppmMin, ppmMax, plotW, marginLeft) {
  const frac = (ppmMax - ppm) / (ppmMax - ppmMin);
  return marginLeft + frac * plotW;
}

/**
 * Draw x-axis baseline, tick marks, and ppm labels.
 */
function drawAxes(ctx, margin, plotW, plotH, ppmMin, ppmMax) {
  const baseY = margin.top + plotH;

  // Baseline
  ctx.beginPath();
  ctx.moveTo(margin.left, baseY);
  ctx.lineTo(margin.left + plotW, baseY);
  ctx.strokeStyle = '#333';
  ctx.lineWidth = 1;
  ctx.stroke();

  // Tick marks & labels — nice round intervals
  const span = ppmMax - ppmMin;
  const rawStep = span / 8;
  const step = niceStep(rawStep);

  ctx.fillStyle = '#333';
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'center';

  const firstTick = Math.ceil(ppmMin / step) * step;
  for (let ppm = firstTick; ppm <= ppmMax; ppm += step) {
    const x = ppmToX(ppm, ppmMin, ppmMax, plotW, margin.left);
    ctx.beginPath();
    ctx.moveTo(x, baseY);
    ctx.lineTo(x, baseY + 5);
    ctx.stroke();
    ctx.fillText(ppm.toFixed(0), x, baseY + 18);
  }

  // Axis label
  ctx.font = '12px sans-serif';
  ctx.fillText('ppm', margin.left + plotW / 2, baseY + 34);
}

/**
 * Round a raw step size to a "nice" number (1, 2, 5, 10, 20, 50, …).
 *
 * Exported for unit-testing.
 */
export function niceStep(raw) {
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  if (norm <= 1.5) return mag;
  if (norm <= 3.5) return 2 * mag;
  if (norm <= 7.5) return 5 * mag;
  return 10 * mag;
}
