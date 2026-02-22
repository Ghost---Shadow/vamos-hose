import { describe, test, expect } from 'bun:test';
import { statSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';
import { MinHeap } from './min-heap.js';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = join(__dirname, '..');
const INVERTED_DIR = join(ROOT, 'inverted-index');

/**
 * Simulate estimateFromSpectra locally, tracking every file access and byte count.
 * New format: PPM files contain [hoseCode, shift] — no chunk loading needed.
 */
async function diagnoseQuery({ peaks, tolerance = 2.0, maxResults = 10 }) {
  const stats = {
    ppmFilesLoaded: new Set(),
    ppmBytesTotal: 0,
    totalEntries: 0,
    entriesWithinTolerance: 0,
    peakDetails: {},
  };

  for (const peak of peaks) {
    const lo = Math.max(0, Math.floor(peak - tolerance));
    const hi = Math.min(249, Math.floor(peak + tolerance));

    const peakStat = {
      ppmRange: [lo, hi],
      ppmFilesNeeded: hi - lo + 1,
      entriesScanned: 0,
      entriesWithinTolerance: 0,
      topResults: [],
    };

    const allEntries = [];
    for (let p = lo; p <= hi; p++) {
      const fileName = `ppm_${String(p).padStart(3, '0')}.js`;
      const filePath = join(INVERTED_DIR, fileName);
      if (!stats.ppmFilesLoaded.has(p)) {
        stats.ppmFilesLoaded.add(p);
        stats.ppmBytesTotal += statSync(filePath).size;
      }
      const mod = await import(filePath);
      allEntries.push(...mod.default);
    }

    peakStat.entriesScanned = allEntries.length;
    stats.totalEntries += allEntries.length;

    const heap = new MinHeap(maxResults);
    for (const [hoseCode, shift] of allEntries) {
      const error = Math.abs(shift - peak);
      if (error <= tolerance) {
        heap.push({ key: error, hose: hoseCode, shift });
        peakStat.entriesWithinTolerance++;
        stats.entriesWithinTolerance++;
      }
    }

    const topEntries = heap.toArray();
    for (const entry of topEntries) {
      peakStat.topResults.push({
        hose: entry.hose,
        shift: entry.shift,
        error: entry.key,
      });
    }

    stats.peakDetails[peak] = peakStat;
  }

  return stats;
}

function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(2)} MB`;
}

describe('Inverted index diagnostic (embedded HOSE)', () => {
  const ASPIRIN_PEAKS = [170.4, 169.8, 151.3, 139.2, 132.8, 125.7, 124.9, 123.6, 21.5];

  test('Aspirin 9 peaks - resource usage analysis', async () => {
    const stats = await diagnoseQuery({ peaks: ASPIRIN_PEAKS, tolerance: 2.0, maxResults: 10 });

    console.log('\n=== ASPIRIN REVERSE LOOKUP (EMBEDDED HOSE) ===');
    console.log(`Peaks: ${ASPIRIN_PEAKS.join(', ')}`);
    console.log(`Tolerance: ±2.0 ppm, Max results/peak: 10`);
    console.log('');
    console.log(`PPM index files loaded: ${stats.ppmFilesLoaded.size}`);
    console.log(`PPM bytes total: ${formatBytes(stats.ppmBytesTotal)}`);
    console.log(`Chunk files loaded: 0 (HOSE embedded in index)`);
    console.log(`TOTAL NETWORK: ${formatBytes(stats.ppmBytesTotal)}`);
    console.log('');
    console.log(`Total entries scanned: ${stats.totalEntries}`);
    console.log(`Entries within tolerance: ${stats.entriesWithinTolerance}`);

    for (const [peak, detail] of Object.entries(stats.peakDetails)) {
      console.log(`\n--- Peak ${peak} ppm ---`);
      console.log(`  PPM files: ${detail.ppmRange[0]}-${detail.ppmRange[1]} (${detail.ppmFilesNeeded} files)`);
      console.log(`  Entries scanned: ${detail.entriesScanned}`);
      console.log(`  Within tolerance: ${detail.entriesWithinTolerance}`);
      console.log(`  Top 3 results:`);
      for (const r of detail.topResults.slice(0, 3)) {
        console.log(`    shift=${r.shift} err=${r.error.toFixed(3)} hose=${r.hose.substring(0, 50)}...`);
      }
    }

    expect(stats.ppmFilesLoaded.size).toBeGreaterThan(0);
    expect(stats.ppmBytesTotal).toBeLessThan(30 * 1024 * 1024); // should be well under 30 MB
  }, 30000);

  test('Single peak (21.5) - minimal resource usage', async () => {
    const stats = await diagnoseQuery({ peaks: [21.5], tolerance: 2.0, maxResults: 10 });

    console.log('\n=== SINGLE PEAK (21.5) EMBEDDED HOSE ===');
    console.log(`PPM files: ${stats.ppmFilesLoaded.size}, bytes: ${formatBytes(stats.ppmBytesTotal)}`);
    console.log(`TOTAL: ${formatBytes(stats.ppmBytesTotal)}`);

    expect(stats.ppmFilesLoaded.size).toBeLessThanOrEqual(5);
  }, 30000);

  test('Morphine 17 peaks - worst case analysis', async () => {
    const MORPHINE_PEAKS = [24.4, 35.3, 40.9, 42.9, 43.1, 59.4, 66.4, 91.3, 112.9, 118.8, 119.5, 128.2, 130.1, 131.7, 142.7, 146.4, 160.4];
    const stats = await diagnoseQuery({ peaks: MORPHINE_PEAKS, tolerance: 2.0, maxResults: 10 });

    console.log('\n=== MORPHINE (17 peaks) EMBEDDED HOSE ===');
    console.log(`PPM files: ${stats.ppmFilesLoaded.size}, bytes: ${formatBytes(stats.ppmBytesTotal)}`);
    console.log(`TOTAL: ${formatBytes(stats.ppmBytesTotal)}`);

    expect(stats.ppmBytesTotal).toBeLessThan(50 * 1024 * 1024);
  }, 60000);
});
