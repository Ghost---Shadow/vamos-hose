/**
 * Build a shift-keyed inverted index from the monolithic hose_shift_lookup.json.
 *
 * Creates 250 PPM files (ppm_000.js through ppm_249.js) in inverted-index/.
 * Each file maps a 1-ppm range to an array of [chunkIndex, exactShift] pairs.
 *
 * This lets the reverse lookup (estimateFromSpectra) load only the PPM files
 * relevant to the user's peaks, instead of all 256 HOSE-code chunks (~211 MB).
 *
 * Usage: node cleaning-scripts/build_inverted_index.mjs
 */
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const NUM_CHUNKS = 256;
const NUM_PPM_BUCKETS = 250;
const INPUT = path.join(__dirname, 'hose_shift_lookup.json');
const OUTPUT_DIR = path.join(__dirname, '..', 'inverted-index');

/** Must match src/database.js hashCode */
function hashCode(str) {
  let h = 0;
  for (let i = 0; i < str.length; i++) {
    h = ((h << 5) - h + str.charCodeAt(i)) | 0;
  }
  return Math.abs(h);
}

/** Must match src/database.js computeWeightedAvg */
function computeWeightedAvg(entry) {
  let totalWeight = 0;
  let weightedSum = 0;
  for (const [key, val] of Object.entries(entry)) {
    if (key === 'n' || key === 's') continue;
    weightedSum += val.avg * val.cnt;
    totalWeight += val.cnt;
  }
  if (totalWeight === 0) return 0;
  return Math.round((weightedSum / totalWeight) * 10) / 10;
}

console.log('Loading database...');
const db = JSON.parse(fs.readFileSync(INPUT, 'utf8'));
const keys = Object.keys(db);
console.log(`Loaded ${keys.length} entries`);

if (!fs.existsSync(OUTPUT_DIR)) {
  fs.mkdirSync(OUTPUT_DIR, { recursive: true });
}

// Build PPM buckets: bucket[ppm] = [[chunkIdx, shift], ...]
const buckets = new Array(NUM_PPM_BUCKETS).fill(null).map(() => []);
let skipped = 0;

for (const hoseCode of keys) {
  const entry = db[hoseCode];
  if (entry.n !== 'C') { skipped++; continue; }

  const shift = computeWeightedAvg(entry);
  const ppm = Math.floor(shift);

  if (ppm < 0 || ppm >= NUM_PPM_BUCKETS) { skipped++; continue; }

  const chunkIdx = hashCode(hoseCode) % NUM_CHUNKS;
  buckets[ppm].push([chunkIdx, shift]);
}

// Write each PPM bucket as a JS file
let totalSize = 0;
let totalEntries = 0;
let maxBucketSize = 0;
let maxBucketPpm = 0;

for (let ppm = 0; ppm < NUM_PPM_BUCKETS; ppm++) {
  const fileName = `ppm_${String(ppm).padStart(3, '0')}.js`;
  const filePath = path.join(OUTPUT_DIR, fileName);
  const content = `export default ${JSON.stringify(buckets[ppm])};\n`;
  fs.writeFileSync(filePath, content);

  const size = Buffer.byteLength(content);
  totalSize += size;
  totalEntries += buckets[ppm].length;

  if (buckets[ppm].length > maxBucketSize) {
    maxBucketSize = buckets[ppm].length;
    maxBucketPpm = ppm;
  }
}

console.log(`\nWrote ${NUM_PPM_BUCKETS} PPM files to ${OUTPUT_DIR}`);
console.log(`Total entries indexed: ${totalEntries}`);
console.log(`Skipped (non-13C or out of range): ${skipped}`);
console.log(`Total size: ${(totalSize / 1024 / 1024).toFixed(1)} MB`);
console.log(`Avg file: ${(totalSize / NUM_PPM_BUCKETS / 1024).toFixed(1)} KB`);
console.log(`Largest bucket: ppm_${String(maxBucketPpm).padStart(3, '0')} with ${maxBucketSize} entries`);

// Distribution summary
console.log('\nDistribution (entries per 10-ppm range):');
for (let start = 0; start < NUM_PPM_BUCKETS; start += 10) {
  const end = Math.min(start + 10, NUM_PPM_BUCKETS);
  let count = 0;
  for (let p = start; p < end; p++) count += buckets[p].length;
  const bar = '#'.repeat(Math.ceil(count / 5000));
  console.log(`  ${String(start).padStart(3)}-${String(end - 1).padStart(3)} ppm: ${String(count).padStart(8)} ${bar}`);
}
