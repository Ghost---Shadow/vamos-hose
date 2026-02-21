/**
 * Shard the monolithic hose_shift_lookup.json into 256 JS chunk files.
 *
 * Each chunk exports a default object mapping HOSE codes to shift data.
 * Chunks are named chunk_000.js through chunk_255.js.
 *
 * The hash function must match the one in src/database.js so lookups
 * route to the correct chunk.
 *
 * Usage: node cleaning-scripts/shard_database.mjs
 */
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const NUM_CHUNKS = 256;
const INPUT = path.join(__dirname, 'hose_shift_lookup.json');
const OUTPUT_DIR = path.join(__dirname, '..', 'chunks');

function hashCode(str) {
  let h = 0;
  for (let i = 0; i < str.length; i++) {
    h = ((h << 5) - h + str.charCodeAt(i)) | 0;
  }
  return Math.abs(h);
}

console.log('Loading database...');
const db = JSON.parse(fs.readFileSync(INPUT, 'utf8'));
const keys = Object.keys(db);
console.log(`Loaded ${keys.length} keys`);

// Create output directory
if (!fs.existsSync(OUTPUT_DIR)) {
  fs.mkdirSync(OUTPUT_DIR, { recursive: true });
}

// Distribute keys into buckets
const buckets = new Array(NUM_CHUNKS).fill(null).map(() => ({}));
for (const k of keys) {
  const bucket = hashCode(k) % NUM_CHUNKS;
  buckets[bucket][k] = db[k];
}

// Write each bucket as a JS file
let totalSize = 0;
for (let i = 0; i < NUM_CHUNKS; i++) {
  const chunkName = `chunk_${String(i).padStart(3, '0')}.js`;
  const chunkPath = path.join(OUTPUT_DIR, chunkName);
  const content = `export default ${JSON.stringify(buckets[i])};\n`;
  fs.writeFileSync(chunkPath, content);
  const size = Buffer.byteLength(content);
  totalSize += size;
  if (i % 32 === 0) {
    const keyCount = Object.keys(buckets[i]).length;
    console.log(`  ${chunkName}: ${keyCount} keys, ${(size / 1024 / 1024).toFixed(2)} MB`);
  }
}

console.log(`\nWrote ${NUM_CHUNKS} chunks to ${OUTPUT_DIR}`);
console.log(`Total size: ${(totalSize / 1024 / 1024).toFixed(1)} MB`);
console.log(`Avg chunk: ${(totalSize / NUM_CHUNKS / 1024 / 1024).toFixed(2)} MB`);
