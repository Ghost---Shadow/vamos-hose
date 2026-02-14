"""
NMR Database Preprocessor

Reads the NMR shift database CSV, converts each HOSE code to a SMILES fragment,
and outputs a JSON lookup table mapping HOSE codes to shift data + SMILES.

Primary key: HOSE code (specific, no information loss)
Secondary field: SMILES fragment (for display/matching in the JS runtime)

Usage:
    python preprocess_database.py [max_entries]
    python preprocess_database.py           # process entire database
    python preprocess_database.py 10000     # process first 10k entries
"""

import json
import sys
import time
from hose_to_smiles import hose_to_smiles, extract_central_atom


def parse_line(line):
    """Parse a database line into its components.

    Format: Solvent_Nucleus;HOSE_min_max_avg_count

    The tricky part: solvent names can contain underscores.
    Strategy:
      1. rsplit on '_' with maxsplit=4 to peel off the 4 numeric fields
      2. Split the front on ';' to separate Solvent_Nucleus from HOSE code
      3. rsplit the before-semicolon part on '_' with maxsplit=1 for nucleus
    """
    line = line.strip()
    if not line:
        return None

    # Step 1: peel off the 4 numeric fields from the right
    parts = line.rsplit("_", 4)
    if len(parts) < 5:
        return None

    front = parts[0]
    try:
        min_shift = float(parts[1])
        max_shift = float(parts[2])
        avg_shift = float(parts[3])
        count = int(parts[4])
    except (ValueError, IndexError):
        return None

    # Step 2: split on ';' to get HOSE code
    if ";" not in front:
        return None

    semi_pos = front.index(";")
    before_semi = front[:semi_pos]
    hose_code = front[semi_pos + 1:]

    # Step 3: last underscore-delimited field before ';' is the nucleus
    nuc_parts = before_semi.rsplit("_", 1)
    if len(nuc_parts) == 2:
        solvent = nuc_parts[0]
        nucleus = nuc_parts[1]
    else:
        solvent = ""
        nucleus = nuc_parts[0]

    return {
        "solvent": solvent,
        "nucleus": nucleus,
        "hose_code": hose_code,
        "min_shift": min_shift,
        "max_shift": max_shift,
        "avg_shift": avg_shift,
        "count": count,
    }


def get_central_atom(nucleus_str):
    """Extract central atom symbol from nucleus string.
    Examples: 'C-4' -> 'C', 'C-3-6' -> 'C', 'H-1' -> 'H', 'N-15' -> 'N'
    """
    if not nucleus_str:
        return "C"
    atom = nucleus_str.split("-")[0].strip()
    return atom if atom else "C"


def preprocess(max_entries=None, db_path="../nmrshiftdb2/nmrshiftdb.csv"):
    """Process the database and build a HOSE->shift lookup table with SMILES."""
    print("NMR Database Preprocessor")
    print("=" * 70)

    # Stats
    total = 0
    parsed = 0
    converted = 0
    failed_parse = 0
    failed_convert = 0
    skipped_no_hose = 0

    # Lookup table: hose_code -> {"n": nucleus, "s": smiles, solvent: {min, max, sum, cnt}}
    lookup = {}

    start_time = time.time()

    with open(db_path, "r", encoding="utf-8") as f:
        for idx, line in enumerate(f):
            if max_entries and idx >= max_entries:
                break

            total += 1
            entry = parse_line(line)

            if entry is None:
                failed_parse += 1
                continue

            parsed += 1
            hose = entry["hose_code"]

            if not hose:
                skipped_no_hose += 1
                continue

            central = get_central_atom(entry["nucleus"])

            # Use HOSE code as the primary key
            key = hose
            if key not in lookup:
                smiles = hose_to_smiles(hose, central)
                if smiles is None:
                    failed_convert += 1
                    continue
                lookup[key] = {"n": central, "s": smiles}

            converted += 1
            solvent = entry["solvent"]
            if solvent not in lookup[key]:
                lookup[key][solvent] = {
                    "min": entry["min_shift"],
                    "max": entry["max_shift"],
                    "sum": entry["avg_shift"] * entry["count"],
                    "cnt": entry["count"],
                }
            else:
                rec = lookup[key][solvent]
                rec["min"] = min(rec["min"], entry["min_shift"])
                rec["max"] = max(rec["max"], entry["max_shift"])
                rec["sum"] += entry["avg_shift"] * entry["count"]
                rec["cnt"] += entry["count"]

            # Progress reporting
            if total % 100000 == 0:
                elapsed = time.time() - start_time
                rate = total / elapsed if elapsed > 0 else 0
                print(
                    f"  Processed {total:,} lines "
                    f"({converted:,} converted, {failed_parse + failed_convert:,} failed) "
                    f"[{rate:.0f} lines/sec]"
                )

    elapsed = time.time() - start_time

    # Finalize: convert sum -> avg, round values
    for hose_key, entry in lookup.items():
        for k, v in entry.items():
            if k in ("n", "s"):
                continue
            v["avg"] = round(v["sum"] / v["cnt"], 4) if v["cnt"] > 0 else 0
            v["min"] = round(v["min"], 4)
            v["max"] = round(v["max"], 4)
            del v["sum"]

    # Write output
    output_path = "hose_shift_lookup.json"
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(lookup, f, separators=(",", ":"))

    # Pretty-printed sample for inspection
    sample_path = "hose_shift_sample.json"
    sample_keys = list(lookup.keys())[:50]
    sample = {k: lookup[k] for k in sample_keys}
    with open(sample_path, "w", encoding="utf-8") as f:
        json.dump(sample, f, indent=2)

    # Stats
    total_solvent_entries = sum(
        sum(1 for k in v if k not in ("n", "s"))
        for v in lookup.values()
    )

    # Report
    print()
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"  Total lines:         {total:,}")
    print(f"  Parsed successfully: {parsed:,}")
    print(f"  Converted to SMILES: {converted:,}")
    print(f"  Failed (parse):      {failed_parse:,}")
    print(f"  Failed (convert):    {failed_convert:,}")
    print(f"  No HOSE code:        {skipped_no_hose:,}")
    print(f"  Unique HOSE keys:    {len(lookup):,}")
    print(f"  Solvent groups:      {total_solvent_entries:,}")
    print(f"  Time elapsed:        {elapsed:.1f} sec")
    print(f"  Output:              {output_path}")
    print(f"  Sample:              {sample_path}")

    conversion_rate = converted / parsed * 100 if parsed > 0 else 0
    print(f"  Conversion rate:     {conversion_rate:.1f}%")

    return lookup


if __name__ == "__main__":
    max_entries = int(sys.argv[1]) if len(sys.argv) > 1 else None
    preprocess(max_entries)
