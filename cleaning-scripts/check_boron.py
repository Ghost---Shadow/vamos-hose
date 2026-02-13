"""
Check if any entries in the NMR database contain Boron (B)
"""

import sys


def check_for_boron(max_entries=None):
    """
    Search the entire database for entries containing Boron

    Returns:
        List of entries containing boron
    """
    print("Searching NMR database for Boron-containing entries...")
    print("=" * 80)

    boron_entries = []
    total_checked = 0

    with open("nmrshiftdb2/nmrshiftdb.csv", "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            total_checked += 1

            # Check if line contains 'B' in HOSE code
            # Need to be careful - 'B' could be in solvent names, etc.
            # So we look specifically in the structure/HOSE code part

            parts = line.split("_")
            if len(parts) >= 3:
                structure = parts[2]  # The structure/HOSE code part

                # Extract HOSE code
                if ";" in structure:
                    hose_code = structure.split(";")[1] if len(structure.split(";")) > 1 else ""
                else:
                    hose_code = structure

                # Check for Boron in HOSE code
                if "B" in hose_code:
                    boron_entries.append(
                        {
                            "line_number": line_num,
                            "full_entry": line,
                            "hose_code": hose_code,
                        }
                    )

                    # Print as we find them
                    print(f"\nFound Boron entry #{len(boron_entries)} (line {line_num}):")
                    print(f"  Full entry: {line[:100]}...")
                    print(f"  HOSE code: {hose_code[:80]}...")

            # Progress update
            if total_checked % 100000 == 0:
                print(f"\rChecked {total_checked:,} entries...", end="", flush=True)

            # Stop if we've checked enough
            if max_entries and total_checked >= max_entries:
                break

    print(f"\n\n" + "=" * 80)
    print(f"Search complete!")
    print(f"Total entries checked: {total_checked:,}")
    print(f"Entries containing Boron: {len(boron_entries)}")

    return boron_entries


def analyze_boron_entries(entries):
    """Analyze the boron entries in detail"""
    if not entries:
        print("\nNo boron entries found.")
        return

    print("\n" + "=" * 80)
    print("DETAILED ANALYSIS OF BORON ENTRIES")
    print("=" * 80)

    for idx, entry in enumerate(entries, 1):
        print(f"\n{'='*80}")
        print(f"Boron Entry #{idx}")
        print(f"{'='*80}")
        print(f"Line number: {entry['line_number']}")
        print(f"\nFull entry:")
        print(f"  {entry['full_entry']}")

        # Parse the entry
        parts = entry["full_entry"].split("_")
        if len(parts) >= 4:
            solvent = parts[0]
            nucleus = parts[1]
            structure = parts[2]

            # Try to extract shift data
            try:
                shift_parts = parts[3:]
                shift_str = "_".join(shift_parts)
                shift_values = shift_str.split("_")

                if len(shift_values) >= 4:
                    min_shift = float(shift_values[0])
                    max_shift = float(shift_values[1])
                    avg_shift = float(shift_values[2])
                    count = int(shift_values[3])

                    print(f"\nParsed data:")
                    print(f"  Solvent: {solvent}")
                    print(f"  Nucleus: {nucleus}")
                    print(f"  Structure: {structure}")
                    print(f"  Chemical shift: {avg_shift} ppm (range: {min_shift}-{max_shift})")
                    print(f"  Number of measurements: {count}")
            except Exception as e:
                print(f"\n  Could not parse shift data: {e}")

        print(f"\nHOSE code: {entry['hose_code']}")

        # Count boron atoms
        boron_count = entry["hose_code"].count("B")
        print(f"Boron occurrences in HOSE code: {boron_count}")


if __name__ == "__main__":
    # Check command line arguments
    if len(sys.argv) > 1:
        max_entries = int(sys.argv[1])
        print(f"Checking first {max_entries:,} entries...\n")
    else:
        max_entries = None
        print("Checking entire database...\n")

    # Search for boron
    boron_entries = check_for_boron(max_entries=max_entries)

    # Analyze findings
    analyze_boron_entries(boron_entries)

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    if boron_entries:
        print(f"✓ Found {len(boron_entries)} entries containing Boron!")
        print(f"\nBoron-containing entries are at lines:")
        for entry in boron_entries[:10]:  # Show first 10
            print(f"  - Line {entry['line_number']}")
        if len(boron_entries) > 10:
            print(f"  ... and {len(boron_entries) - 10} more")
    else:
        print("✗ No Boron-containing entries found in the database.")
