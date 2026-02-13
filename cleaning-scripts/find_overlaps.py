"""
Find entries with 100% overlapping chemical shift ranges
These are problematic because we can't distinguish between them based on shift alone
"""

import sys
from collections import defaultdict
from tqdm import tqdm
from hose_decoder import describe_environment


def ranges_overlap_100_percent(range1, range2):
    """
    Check if two ranges have 100% overlap
    This means one range is completely contained within the other (or they're identical)

    Args:
        range1: tuple (min1, max1)
        range2: tuple (min2, max2)

    Returns:
        True if one range completely contains the other
    """
    min1, max1 = range1
    min2, max2 = range2

    # Check if range1 completely contains range2
    if min1 <= min2 and max1 >= max2:
        return True

    # Check if range2 completely contains range1
    if min2 <= min1 and max2 >= max1:
        return True

    return False


def load_and_parse_database(max_entries=None):
    """Load and parse the NMR database with progress bar"""
    print("Loading NMR database...")

    entries = []
    errors = 0

    # First count total lines for progress bar
    with open("nmrshiftdb2/nmrshiftdb.csv", "r", encoding="utf-8") as f:
        total_lines = sum(1 for _ in f)

    if max_entries:
        total_lines = min(total_lines, max_entries)

    with open("nmrshiftdb2/nmrshiftdb.csv", "r", encoding="utf-8") as f:
        for idx, line in enumerate(tqdm(f, total=total_lines, desc="Parsing entries")):
            if max_entries and idx >= max_entries:
                break

            line = line.strip()
            if not line:
                continue

            try:
                parts = line.split("_")
                if len(parts) >= 4:
                    solvent = parts[0]
                    nucleus = parts[1]
                    structure = parts[2]

                    # Extract chemical shifts
                    shift_parts = parts[3:]
                    shift_str = "_".join(shift_parts)
                    shift_values = shift_str.split("_")

                    if len(shift_values) >= 4:
                        min_shift = float(shift_values[0])
                        max_shift = float(shift_values[1])
                        avg_shift = float(shift_values[2])
                        count = int(shift_values[3])

                        # Extract HOSE code
                        hose_code = ""
                        if ";" in structure:
                            hose_code = structure.split(";")[1] if len(structure.split(";")) > 1 else ""

                        entries.append({
                            "line_number": idx + 1,
                            "full_entry": line,
                            "solvent": solvent,
                            "nucleus": nucleus,
                            "structure": structure,
                            "hose_code": hose_code,
                            "min_shift": min_shift,
                            "max_shift": max_shift,
                            "avg_shift": avg_shift,
                            "count": count,
                        })
            except Exception as e:
                errors += 1
                continue

    print(f"\nSuccessfully parsed: {len(entries):,} entries")
    print(f"Parsing errors: {errors:,}")

    return entries


def find_overlapping_ranges(entries, tolerance=0.0):
    """
    Find all pairs of entries with 100% overlapping ranges

    Args:
        entries: List of database entries
        tolerance: Additional tolerance in ppm (default: 0.0 for exact overlaps)

    Returns:
        List of overlapping groups
    """
    print(f"\nSearching for overlapping ranges (tolerance: {tolerance} ppm)...")

    overlaps = []

    # Compare all pairs
    for i in tqdm(range(len(entries)), desc="Finding overlaps"):
        entry1 = entries[i]
        range1 = (entry1["min_shift"] - tolerance, entry1["max_shift"] + tolerance)

        for j in range(i + 1, len(entries)):
            entry2 = entries[j]
            range2 = (entry2["min_shift"] - tolerance, entry2["max_shift"] + tolerance)

            # Check if ranges overlap 100%
            if ranges_overlap_100_percent(range1, range2):
                overlaps.append((entry1, entry2))

    return overlaps


def group_overlaps(overlaps):
    """
    Group overlapping entries together
    Creates groups where all members have overlapping ranges with at least one other member
    """
    # Build adjacency list
    graph = defaultdict(set)
    all_entries = set()

    for entry1, entry2 in overlaps:
        id1 = entry1["line_number"]
        id2 = entry2["line_number"]
        graph[id1].add(id2)
        graph[id2].add(id1)
        all_entries.add(id1)
        all_entries.add(id2)

    # Find connected components using iterative BFS (avoid recursion depth issues)
    visited = set()
    groups = []

    for entry_id in all_entries:
        if entry_id not in visited:
            group = set()
            queue = [entry_id]

            while queue:
                node = queue.pop(0)
                if node in visited:
                    continue

                visited.add(node)
                group.add(node)

                for neighbor in graph[node]:
                    if neighbor not in visited:
                        queue.append(neighbor)

            groups.append(group)

    return groups


def analyze_overlaps(overlaps, entries_dict):
    """Analyze and display overlap statistics"""
    if not overlaps:
        print("\nNo overlapping ranges found!")
        return

    print(f"\n{'='*80}")
    print(f"OVERLAP ANALYSIS")
    print(f"{'='*80}")
    print(f"Total overlapping pairs: {len(overlaps):,}")

    # Group overlaps
    print("\nGrouping overlapping entries...")
    groups = group_overlaps(overlaps)
    print(f"Number of overlap groups: {len(groups):,}")

    # Analyze group sizes
    group_sizes = [len(g) for g in groups]
    print(f"\nGroup size statistics:")
    print(f"  Smallest group: {min(group_sizes)} entries")
    print(f"  Largest group: {max(group_sizes)} entries")
    print(f"  Average group size: {sum(group_sizes) / len(group_sizes):.1f} entries")

    # Show largest groups
    print(f"\n{'='*80}")
    print(f"TOP 10 LARGEST OVERLAP GROUPS")
    print(f"{'='*80}")

    sorted_groups = sorted(groups, key=lambda g: len(g), reverse=True)

    for idx, group in enumerate(sorted_groups[:10], 1):
        print(f"\n{'='*80}")
        print(f"Group #{idx}: {len(group)} entries with overlapping ranges")
        print(f"{'='*80}")

        group_entries = [entries_dict[line_num] for line_num in sorted(group)]

        # Find shift range of entire group
        all_mins = [e["min_shift"] for e in group_entries]
        all_maxs = [e["max_shift"] for e in group_entries]
        group_min = min(all_mins)
        group_max = max(all_maxs)

        print(f"Overall range: {group_min:.2f} - {group_max:.2f} ppm")

        # Show first 5 entries in group
        for i, entry in enumerate(group_entries[:5], 1):
            print(f"\n  Entry {i} (line {entry['line_number']}):")
            print(f"    Range: {entry['min_shift']:.2f} - {entry['max_shift']:.2f} ppm (avg: {entry['avg_shift']:.2f})")
            print(f"    Nucleus: {entry['nucleus'][:50]}")
            print(f"    Solvent: {entry['solvent']}")

            if entry["hose_code"]:
                env_desc = describe_environment(entry["hose_code"])
                print(f"    Environment: {env_desc[:70]}")
                print(f"    HOSE: {entry['hose_code'][:60]}...")

        if len(group_entries) > 5:
            print(f"\n  ... and {len(group_entries) - 5} more entries in this group")


def export_overlaps(overlaps, entries_dict, output_file="overlapping_ranges.txt"):
    """Export overlapping entries to a file"""
    print(f"\nExporting overlaps to {output_file}...")

    groups = group_overlaps(overlaps)
    sorted_groups = sorted(groups, key=lambda g: len(g), reverse=True)

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("NMR Database Overlapping Ranges Report\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total overlapping pairs: {len(overlaps):,}\n")
        f.write(f"Number of overlap groups: {len(groups):,}\n\n")

        for idx, group in enumerate(sorted_groups, 1):
            f.write("\n" + "="*80 + "\n")
            f.write(f"Group #{idx}: {len(group)} entries\n")
            f.write("="*80 + "\n")

            group_entries = [entries_dict[line_num] for line_num in sorted(group)]

            for entry in group_entries:
                f.write(f"\nLine {entry['line_number']}:\n")
                f.write(f"  {entry['full_entry']}\n")
                f.write(f"  Range: {entry['min_shift']:.2f} - {entry['max_shift']:.2f} ppm\n")

    print(f"Export complete: {output_file}")


if __name__ == "__main__":
    # Parse arguments
    max_entries = int(sys.argv[1]) if len(sys.argv) > 1 else None
    tolerance = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0

    if max_entries:
        print(f"Analyzing first {max_entries:,} entries...")
    else:
        print("Analyzing entire database...")

    print(f"Tolerance: {tolerance} ppm\n")

    # Load database
    entries = load_and_parse_database(max_entries=max_entries)

    # Create lookup dict
    entries_dict = {e["line_number"]: e for e in entries}

    # Find overlaps
    overlaps = find_overlapping_ranges(entries, tolerance=tolerance)

    # Analyze
    analyze_overlaps(overlaps, entries_dict)

    # Export
    if overlaps:
        export_overlaps(overlaps, entries_dict)

    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
