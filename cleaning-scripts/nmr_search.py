"""
NMR Database Search Tool
Search for peaks by chemical shift values
"""

import sys
from hose_decoder import describe_environment


def load_database(max_entries=None):
    """Load the NMR database"""
    print("Loading NMR database...")
    data = []

    with open("nmrshiftdb2/nmrshiftdb.csv", "r", encoding="utf-8") as f:
        for idx, line in enumerate(f):
            if max_entries and idx >= max_entries:
                break

            line = line.strip()
            if not line:
                continue

            try:
                parts = line.split("_")
                if len(parts) >= 4:
                    solvent = parts[0]
                    nucleus_info = parts[1]
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
                            hose_code = (
                                structure.split(";")[1]
                                if len(structure.split(";")) > 1
                                else ""
                            )

                        data.append(
                            {
                                "solvent": solvent,
                                "nucleus": nucleus_info,
                                "structure": structure,
                                "hose_code": hose_code,
                                "min_shift": min_shift,
                                "max_shift": max_shift,
                                "avg_shift": avg_shift,
                                "count": count,
                            }
                        )
            except Exception as e:
                continue

    print(f"Loaded {len(data)} entries")
    return data


def search_by_shift(database, target_shift, tolerance=0.5, max_results=20):
    """
    Search for peaks matching a chemical shift value

    Args:
        database: List of NMR data entries
        target_shift: The chemical shift to search for (ppm)
        tolerance: Search tolerance in ppm (default: 0.5)
        max_results: Maximum number of results to return

    Returns:
        List of matching entries
    """
    matches = []

    for entry in database:
        # Check if target shift falls within the range or is close to average
        if abs(entry["avg_shift"] - target_shift) <= tolerance:
            diff = abs(entry["avg_shift"] - target_shift)
            entry["shift_difference"] = diff
            matches.append(entry)

    # Sort by difference (closest matches first)
    matches.sort(key=lambda x: x["shift_difference"])

    return matches[:max_results]


def search_by_peak_list(database, peak_list, tolerance=0.5, min_matches=3):
    """
    Search for compounds matching a list of peaks

    Args:
        database: List of NMR data entries
        peak_list: List of chemical shift values
        tolerance: Search tolerance in ppm
        min_matches: Minimum number of peaks that must match

    Returns:
        Dictionary of compounds with their matched peaks
    """
    from collections import defaultdict

    # Group database by structure
    compounds = defaultdict(list)
    for entry in database:
        compounds[entry["structure"]].append(entry)

    matches = []

    for structure, peaks in compounds.items():
        matched_peaks = []
        db_shifts = [p["avg_shift"] for p in peaks]

        # For each query peak, find if there's a matching DB peak
        for query_shift in peak_list:
            for idx, db_shift in enumerate(db_shifts):
                if abs(query_shift - db_shift) <= tolerance:
                    matched_peaks.append(
                        {
                            "query_shift": query_shift,
                            "db_shift": db_shift,
                            "db_entry": peaks[idx],
                            "difference": abs(query_shift - db_shift),
                        }
                    )
                    break

        # If enough peaks matched, add this compound to results
        if len(matched_peaks) >= min_matches:
            matches.append(
                {
                    "structure": structure,
                    "all_peaks": peaks,
                    "matched_peaks": matched_peaks,
                    "match_count": len(matched_peaks),
                    "total_peaks": len(peaks),
                    "match_ratio": len(matched_peaks) / len(peak_list),
                }
            )

    # Sort by number of matches and match ratio
    matches.sort(key=lambda x: (x["match_count"], x["match_ratio"]), reverse=True)

    return matches


def display_search_results(results, show_hose=True):
    """Display search results in a readable format"""
    if not results:
        print("No matches found!")
        return

    print(f"\nFound {len(results)} matches:\n")
    print("=" * 90)

    for idx, entry in enumerate(results, 1):
        print(f"\n{idx}. Chemical Shift: {entry['avg_shift']:.2f} ppm")
        print(f"   Difference: +/-{entry['shift_difference']:.2f} ppm")
        print(f"   Nucleus: {entry['nucleus']}")
        print(f"   Solvent: {entry['solvent']}")

        if show_hose and entry["hose_code"]:
            env_desc = describe_environment(entry["hose_code"])
            print(f"   Environment: {env_desc}")
            print(f"   HOSE code: {entry['hose_code'][:60]}...")

        print("-" * 90)


def display_compound_matches(matches, top_n=10):
    """Display compound matching results"""
    if not matches:
        print("No matching compounds found!")
        return

    print(f"\nFound {len(matches)} matching compounds:\n")
    print("=" * 90)

    for idx, match in enumerate(matches[:top_n], 1):
        print(f"\n{idx}. COMPOUND MATCH")
        print(f"   Total peaks in compound: {match['total_peaks']}")
        print(f"   Matched peaks: {match['match_count']}")
        print(f"   Match ratio: {match['match_ratio']:.1%}")

        # Get environment descriptions
        if match["all_peaks"][0]["hose_code"]:
            env_desc = describe_environment(match["all_peaks"][0]["hose_code"])
            print(f"   Environment: {env_desc}")

        print(f"\n   Matched Peaks:")
        for peak in match["matched_peaks"]:
            print(
                f"      Query: {peak['query_shift']:.1f} ppm -> "
                f"Database: {peak['db_shift']:.1f} ppm "
                f"(diff = {peak['difference']:.2f})"
            )

        # Show all compound peaks
        all_shifts = [p["avg_shift"] for p in match["all_peaks"]]
        print(
            f"\n   All peaks in this compound: {', '.join([f'{s:.1f}' for s in sorted(all_shifts, reverse=True)])}"
        )
        print("-" * 90)


# Main interface
if __name__ == "__main__":
    print("NMR Database Search Tool")
    print("=" * 90)

    # Load database (use subset for speed, or None for full database)
    database = load_database(max_entries=100000)  # Load first 100k entries

    print("\nSearch Options:")
    print("1. Search by single chemical shift")
    print("2. Search by list of peaks (find matching compounds)")
    print()

    choice = input("Enter choice (1 or 2): ").strip()

    if choice == "1":
        # Single peak search
        shift = float(input("\nEnter chemical shift (ppm): "))
        tolerance = float(input("Enter tolerance (ppm) [default: 0.5]: ") or "0.5")

        print(f"\nSearching for peaks near {shift} ppm (±{tolerance} ppm)...")
        results = search_by_shift(database, shift, tolerance, max_results=20)
        display_search_results(results)

    elif choice == "2":
        # Multiple peak search
        print("\nEnter your peak list (comma-separated, e.g., 128.5, 77.2, 42.1):")
        peak_input = input("Peaks: ")
        peaks = [float(x.strip()) for x in peak_input.split(",")]

        tolerance = float(input("Enter tolerance (ppm) [default: 1.0]: ") or "1.0")
        min_matches = int(
            input(f"Minimum matches required [default: {min(3, len(peaks))}]: ")
            or str(min(3, len(peaks)))
        )

        print(f"\nSearching for compounds matching {len(peaks)} peaks...")
        matches = search_by_peak_list(
            database, peaks, tolerance=tolerance, min_matches=min_matches
        )
        display_compound_matches(matches)

    else:
        print("Invalid choice!")
