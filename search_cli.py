"""
Command-line NMR Database Search Tool
Usage: python search_cli.py <shift> [tolerance] [max_results]
"""

import sys
from nmr_search import load_database, search_by_shift
from hose_decoder import describe_environment


def format_result(entry, index=1):
    """Format a single search result as a string"""
    result_str = f"\n[{index}] Chemical Shift: {entry['avg_shift']:.2f} ppm (diff: +/-{entry['shift_difference']:.2f})\n"
    result_str += f"    Nucleus: {entry['nucleus']}\n"
    result_str += f"    Solvent: {entry['solvent']}\n"

    if entry.get("hose_code"):
        env_desc = describe_environment(entry["hose_code"])
        result_str += f"    Environment: {env_desc}\n"
        result_str += f"    HOSE: {entry['hose_code'][:60]}...\n"

    return result_str


def format_all_results(results, query_shift, tolerance):
    """Format all search results as a single string"""
    if not results:
        return f"No matches found for {query_shift} ppm (tolerance: +/-{tolerance} ppm)"

    output = f"Search Results for {query_shift} ppm (tolerance: +/-{tolerance} ppm)\n"
    output += "=" * 80 + "\n"
    output += f"Found {len(results)} matches\n"

    for idx, entry in enumerate(results, 1):
        output += format_result(entry, idx)

    return output


def main():
    if len(sys.argv) < 2:
        print("Usage: python search_cli.py <shift> [tolerance] [max_results] [max_db_entries]")
        print("\nExamples:")
        print("  python search_cli.py 77.0")
        print("  python search_cli.py 77.0 1.0")
        print("  python search_cli.py 77.0 1.0 10")
        print("  python search_cli.py 77.0 1.0 10 50000")
        sys.exit(1)

    # Parse arguments
    target_shift = float(sys.argv[1])
    tolerance = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
    max_results = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    max_db_entries = int(sys.argv[4]) if len(sys.argv) > 4 else 100000

    # Load database
    print(f"Loading database (up to {max_db_entries} entries)...", file=sys.stderr)
    database = load_database(max_entries=max_db_entries)

    # Search
    print(f"Searching for {target_shift} ppm...", file=sys.stderr)
    results = search_by_shift(database, target_shift, tolerance, max_results)

    # Format and print results
    output = format_all_results(results, target_shift, tolerance)
    print(output)


if __name__ == "__main__":
    main()
