"""
Quick search function that returns formatted string
"""

from nmr_search import load_database, search_by_shift
from hose_decoder import describe_environment


def search_nmr(shift, tolerance=0.5, max_results=10, max_db_entries=100000):
    """
    Search NMR database and return formatted string

    Args:
        shift: Chemical shift in ppm
        tolerance: Search tolerance in ppm (default: 0.5)
        max_results: Maximum number of results (default: 10)
        max_db_entries: Maximum database entries to load (default: 100000)

    Returns:
        Formatted string with search results
    """
    # Load database
    database = load_database(max_entries=max_db_entries)

    # Search
    results = search_by_shift(database, shift, tolerance, max_results)

    # Format results
    if not results:
        return f"No matches found for {shift} ppm (tolerance: +/-{tolerance} ppm)"

    output = f"Search Results for {shift} ppm (tolerance: +/-{tolerance} ppm)\n"
    output += "=" * 80 + "\n"
    output += f"Found {len(results)} matches\n"

    for idx, entry in enumerate(results, 1):
        output += f"\n[{idx}] Chemical Shift: {entry['avg_shift']:.2f} ppm (diff: +/-{entry['shift_difference']:.2f})\n"
        output += f"    Nucleus: {entry['nucleus']}\n"
        output += f"    Solvent: {entry['solvent']}\n"

        if entry.get("hose_code"):
            env_desc = describe_environment(entry["hose_code"])
            output += f"    Environment: {env_desc}\n"
            output += f"    HOSE: {entry['hose_code'][:60]}...\n"

    return output


# Example usage
if __name__ == "__main__":
    # Test with the compound you provided
    result = search_nmr(26.3, tolerance=1.0, max_results=3)
    print(result)

    print("\n" + "=" * 80 + "\n")

    # Another test
    result = search_nmr(77.0, tolerance=0.5, max_results=5)
    print(result)
