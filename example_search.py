"""
Example NMR Database Search
Demonstrates how to search for peaks
"""

from nmr_search import load_database, search_by_shift, search_by_peak_list
from nmr_search import display_search_results, display_compound_matches

print("Example NMR Database Search")
print("=" * 90)

# Load database (subset for speed)
database = load_database(max_entries=100000)

# ============================================================================
# EXAMPLE 1: Search for a single peak
# ============================================================================
print("\n" + "=" * 90)
print("EXAMPLE 1: Searching for a single peak at ~77 ppm")
print("(This is typical for CDCl3 or alcohol/ether carbons)")
print("=" * 90)

results = search_by_shift(database, target_shift=77.0, tolerance=1.0, max_results=10)
display_search_results(results, show_hose=True)

# ============================================================================
# EXAMPLE 2: Search for aromatic peaks
# ============================================================================
print("\n\n" + "=" * 90)
print("EXAMPLE 2: Searching for aromatic peak at ~128 ppm")
print("(Typical aromatic carbon region)")
print("=" * 90)

results = search_by_shift(database, target_shift=128.0, tolerance=2.0, max_results=10)
display_search_results(results, show_hose=True)

# ============================================================================
# EXAMPLE 3: Search for a compound by multiple peaks
# ============================================================================
print("\n\n" + "=" * 90)
print("EXAMPLE 3: Searching for compounds matching a peak list")
print("Example peaks: 128.5, 77.2, 42.1, 25.3 ppm")
print("(Could be a compound with aromatic ring + aliphatic chain)")
print("=" * 90)

query_peaks = [128.5, 77.2, 42.1, 25.3]
print(f"\nQuery peaks: {', '.join([f'{p:.1f}' for p in query_peaks])} ppm")

matches = search_by_peak_list(
    database, query_peaks, tolerance=2.0, min_matches=2  # At least 2 peaks must match
)

display_compound_matches(matches, top_n=5)

# ============================================================================
# EXAMPLE 4: Search for aliphatic carbons
# ============================================================================
print("\n\n" + "=" * 90)
print("EXAMPLE 4: Searching for aliphatic methyl peak at ~20 ppm")
print("(Typical for CH3 groups)")
print("=" * 90)

results = search_by_shift(database, target_shift=20.0, tolerance=2.0, max_results=10)
display_search_results(results, show_hose=True)

print("\n" + "=" * 90)
print("Search complete!")
print("=" * 90)
print(
    "\nTip: Run 'python nmr_search.py' for interactive search with your own peaks!"
)
