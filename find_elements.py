"""
Find all chemical elements present in the NMR database
Scans HOSE codes to identify all atoms
"""

import sys
from collections import Counter
from tqdm import tqdm


def extract_elements_from_hose(hose_code):
    """
    Extract all elements from a HOSE code

    Returns:
        set of element symbols found
    """
    elements = set()

    # List of possible elements (ordered by length for proper matching)
    # Two-letter elements first to avoid confusion (e.g., Cl before C)
    element_list = [
        'Cl', 'Br', 'Si', 'Se', 'Sn', 'Pb', 'Bi',  # Two-letter elements
        'C', 'H', 'O', 'N', 'S', 'P', 'F', 'B', 'I',  # Single-letter common
        'K', 'V', 'W', 'X', 'Y', 'Z', 'Q',  # Other single letters
        'Al', 'As', 'Au', 'Ag', 'Ca', 'Cd', 'Co', 'Cr', 'Cu',  # More two-letter
        'Fe', 'Ga', 'Ge', 'Hg', 'Li', 'Mg', 'Mn', 'Mo', 'Na',
        'Ni', 'Pt', 'Rb', 'Re', 'Ru', 'Sb', 'Sc', 'Sr', 'Te',
        'Ti', 'Tl', 'Zn', 'Zr'
    ]

    i = 0
    while i < len(hose_code):
        # Try two-letter elements first
        if i + 1 < len(hose_code):
            two_char = hose_code[i:i+2]
            if two_char in element_list:
                elements.add(two_char)
                i += 2
                continue

        # Try single-letter elements
        one_char = hose_code[i]
        if one_char in element_list:
            elements.add(one_char)
            i += 1
        else:
            i += 1

    return elements


def find_all_elements(max_entries=None):
    """
    Scan the entire database and find all elements

    Returns:
        Counter object with element frequencies
    """
    print("Scanning NMR database for chemical elements...")
    print("=" * 80)

    element_counter = Counter()
    total_entries = 0
    entries_with_elements = 0

    # Count total lines for progress bar
    with open("nmrshiftdb2/nmrshiftdb.csv", "r", encoding="utf-8") as f:
        total_lines = sum(1 for _ in f)

    if max_entries:
        total_lines = min(total_lines, max_entries)

    # Scan database
    with open("nmrshiftdb2/nmrshiftdb.csv", "r", encoding="utf-8") as f:
        for line_num, line in enumerate(tqdm(f, total=total_lines, desc="Scanning entries"), 1):
            if max_entries and line_num > max_entries:
                break

            line = line.strip()
            if not line:
                continue

            total_entries += 1

            # Extract HOSE code from the line
            parts = line.split("_")
            if len(parts) >= 3:
                structure = parts[2]

                # Extract HOSE code (after semicolon)
                if ";" in structure:
                    hose_code = structure.split(";")[1] if len(structure.split(";")) > 1 else ""
                else:
                    hose_code = structure

                # Find elements in this HOSE code
                elements = extract_elements_from_hose(hose_code)

                if elements:
                    entries_with_elements += 1
                    for elem in elements:
                        element_counter[elem] += 1

    print(f"\n{'='*80}")
    print("Scan complete!")
    print(f"Total entries scanned: {total_entries:,}")
    print(f"Entries with identified elements: {entries_with_elements:,}")

    return element_counter


def categorize_elements(element_counter):
    """Categorize elements by type"""

    categories = {
        'Organic (CHNOPS)': ['C', 'H', 'N', 'O', 'P', 'S'],
        'Halogens': ['F', 'Cl', 'Br', 'I'],
        'Alkali/Alkaline Earth': ['Li', 'Na', 'K', 'Rb', 'Mg', 'Ca', 'Sr', 'Ba'],
        'Transition Metals': ['Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ag', 'Au', 'Pt', 'Pd', 'Ru', 'Rh', 'Ir', 'Mn', 'Cr', 'Mo', 'W', 'V', 'Ti', 'Zr', 'Hf', 'Sc', 'Y'],
        'Metalloids': ['B', 'Si', 'Ge', 'As', 'Sb', 'Te'],
        'Post-transition Metals': ['Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi'],
        'Special/Other': ['Se', 'Q', 'X', 'Y', 'Z']  # Q, X, Y, Z might be placeholders
    }

    categorized = {}
    for category, elements in categories.items():
        found = {elem: element_counter[elem] for elem in elements if elem in element_counter}
        if found:
            categorized[category] = found

    # Find uncategorized elements
    all_categorized = set()
    for elems in categories.values():
        all_categorized.update(elems)

    uncategorized = {elem: count for elem, count in element_counter.items()
                     if elem not in all_categorized}

    if uncategorized:
        categorized['Uncategorized'] = uncategorized

    return categorized


def display_results(element_counter):
    """Display element statistics"""

    if not element_counter:
        print("\nNo elements found!")
        return

    print(f"\n{'='*80}")
    print("ELEMENTS FOUND IN DATABASE")
    print(f"{'='*80}")
    print(f"\nTotal unique elements: {len(element_counter)}")
    print()

    # Sort by frequency
    sorted_elements = element_counter.most_common()

    print("All elements (sorted by frequency):")
    print(f"{'='*80}")
    print(f"{'Element':<10} {'Count':>15} {'Percentage':>12}")
    print(f"{'-'*80}")

    total_count = sum(element_counter.values())
    for elem, count in sorted_elements:
        percentage = (count / total_count) * 100
        print(f"{elem:<10} {count:>15,} {percentage:>11.2f}%")

    # Categorized view
    print(f"\n{'='*80}")
    print("ELEMENTS BY CATEGORY")
    print(f"{'='*80}")

    categorized = categorize_elements(element_counter)

    for category, elements in categorized.items():
        print(f"\n{category}:")
        print(f"{'-'*80}")
        for elem, count in sorted(elements.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_count) * 100
            print(f"  {elem:<8} {count:>15,} ({percentage:>6.2f}%)")


def export_results(element_counter, output_file="elements_found.txt"):
    """Export results to file"""
    print(f"\nExporting results to {output_file}...")

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("NMR Database - Chemical Elements Found\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total unique elements: {len(element_counter)}\n\n")

        # All elements
        f.write("All Elements (sorted by frequency):\n")
        f.write("="*80 + "\n")
        f.write(f"{'Element':<10} {'Count':>15} {'Percentage':>12}\n")
        f.write("-"*80 + "\n")

        total_count = sum(element_counter.values())
        for elem, count in element_counter.most_common():
            percentage = (count / total_count) * 100
            f.write(f"{elem:<10} {count:>15,} {percentage:>11.2f}%\n")

        # Categorized
        f.write("\n" + "="*80 + "\n")
        f.write("Elements by Category\n")
        f.write("="*80 + "\n")

        categorized = categorize_elements(element_counter)
        for category, elements in categorized.items():
            f.write(f"\n{category}:\n")
            f.write("-"*80 + "\n")
            for elem, count in sorted(elements.items(), key=lambda x: x[1], reverse=True):
                percentage = (count / total_count) * 100
                f.write(f"  {elem:<8} {count:>15,} ({percentage:>6.2f}%)\n")

    print(f"Export complete: {output_file}")


if __name__ == "__main__":
    # Parse arguments
    if len(sys.argv) > 1:
        max_entries = int(sys.argv[1])
        print(f"Scanning first {max_entries:,} entries...\n")
    else:
        max_entries = None
        print("Scanning entire database...\n")

    # Find elements
    element_counter = find_all_elements(max_entries=max_entries)

    # Display results
    display_results(element_counter)

    # Export
    export_results(element_counter)

    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
