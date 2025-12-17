# NMR Database Exploration Project

## Project Overview

This project provides tools for exploring and searching a large NMR (Nuclear Magnetic Resonance) spectroscopy database containing 1,721,704 chemical shift entries. The database uses HOSE (Hierarchically Ordered Spherical Environment) codes to describe the local chemical environment around atoms.

## Database Information

### Statistics
- **Total Entries**: 1,721,704
- **Unique Elements**: 30
- **Boron-Containing Entries**: 5,740 (0.33%)
- **Format**: CSV with structure `Solvent_Nucleus_Structure;HOSE-code_min_shift_max_shift_avg_shift_count`

### Database Format

Each entry contains:
- **Solvent**: The solvent used for measurement (e.g., `C_Pyridin-D5 (C5D5N) + Methanol-D4 (CD3OD)`)
- **Nucleus**: The atom being observed (e.g., `C-4` for carbon-4)
- **Structure**: SMILES-like structure with HOSE code after semicolon
- **HOSE Code**: Hierarchical encoding of chemical environment (spheres separated by `/`)
- **Shift Data**: Four numbers representing min, max, avg shifts (ppm) and count

Example entry:
```
C_Pyridin-D5 (C5D5N) + Methanol-D4 (CD3OD)_C-4;C(=CC/HC,HHH/HHC)_17.6_21.975_26.3_4
```

This means:
- Average shift: 21.975 ppm
- Range: 17.6 - 26.3 ppm
- Based on 4 measurements

## Scripts Created

### 1. `hose_decoder.py`
**Purpose**: Decode HOSE codes into human-readable descriptions

**Key Functions**:
- `decode_hose_code(hose_code)` - Extracts atoms, bonds, and structural features
- `describe_environment(hose_code)` - Returns plain English description

**Example**:
```python
from hose_decoder import describe_environment

hose = "C(=CC/HC,HHH/HHC)"
print(describe_environment(hose))
# Output: "aromatic, bonded to C=C (5C, 6H nearby)"
```

### 2. `plot_nmr_spectra.py` & `plot_nmr_with_labels.py`
**Purpose**: Visualize NMR spectra with compound labels

**Features**:
- Loads database entries
- Groups peaks by structure
- Generates Lorentzian peak shapes for realistic spectra
- Labels peaks with decoded HOSE descriptions

**Usage**:
```bash
python plot_nmr_with_labels.py
```

### 3. `nmr_search.py`
**Purpose**: Core search functionality for the database

**Key Functions**:

#### Single Peak Search
```python
def search_by_shift(database, target_shift, tolerance=0.5, max_results=20)
```
Find entries matching a single chemical shift value.

#### Multi-Peak Search
```python
def search_by_peak_list(database, peak_list, tolerance=0.5, min_matches=3)
```
Search for compounds matching multiple peaks simultaneously.

**Example**:
```python
from nmr_search import load_database, search_by_shift

db = load_database(max_entries=100000)
results = search_by_shift(db, 77.0, tolerance=0.5, max_results=10)
```

### 4. `search_cli.py`
**Purpose**: Command-line interface for NMR searching

**Usage**:
```bash
# Basic search
python search_cli.py 77.0

# With tolerance
python search_cli.py 77.0 1.0

# With max results
python search_cli.py 77.0 1.0 10

# With database size limit
python search_cli.py 77.0 1.0 10 50000
```

**Example Output**:
```
Search Results for 77.0 ppm (tolerance: +/-0.5 ppm)
================================================================================
Found 5 matches

[1] Chemical Shift: 77.16 ppm (diff: +/-0.16)
    Nucleus: C-13
    Solvent: CDCl3
    Environment: bonded to C (3H nearby)
    HOSE: C(*C,C,C,H,H,H/)...
```

### 5. `quick_search.py`
**Purpose**: Python API wrapper for quick searches

**Usage**:
```python
from quick_search import search_nmr

result = search_nmr(26.3, tolerance=1.0, max_results=3)
print(result)
```

### 6. `example_search.py`
**Purpose**: Demonstration of search functionality with examples

### 7. `check_boron.py`
**Purpose**: Find all Boron-containing entries in the database

**Usage**:
```bash
# Search entire database
python check_boron.py

# Search first 100k entries
python check_boron.py 100000
```

**Results**: Found 5,740 entries containing Boron

### 8. `find_elements.py`
**Purpose**: Catalog all chemical elements present in the database

**Usage**:
```bash
# Scan entire database
python find_elements.py

# Scan first 100k entries
python find_elements.py 100000
```

**Results**: Found 30 unique elements

#### Element Distribution (Top 10):
| Element | Count | Percentage |
|---------|-------|------------|
| C | 1,720,268 | 28.02% |
| H | 1,696,077 | 27.63% |
| O | 1,312,504 | 21.38% |
| N | 807,243 | 13.15% |
| S | 208,699 | 3.40% |
| Cl | 151,179 | 2.46% |
| F | 105,063 | 1.71% |
| Br | 53,823 | 0.88% |
| P | 51,634 | 0.84% |
| Si | 18,082 | 0.29% |

#### Elements by Category:

**Organic (CHNOPS)**:
- C, H, N, O, P, S

**Halogens**:
- F, Cl, Br, I

**Metalloids**:
- B, Si, Ge, As, Sb, Te, Se

**Transition Metals**:
- Fe, Co, Ni, Cu, Zn, Ag, Au, Pt, Mn, Cr, Mo, V, Ti, Zr, Sc

**Post-transition Metals**:
- Al, Ga, Sn, Tl, Pb, Bi

### 9. `find_overlaps.py`
**Purpose**: Find entries with 100% overlapping chemical shift ranges

These are problematic because we cannot distinguish between them based on shift alone. One range completely contains another (or they're identical).

**Usage**:
```bash
# Search for exact overlaps
python find_overlaps.py 50000 0.0

# With tolerance
python find_overlaps.py 50000 0.5
```

**Results** (50,000 entries analyzed):
- **Total overlapping pairs**: 99,914
- **Number of overlap groups**: 26
- **Largest group**: 2,005 entries with overlapping ranges

**Algorithm**:
- Compares all pairs of entries
- Checks if one range completely contains another
- Groups overlapping entries using BFS (Breadth-First Search)
- Reports largest overlap groups

**Key Function**:
```python
def ranges_overlap_100_percent(range1, range2):
    """Check if one range completely contains the other"""
    min1, max1 = range1
    min2, max2 = range2

    # Check if range1 completely contains range2
    if min1 <= min2 and max1 >= max2:
        return True

    # Check if range2 completely contains range1
    if min2 <= min1 and max2 >= max1:
        return True

    return False
```

**Implications**: The high number of overlaps (99,914 pairs in just 50,000 entries) indicates that single chemical shift values are often insufficient for unambiguous compound identification. Multi-peak matching or HOSE code pattern matching is necessary for reliable compound identification.

## Technical Details

### HOSE Code Structure

HOSE codes encode the chemical environment in spheres (distance layers) separated by `/`:

```
C(=CC/HC,HHH/HHC)
│  ││ │      └── Sphere 3: atoms 3 bonds away
│  ││ └──────── Sphere 2: atoms 2 bonds away
│  │└────────── Sphere 1: atoms 1 bond away
│  └─────────── Bond types: = (double), # (triple), * (aromatic)
└────────────── Central atom
```

**Special Characters**:
- `=` - Double bond
- `#` - Triple bond
- `*` - Aromatic bond (delocalized electrons in aromatic rings like benzene)
- `@` - Ring indicator
- `,` - Separates branches within same sphere
- `/` - Separates spheres (distance layers)

**Example**: `*C*C` represents two aromatic carbons connected by aromatic bonds (part of an aromatic ring)

### Parsing Notes

The database uses **US number format**:
- `.` (period) for decimal separator
- `_` (underscore) for field separator
- NOT German format (which would use `,` for decimals)

Verified by examining raw bytes: `6 7 . 0 _ 6 7 . 0`

## Common Workflows

### Searching for a Compound by Peak

If you have an unknown compound with a peak at 77.0 ppm:

```bash
python search_cli.py 77.0 0.5 10
```

This searches with ±0.5 ppm tolerance and returns top 10 matches.

### Searching with Multiple Peaks

If you have multiple peaks (e.g., 77.0, 26.3, 130.5 ppm):

```python
from nmr_search import load_database, search_by_peak_list

db = load_database(max_entries=100000)
peaks = [77.0, 26.3, 130.5]
results = search_by_peak_list(db, peaks, tolerance=0.5, min_matches=2)

for compound in results:
    print(f"Match: {compound['structure']}")
    print(f"Matched peaks: {compound['matched_peak_count']}")
```

### Finding Specific Elements

To find entries containing a specific element:

```python
# Check for Boron
python check_boron.py

# For other elements, modify check_boron.py
# or use grep on the CSV file
```

## Key Findings

1. **Database Coverage**: 30 elements present, dominated by organic elements (CHNOPS = 94.4%)

2. **Boron Presence**: 5,740 entries contain Boron (0.33% of database)

3. **Overlapping Ranges**: Extensive overlap problem detected
   - 99,914 overlapping pairs found in just 50,000 entries analyzed
   - 26 distinct overlap groups identified
   - Largest group contains 2,005 entries with overlapping ranges
   - Single chemical shift values are often insufficient for unambiguous compound identification
   - Multi-peak matching or HOSE code pattern analysis required for reliable identification

4. **Measurement Uncertainty**: Each peak has min/max values representing measurement variability across multiple observations

5. **HOSE Code Richness**: HOSE codes provide detailed structural information including:
   - Atom types and counts
   - Bond types (single, double, triple, aromatic)
   - Ring structures
   - Branching patterns
   - Distance relationships

6. **Aromatic Compounds**: Aromatic carbonyl carbons (C=O attached to aromatic rings) cluster around 192 ppm with many identical shift values, making them particularly problematic for identification

## Error Fixes Applied

### Unicode Encoding
Fixed Windows compatibility by replacing Unicode symbols:
- `→` changed to `->`
- `±` changed to `+/-`
- `Δ` changed to `diff`

### Recursion Depth
Replaced recursive DFS with iterative BFS in overlap finder to handle large connected components.

## Dependencies

```bash
pip install matplotlib pandas numpy tqdm
```

## Files Generated

- `elements_found.txt` - Complete list of all 30 elements with frequency statistics
- `overlapping_ranges.txt` - Detailed report of 99,914 overlapping pairs grouped into 26 overlap groups

## Future Enhancements

Potential improvements:
1. Web interface for searching
2. Machine learning for compound prediction
3. Integration with structure drawing tools
4. API endpoint for programmatic access
5. Fuzzy matching for HOSE patterns
6. Export to standard NMR formats (JCAMP-DX)

## References

- **HOSE Code**: Bremser, W. (1978). "HOSE - A Novel Substructure Code". Analytica Chimica Acta
- **NMRShiftDB**: Open database for organic structures and their NMR spectra
- **Database Location**: `nmrshiftdb2/nmrshiftdb.csv`

## Contact & Usage

This toolset was created for exploring the NMRShiftDB2 database. All scripts are standalone Python files with minimal dependencies. Feel free to modify and extend for your specific use case.

---

**Last Updated**: 2025-12-17
**Database Version**: nmrshiftdb2
**Total Entries Analyzed**: 1,721,704
