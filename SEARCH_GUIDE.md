# NMR Database Search Guide

## Overview
This toolkit allows you to search the NMR database by chemical shift values to identify potential structures.

## Files
- `nmr_search.py` - Main search tool (interactive)
- `example_search.py` - Example searches with demonstration
- `hose_decoder.py` - Decodes HOSE codes to structural descriptions

## Usage

### 1. Interactive Search

Run the interactive search tool:
```bash
python nmr_search.py
```

You'll be prompted to choose:
- **Option 1**: Search by single chemical shift
- **Option 2**: Search by list of peaks (find matching compounds)

#### Example Session (Option 1 - Single Peak):
```
Enter choice (1 or 2): 1
Enter chemical shift (ppm): 77.5
Enter tolerance (ppm) [default: 0.5]: 1.0
```

This will find all peaks near 77.5 ppm (±1.0 ppm) and show:
- Chemical shift values
- Solvent and nucleus type
- Environment description (decoded from HOSE code)
- HOSE code snippet

#### Example Session (Option 2 - Multiple Peaks):
```
Enter choice (1 or 2): 2
Enter your peak list (comma-separated): 128.5, 77.2, 42.1, 25.3
Enter tolerance (ppm) [default: 1.0]: 2.0
Minimum matches required [default: 3]: 2
```

This will find compounds that match your peak list and show:
- Total peaks in the compound
- How many peaks matched
- Match ratio (percentage)
- All chemical shifts for that compound
- Environment descriptions

### 2. Run Examples

See pre-configured search examples:
```bash
python example_search.py
```

This demonstrates:
- Searching for alcohol/ether carbons (~77 ppm)
- Searching for aromatic carbons (~128 ppm)
- Compound matching with multiple peaks
- Searching for aliphatic methyls (~20 ppm)

## Understanding Results

### Chemical Shift Ranges (¹³C NMR):
- **0-50 ppm**: Aliphatic carbons (sp³)
  - 10-30 ppm: Methyl groups (CH₃)
  - 20-40 ppm: Methylene groups (CH₂)
  - 30-50 ppm: Methine groups (CH)

- **50-90 ppm**: Carbons bonded to oxygen/nitrogen
  - 50-65 ppm: C-N (amines)
  - 60-90 ppm: C-O (alcohols, ethers)

- **100-150 ppm**: Aromatic/sp² carbons
  - 110-140 ppm: Aromatic C-H
  - 140-165 ppm: Aromatic C-O, C-N

- **160-220 ppm**: Carbonyl carbons
  - 160-185 ppm: Carboxylic acids, esters, amides
  - 185-220 ppm: Ketones, aldehydes

### HOSE Code Features:
The environment descriptions show:
- **aromatic**: Benzene rings or aromatic systems
- **alcohol/ether**: C-O bonds
- **carbonyl/carboxyl**: C=O groups
- **amine/amide**: C-N bonds
- **alkene**: C=C double bonds
- **ring**: Cyclic structures
- **aliphatic**: Saturated alkyl chains

## Tips for Effective Searching

1. **Start with wider tolerance** (±2.0 ppm) then narrow down
2. **For compound matching**, require at least 50% of peaks to match
3. **Use chemical shift regions** to filter by functional group
4. **Compare HOSE descriptions** to understand structural environment
5. **Consider solvent effects** - different solvents can shift peaks by ~0.5-2 ppm

## Python API Usage

You can also use the search functions programmatically:

```python
from nmr_search import load_database, search_by_shift, search_by_peak_list

# Load database
database = load_database(max_entries=100000)

# Search for single peak
results = search_by_shift(database, target_shift=77.0, tolerance=1.0)

# Search for compound by multiple peaks
peak_list = [128.5, 77.2, 42.1, 25.3]
matches = search_by_peak_list(database, peak_list, tolerance=2.0, min_matches=2)

# Display results
from nmr_search import display_search_results, display_compound_matches
display_search_results(results)
display_compound_matches(matches)
```

## Example Queries

### Find aromatic compounds:
- Search for peaks in 120-140 ppm range
- Look for "aromatic" in environment description

### Find sugars/carbohydrates:
- Multiple peaks in 60-90 ppm range
- Look for "alcohol/ether + ring" in environment

### Find aliphatic chains:
- Multiple peaks in 10-40 ppm range
- Look for "aliphatic" in environment

### Find carbonyl compounds:
- Peaks in 160-220 ppm range
- Look for "carbonyl/carboxyl" in environment
