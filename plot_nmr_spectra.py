import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import re

# Read the NMR database
print("Reading NMR database...")
data = []
with open('nmrshiftdb2/nmrshiftdb.csv', 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        data.append(line)

print(f"Total entries: {len(data)}")

# Parse the data
parsed_data = []
for entry in data[:10000]:  # Parse first 10000 entries for speed
    try:
        # Split by underscore to get components
        parts = entry.split('_')
        if len(parts) >= 4:
            solvent = parts[0]
            nucleus_info = parts[1]
            structure = parts[2]

            # Extract chemical shifts
            shift_parts = parts[3:]
            shift_str = '_'.join(shift_parts)
            shift_values = shift_str.split('_')

            if len(shift_values) >= 4:
                min_shift = float(shift_values[0])
                max_shift = float(shift_values[1])
                avg_shift = float(shift_values[2])
                count = int(shift_values[3])

                parsed_data.append({
                    'solvent': solvent,
                    'nucleus': nucleus_info,
                    'structure': structure,
                    'min_shift': min_shift,
                    'max_shift': max_shift,
                    'avg_shift': avg_shift,
                    'count': count
                })
    except Exception as e:
        continue

print(f"Successfully parsed: {len(parsed_data)} entries")

# Group spectra by structure (compound)
compounds = defaultdict(list)
for entry in parsed_data:
    compounds[entry['structure']].append(entry)

# Select a few compounds with multiple peaks
interesting_compounds = []
for structure, peaks in compounds.items():
    if 3 <= len(peaks) <= 20:  # Compounds with 3-20 peaks
        interesting_compounds.append((structure, peaks))
    if len(interesting_compounds) >= 5:
        break

print(f"\nSelected {len(interesting_compounds)} compounds to plot")

# Plot the spectra
fig, axes = plt.subplots(len(interesting_compounds), 1, figsize=(12, 4*len(interesting_compounds)))
if len(interesting_compounds) == 1:
    axes = [axes]

for idx, (structure, peaks) in enumerate(interesting_compounds):
    ax = axes[idx]

    # Get all chemical shifts for this compound
    shifts = [p['avg_shift'] for p in peaks]
    intensities = [p['count'] for p in peaks]

    # Normalize intensities
    max_intensity = max(intensities)
    intensities = [i / max_intensity for i in intensities]

    # Create a simulated NMR spectrum
    ppm_range = np.linspace(min(shifts) - 20, max(shifts) + 20, 2000)
    spectrum = np.zeros_like(ppm_range)

    # Add Lorentzian peaks
    for shift, intensity in zip(shifts, intensities):
        # Lorentzian line shape
        width = 2.0  # Peak width
        spectrum += intensity / (1 + ((ppm_range - shift) / width)**2)

    # Plot
    ax.plot(ppm_range, spectrum, 'b-', linewidth=1)
    ax.fill_between(ppm_range, 0, spectrum, alpha=0.3)

    # Mark peak positions
    for shift, intensity in zip(shifts, intensities):
        ax.axvline(x=shift, color='r', linestyle='--', alpha=0.3, linewidth=0.5)

    # Format plot
    ax.set_xlim(max(shifts) + 20, min(shifts) - 20)  # Reverse x-axis (NMR convention)
    ax.set_ylim(0, max(spectrum) * 1.1)
    ax.set_xlabel('Chemical Shift (ppm)', fontsize=10)
    ax.set_ylabel('Intensity', fontsize=10)

    # Create compound label
    nucleus = peaks[0]['nucleus'] if peaks else 'Unknown'
    solvent = peaks[0]['solvent'] if peaks else 'Unknown'
    compound_label = f"{nucleus} NMR\nSolvent: {solvent}\nPeaks: {len(peaks)}"

    ax.set_title(compound_label, fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Add text with peak positions
    peak_text = f"Chemical shifts: {', '.join([f'{s:.1f}' for s in sorted(shifts, reverse=True)])}"
    ax.text(0.02, 0.95, peak_text, transform=ax.transAxes,
            fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('nmr_spectra.png', dpi=300, bbox_inches='tight')
print("\nPlot saved as 'nmr_spectra.png'")
plt.show()

# Print summary
print("\n" + "="*60)
print("SUMMARY OF PLOTTED COMPOUNDS")
print("="*60)
for idx, (structure, peaks) in enumerate(interesting_compounds):
    print(f"\nCompound {idx+1}:")
    print(f"  Structure code: {structure[:50]}...")
    print(f"  Number of peaks: {len(peaks)}")
    print(f"  Nucleus: {peaks[0]['nucleus']}")
    print(f"  Solvent: {peaks[0]['solvent']}")
    print(f"  Chemical shift range: {min(p['avg_shift'] for p in peaks):.1f} - {max(p['avg_shift'] for p in peaks):.1f} ppm")
