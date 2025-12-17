import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import re
from hose_decoder import decode_hose_code, describe_environment

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
for entry in data[:20000]:  # Parse first 20000 entries for better variety
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

                # Extract HOSE code (after the semicolon)
                hose_code = ""
                if ';' in structure:
                    hose_code = structure.split(';')[1] if len(structure.split(';')) > 1 else ""

                parsed_data.append({
                    'solvent': solvent,
                    'nucleus': nucleus_info,
                    'structure': structure,
                    'hose_code': hose_code,
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

# Select diverse compounds with multiple peaks
interesting_compounds = []
seen_types = set()

for structure, peaks in compounds.items():
    if 4 <= len(peaks) <= 15:  # Compounds with 4-15 peaks
        # Try to get diverse compounds
        if peaks[0]['hose_code']:
            env_desc = describe_environment(peaks[0]['hose_code'])

            # Check for diversity
            env_key = env_desc.split('(')[0].strip()  # Get main description
            if env_key not in seen_types or len(interesting_compounds) < 3:
                interesting_compounds.append((structure, peaks))
                seen_types.add(env_key)

    if len(interesting_compounds) >= 6:
        break

print(f"\nSelected {len(interesting_compounds)} compounds to plot")

# Plot the spectra
fig, axes = plt.subplots(len(interesting_compounds), 1, figsize=(14, 4.5*len(interesting_compounds)))
if len(interesting_compounds) == 1:
    axes = [axes]

for idx, (structure, peaks) in enumerate(interesting_compounds):
    ax = axes[idx]

    # Get all chemical shifts for this compound
    shifts = [p['avg_shift'] for p in peaks]
    intensities = [p['count'] for p in peaks]

    # Normalize intensities
    max_intensity = max(intensities) if intensities else 1
    intensities = [i / max_intensity for i in intensities]

    # Create a simulated NMR spectrum
    ppm_range = np.linspace(min(shifts) - 20, max(shifts) + 20, 2000)
    spectrum = np.zeros_like(ppm_range)

    # Add Lorentzian peaks
    for shift, intensity in zip(shifts, intensities):
        # Lorentzian line shape
        width = 1.5  # Peak width
        spectrum += intensity / (1 + ((ppm_range - shift) / width)**2)

    # Plot
    ax.plot(ppm_range, spectrum, 'b-', linewidth=1.2)
    ax.fill_between(ppm_range, 0, spectrum, alpha=0.3)

    # Mark peak positions
    for shift, intensity in zip(shifts, intensities):
        ax.axvline(x=shift, color='r', linestyle='--', alpha=0.3, linewidth=0.5)

    # Format plot
    ax.set_xlim(max(shifts) + 20, min(shifts) - 20)  # Reverse x-axis (NMR convention)
    ax.set_ylim(0, max(spectrum) * 1.15)
    ax.set_xlabel('Chemical Shift (ppm)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Intensity', fontsize=11)

    # Decode HOSE codes to create meaningful labels
    if peaks[0]['hose_code']:
        # Get primary environment description
        primary_env = describe_environment(peaks[0]['hose_code'])

        # Decode all environments in the compound
        all_envs = []
        for peak in peaks[:5]:  # Check first 5 peaks
            if peak['hose_code']:
                env = describe_environment(peak['hose_code'])
                env_type = env.split('(')[0].strip()
                if env_type not in all_envs:
                    all_envs.append(env_type)

        compound_desc = " | ".join(all_envs[:3])  # Show up to 3 types
    else:
        compound_desc = "Unknown structure"

    # Create title
    nucleus = peaks[0]['nucleus'] if peaks else 'Unknown'
    solvent = peaks[0]['solvent'] if peaks else 'Unknown'

    title_text = f"Compound {idx+1}: {compound_desc}\n{nucleus} NMR in {solvent} | {len(peaks)} peaks"
    ax.set_title(title_text, fontsize=12, fontweight='bold', pad=10)
    ax.grid(True, alpha=0.3, linestyle=':')

    # Add text with key peak information
    shift_range = f"{min(shifts):.1f} - {max(shifts):.1f} ppm"
    info_text = f"Chemical shift range: {shift_range}\nKey peaks: {', '.join([f'{s:.1f}' for s in sorted(shifts, reverse=True)[:5]])}"

    ax.text(0.02, 0.97, info_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8, edgecolor='gray'))

plt.tight_layout()
plt.savefig('nmr_spectra_labeled.png', dpi=300, bbox_inches='tight')
print("\nPlot saved as 'nmr_spectra_labeled.png'")
plt.show()

# Print detailed summary
print("\n" + "="*80)
print("DETAILED COMPOUND ANALYSIS")
print("="*80)

for idx, (structure, peaks) in enumerate(interesting_compounds):
    print(f"\n{'='*80}")
    print(f"COMPOUND {idx+1}")
    print(f"{'='*80}")

    # Analyze all peaks in this compound
    print(f"Total peaks: {len(peaks)}")
    print(f"Nucleus: {peaks[0]['nucleus']}")
    print(f"Solvent: {peaks[0]['solvent']}")
    print(f"Chemical shift range: {min(p['avg_shift'] for p in peaks):.1f} - {max(p['avg_shift'] for p in peaks):.1f} ppm")

    print("\nPeak Analysis:")
    print(f"{'Peak #':<8} {'Shift (ppm)':<15} {'Environment':<50}")
    print("-" * 73)

    for i, peak in enumerate(sorted(peaks, key=lambda x: x['avg_shift'], reverse=True)[:10], 1):
        if peak['hose_code']:
            env_desc = describe_environment(peak['hose_code'])
        else:
            env_desc = "Unknown"
        print(f"{i:<8} {peak['avg_shift']:<15.1f} {env_desc[:48]:<50}")

    if len(peaks) > 10:
        print(f"... and {len(peaks)-10} more peaks")

print("\n" + "="*80)
