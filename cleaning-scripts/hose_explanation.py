"""
HOSE Code Sphere Explanation
"""

def explain_hose_spheres(hose_code):
    """Break down a HOSE code sphere by sphere"""

    print(f"\nHOSE Code: {hose_code}")
    print("="*70)

    if '/' in hose_code:
        spheres = hose_code.split('/')
    else:
        spheres = [hose_code]

    print(f"\nThis HOSE code has {len(spheres)} spheres\n")

    for i, sphere in enumerate(spheres):
        print(f"Sphere {i} (distance {i} bonds from central atom):")
        print(f"  Raw content: {sphere}")

        # Count atoms in this sphere
        atoms = {'C': 0, 'H': 0, 'O': 0, 'N': 0, 'S': 0}
        special = []

        for char in sphere:
            if char in atoms:
                atoms[char] += 1
            elif char == '=':
                special.append('double bond')
            elif char == '#':
                special.append('triple bond')
            elif char == '@':
                special.append('ring closure')
            elif char == '&':
                special.append('heteroatom')

        # Print atom counts
        atom_list = [f"{count}{atom}" for atom, count in atoms.items() if count > 0]
        if atom_list:
            print(f"  Atoms: {', '.join(atom_list)}")
        if special:
            print(f"  Special features: {', '.join(set(special))}")
        print()

# Example HOSE codes with different sphere counts
examples = [
    # Simple example - methyl group
    ("HHHC", "Methyl carbon (CH3)"),

    # Two spheres - ethyl group
    ("HHHC(HHC/", "Ethyl carbon (CH3-CH2-)"),

    # Three spheres - propyl group
    ("HHHC(HHC/HHC/", "Propyl chain (CH3-CH2-CH2-)"),

    # Aromatic example with 3 spheres
    ("H=CC(HC,=CC/H=C,\\H|C,HHH/", "Aromatic carbon in benzene ring"),

    # Oxygen-containing with rings
    ("@HCOC(@HOC,H,@OCH/C,HHO,H,@HO&/@H&O,H,H)", "Sugar-like structure with ring"),
]

print("HOSE CODE SPHERE BREAKDOWN")
print("="*70)
print("\nThe '/' character separates spheres (distances from central atom)")
print("Sphere 0 = central atom and its direct neighbors")
print("Sphere 1 = atoms 1 bond away")
print("Sphere 2 = atoms 2 bonds away")
print("... and so on\n")

for hose_code, description in examples:
    print("\n" + "="*70)
    print(f"Example: {description}")
    explain_hose_spheres(hose_code)

# Visual example
print("\n" + "="*70)
print("VISUAL EXAMPLE: Propane (CH3-CH2-CH3)")
print("="*70)
print("""
Looking at the MIDDLE carbon (CH2):

         H   H   H
         |   |   |
     H - C - C - C - H
         |   |   |
         H   H   H
             ^
             |
        Central atom

Sphere 0: The central carbon and what's directly attached
          Content: HHC(
          2 Hydrogens (HH) and 1 Carbon (C) directly bonded

Sphere 1: What's attached to those neighbors (1 bond away)
          Content: HHC,HHC/
          Two CH3 groups (HHC and HHC) separated by comma

So the HOSE code for the middle carbon: HHC(HHC,HHC/

The '/' at the end indicates we could continue to sphere 2,
but for propane there's nothing further.
""")
