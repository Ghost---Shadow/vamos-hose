"""
HOSE Code Decoder
HOSE = Hierarchically Ordered Spherical Environment
Describes the chemical environment around a central atom in spheres
"""


def decode_hose_code(hose_code):
    """
    Decode a HOSE code to extract structural information

    HOSE code format: Each sphere (distance from central atom) is separated by /
    Within each sphere, atoms are listed with their properties:
    - H = Hydrogen
    - C = Carbon
    - O = Oxygen
    - N = Nitrogen
    - = indicates double bond
    - # indicates triple bond
    - @ indicates ring closure
    - Numbers indicate connectivity
    """

    info = {
        "spheres": [],
        "atoms": {
            "C": 0,
            "H": 0,
            "O": 0,
            "N": 0,
            "S": 0,
            "P": 0,
            "F": 0,
            "Cl": 0,
            "Br": 0,
            "I": 0,
        },
        "double_bonds": 0,
        "triple_bonds": 0,
        "aromatic": False,
        "ring": False,
    }

    # Split by spheres (separated by /)
    if "/" in hose_code:
        spheres = hose_code.split("/")
    else:
        spheres = [hose_code]

    info["num_spheres"] = len(spheres)

    for sphere_idx, sphere in enumerate(spheres):
        sphere_info = {"level": sphere_idx + 1, "content": sphere, "atoms": []}

        # Parse the sphere content
        i = 0
        while i < len(sphere):
            char = sphere[i]

            # Check for atom symbols
            if char in ["C", "H", "O", "N", "S", "P", "F", "B"]:
                info["atoms"][char] += 1
                sphere_info["atoms"].append(char)
            # Check for bond types
            elif char == "=":
                info["double_bonds"] += 1
                if i + 1 < len(sphere) and sphere[i + 1] == "C":
                    # Could be aromatic or alkene
                    info["aromatic"] = True
            elif char == "#":
                info["triple_bonds"] += 1
            elif char == "@":
                info["ring"] = True

            i += 1

        info["spheres"].append(sphere_info)

    return info


def describe_environment(hose_code):
    """Generate a human-readable description of the HOSE code"""
    info = decode_hose_code(hose_code)

    descriptions = []

    # Count total heavy atoms (non-H)
    heavy_atoms = sum(
        count for atom, count in info["atoms"].items() if atom != "H" and count > 0
    )

    # Describe the environment
    if info["atoms"]["O"] > 0:
        if info["double_bonds"] > 0:
            descriptions.append("carbonyl/carboxyl")
        else:
            descriptions.append("alcohol/ether")

    if info["aromatic"] or (info["double_bonds"] >= 2 and info["atoms"]["C"] > 4):
        descriptions.append("aromatic")
    elif info["double_bonds"] > 0:
        descriptions.append("alkene")

    if info["triple_bonds"] > 0:
        descriptions.append("alkyne")

    if info["ring"]:
        descriptions.append("ring")

    if info["atoms"]["N"] > 0:
        descriptions.append("amine/amide")

    if info["atoms"]["S"] > 0:
        descriptions.append("sulfur")

    # If no special features, classify by saturation
    if not descriptions:
        if info["atoms"]["H"] > 2:
            descriptions.append("aliphatic")
        elif info["atoms"]["C"] > 0:
            descriptions.append("quaternary")

    # Build description string
    if descriptions:
        env_desc = " + ".join(descriptions)
    else:
        env_desc = "hydrocarbon"

    # Add atom count info
    atom_summary = []
    for atom in ["C", "O", "N", "S"]:
        if info["atoms"][atom] > 0:
            atom_summary.append(f"{info['atoms'][atom]}{atom}")

    if atom_summary:
        env_desc += f" ({', '.join(atom_summary)} nearby)"

    return env_desc


# Test the decoder
if __name__ == "__main__":
    test_codes = [
        "CC(HHC,H=C/=CC,CC/HC,HHH,HHH,HHH)H=C/HC/",
        "@HCOC(@HOC,H,@OCH/C,HHO,H,@HO&/@H&O,H,H)",
        "H=CC(HC,=CC/H=C,\\H|C,HHH/",
        "=OCC(,@CCH,=CC/HHH,@HCC,|H\\&,HHC/@H&O,HCC,H=C)",
        "HHHC(HCC/H=C,HHH/\\H|C)",
    ]

    print("HOSE Code Decoder Test\n" + "=" * 60)
    for code in test_codes:
        print(f"\nHOSE Code: {code[:50]}...")
        info = decode_hose_code(code)
        print(f"Spheres: {info['num_spheres']}")
        print(
            f"Atoms: C={info['atoms']['C']}, H={info['atoms']['H']}, O={info['atoms']['O']}, N={info['atoms']['N']}"
        )
        print(
            f"Features: Double bonds={info['double_bonds']}, Ring={info['ring']}, Aromatic={info['aromatic']}"
        )
        print(f"Description: {describe_environment(code)}")
