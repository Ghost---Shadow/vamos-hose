"""
Tests for HOSE-to-SMILES converter.
Validates the tokenizer, parser, tree builder, and SMILES serializer.
"""

from hose_to_smiles import (
    tokenize,
    parse_hose,
    hose_to_smiles,
    extract_central_atom,
    tree_to_smiles,
    _find_ring_pairs,
)
from preprocess_database import parse_line


def test_tokenizer():
    """Test that the tokenizer produces correct token sequences."""
    tokens = tokenize("HHHC")
    assert len(tokens) == 4
    assert all(t.type == "ATOM" for t in tokens)
    assert [t.value for t in tokens] == ["H", "H", "H", "C"]

    tokens = tokenize("=CC")
    assert tokens[0].type == "DOUBLE"
    assert tokens[1].type == "ATOM" and tokens[1].value == "C"
    assert tokens[2].type == "ATOM" and tokens[2].value == "C"

    tokens = tokenize("*C*C")
    assert tokens[0].type == "AROMATIC"
    assert tokens[1].value == "C"
    assert tokens[2].type == "AROMATIC"
    assert tokens[3].value == "C"

    # Bremser substitutions
    tokens = tokenize("XYQ")
    assert [t.value for t in tokens] == ["Cl", "Br", "Si"]

    # Charges
    tokens = tokenize("N+")
    assert tokens[0].type == "ATOM" and tokens[0].value == "N"
    assert tokens[1].type == "CHARGE_POS"

    # Ring closure
    tokens = tokenize("@HC")
    assert tokens[0].type == "RING"
    assert tokens[1].value == "H"
    assert tokens[2].value == "C"

    print("  tokenizer: PASS")


def test_tree_simple():
    """Test tree building for simple HOSE codes."""
    # HHHC: central C bonded to H, H, H, C
    root = parse_hose("HHHC", "C")
    assert root.atom == "C"
    assert len(root.children) == 4
    atoms = sorted([c.atom for c in root.children])
    assert atoms == ["C", "H", "H", "H"]

    # C(=CC/...): central atom with 1 neighbor C, which has children =C and C
    root = parse_hose("C(=CC/HC,HHH/HHC)", "C")
    assert root.atom == "C"
    assert len(root.children) == 1  # one sphere-0 neighbor: C
    c_child = root.children[0]
    assert c_child.atom == "C"
    # Sphere 1: =C and C (two children of the C neighbor)
    assert len(c_child.children) == 2
    bonds = sorted([c.bond for c in c_child.children])
    assert "=" in bonds

    print("  tree (simple): PASS")


def test_smiles_output():
    """Test that SMILES output is non-empty and contains expected atoms."""
    # Simple methyl: HHHC -> CC (implicit H, two carbons)
    s = hose_to_smiles("HHHC", "C")
    assert s is not None
    assert "C" in s
    assert s.count("C") >= 2  # central C + neighbor C
    # H atoms should be implicit (suppressed)
    assert "[H]" not in s

    # Aromatic
    s = hose_to_smiles("*C*C", "C")
    assert s is not None
    assert "c" in s  # lowercase aromatic c

    # Double bond
    s = hose_to_smiles("=CC", "C")
    assert s is not None
    assert "=" in s

    # Halogens via Bremser
    s = hose_to_smiles("XC", "C")
    assert s is not None
    assert "Cl" in s

    s = hose_to_smiles("YC", "C")
    assert s is not None
    assert "Br" in s

    # Charged nitrogen
    s = hose_to_smiles("N+C", "C")
    assert s is not None
    assert "N+" in s or "[N+]" in s

    print("  smiles output: PASS")


def test_empty_and_edge_cases():
    """Test edge cases."""
    # Empty HOSE code
    assert hose_to_smiles("", "C") is None
    assert hose_to_smiles(None, "C") is None

    # Truncated (missing closing paren)
    s = hose_to_smiles("HHHC(HCC/H=C,HHH/", "C")
    assert s is not None  # should not crash

    # Double commas (empty branches)
    s = hose_to_smiles("CC(HH,,HH/", "C")
    assert s is not None

    # Trailing slash
    s = hose_to_smiles("HHHC/", "C")
    assert s is not None

    print("  edge cases: PASS")


def test_nucleus_extraction():
    """Test central atom extraction from nucleus strings."""
    assert extract_central_atom("C-4") == "C"
    assert extract_central_atom("C-3-6") == "C"
    assert extract_central_atom("H-1") == "H"
    assert extract_central_atom("N-15") == "N"
    assert extract_central_atom("") == "C"
    assert extract_central_atom(None) == "C"

    print("  nucleus extraction: PASS")


def test_line_parsing():
    """Test database line parsing."""
    line = "C_Pyridin-D5 (C5D5N) + Methanol-D4 (CD3OD)_C-4-6;@HCOC(HHC,H,HHC/@CCC,=&C/=&C,H=C,HHH,HHH),@HCO,\\H|C/_67.0_67.0_67.0_1"
    entry = parse_line(line)
    assert entry is not None
    assert entry["nucleus"] == "C-4-6"
    assert entry["solvent"] == "C_Pyridin-D5 (C5D5N) + Methanol-D4 (CD3OD)"
    assert entry["min_shift"] == 67.0
    assert entry["max_shift"] == 67.0
    assert entry["avg_shift"] == 67.0
    assert entry["count"] == 1
    assert entry["hose_code"].startswith("@HCOC")

    # Simple solvent
    line2 = "C_Unreported_C-3;H=CC(HC,HHC/H=C,@OCC/_131.4_131.4_131.4_1"
    entry2 = parse_line(line2)
    assert entry2 is not None
    assert entry2["nucleus"] == "C-3"
    assert entry2["solvent"] == "C_Unreported"

    # CDCl3 solvent
    line3 = "C_Chloroform-D1 (CDCl3)_C-4;CC(HHC,%C/HHC,C/_18.1607_18.665174999999998_19.2_4"
    entry3 = parse_line(line3)
    assert entry3 is not None
    assert entry3["nucleus"] == "C-4"
    assert abs(entry3["avg_shift"] - 19.2) < 0.001

    print("  line parsing: PASS")


def test_real_database_entries():
    """Test conversion of actual database entries."""
    test_lines = [
        "C_Pyridin-D5 (C5D5N) + Methanol-D4 (CD3OD)_C-4;HHHC(HCC/H=C,HHH/\\H|C)@HCO/@HCO,C/_21.3_21.3_21.3_1",
        "C_Unreported_C-3-6;*C*CC(*C,*C,C,C,=OO/H,H,*C,*&,HHH,HHH,,H/*&C)_129.2_129.2_129.2_1",
        "C_Unreported_C-4;HHHC(HHN/CC/HH,HHC,C)HCO,HHH/HHO,H/_11.77_11.77_11.77_1",
        "C_Chloroform-D1 (CDCl3)_C-4;CC(HHC,%C/HHC,C/_18.1607_18.665174999999998_19.2_4",
    ]

    for line in test_lines:
        entry = parse_line(line)
        assert entry is not None, f"Failed to parse: {line[:60]}"
        atom = entry["nucleus"].split("-")[0]
        smiles = hose_to_smiles(entry["hose_code"], atom)
        assert smiles is not None, f"Failed to convert HOSE: {entry['hose_code'][:40]}"
        assert len(smiles) > 0

    print("  real database entries: PASS")


if __name__ == "__main__":
    print("HOSE-to-SMILES Converter Tests")
    print("=" * 50)
    test_tokenizer()
    test_tree_simple()
    test_smiles_output()
    test_empty_and_edge_cases()
    test_nucleus_extraction()
    test_line_parsing()
    test_real_database_entries()
    print()
    print("ALL TESTS PASSED")
