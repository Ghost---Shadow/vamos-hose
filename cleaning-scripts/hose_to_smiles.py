"""
HOSE Code to SMILES Fragment Converter

Converts HOSE (Hierarchically Ordered Spherical Environment) codes
into SMILES fragment strings representing the local chemical environment.

HOSE code structure:
  - Atoms before '(' = sphere 0 (direct neighbors of central atom)
  - Content inside '(...)' = spheres 1,2,3... for the LAST sphere-0 atom
  - Content after ')' = spheres 1,2,3... for the REMAINING sphere-0 atoms
  - '/' separates spheres, ',' separates branches within a sphere
"""

# Bremser single-letter element substitutions used in HOSE codes
BREMSER_MAP = {
    "X": "Cl",
    "Y": "Br",
    "Q": "Si",
    "G": "Ge",
    "T": "Te",
    "A": "As",
    "L": "Tl",
    "M": "Sn",
    "K": "Sb",
}

# Atoms that appear in HOSE codes (single uppercase letter)
HOSE_ATOMS = set("CHONSPFBIXYQGTALMKZRVWUDEe")

# SMILES organic subset (atoms that don't need brackets when uncharged)
ORGANIC_SUBSET = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}

# Token types
TOK_ATOM = "ATOM"
TOK_DOUBLE = "DOUBLE"
TOK_AROMATIC = "AROMATIC"
TOK_TRIPLE = "TRIPLE"
TOK_RING = "RING"
TOK_SPHERE_SEP = "SPHERE_SEP"
TOK_BRANCH_SEP = "BRANCH_SEP"
TOK_OPEN = "OPEN"
TOK_CLOSE = "CLOSE"
TOK_CHARGE_POS = "CHARGE_POS"
TOK_CHARGE_NEG = "CHARGE_NEG"
TOK_DELOCALIZED = "DELOCALIZED"
TOK_STEREO = "STEREO"


class Token:
    __slots__ = ("type", "value")

    def __init__(self, type, value=None):
        self.type = type
        self.value = value

    def __repr__(self):
        if self.value:
            return f"Token({self.type}, {self.value})"
        return f"Token({self.type})"


class TreeNode:
    __slots__ = (
        "atom",
        "bond",
        "children",
        "sphere",
        "is_ring_atom",
        "charge",
        "is_aromatic",
        "parent",
        "_ring_close_to",
    )

    def __init__(self, atom, bond="", sphere=0):
        self.atom = atom
        self.bond = bond  # "", "=", "#", "aromatic"
        self.children = []
        self.sphere = sphere
        self.is_ring_atom = False
        self.charge = 0
        self.is_aromatic = False
        self.parent = None
        self._ring_close_to = None  # set when this node closes a ring

    def add_child(self, child):
        child.parent = self
        self.children.append(child)


def tokenize(hose_code):
    """Convert a HOSE code string into a list of tokens."""
    tokens = []
    i = 0
    n = len(hose_code)

    while i < n:
        ch = hose_code[i]

        if ch == "(":
            tokens.append(Token(TOK_OPEN))
        elif ch == ")":
            tokens.append(Token(TOK_CLOSE))
        elif ch == "/":
            tokens.append(Token(TOK_SPHERE_SEP))
        elif ch == ",":
            tokens.append(Token(TOK_BRANCH_SEP))
        elif ch == "=":
            tokens.append(Token(TOK_DOUBLE))
        elif ch == "*":
            tokens.append(Token(TOK_AROMATIC))
        elif ch == "%":
            tokens.append(Token(TOK_TRIPLE))
        elif ch == "#":
            tokens.append(Token(TOK_TRIPLE))
        elif ch == "@":
            tokens.append(Token(TOK_RING))
        elif ch == "&":
            tokens.append(Token(TOK_DELOCALIZED))
        elif ch == "+":
            tokens.append(Token(TOK_CHARGE_POS))
        elif ch == "-":
            tokens.append(Token(TOK_CHARGE_NEG))
        elif ch == "|" or ch == "\\":
            tokens.append(Token(TOK_STEREO))
        elif ch.isupper():
            resolved = BREMSER_MAP.get(ch, ch)
            tokens.append(Token(TOK_ATOM, resolved))
        # Skip anything else (lowercase, digits, whitespace, unknown)

        i += 1

    return tokens


def _is_bond_or_modifier(tok):
    return tok.type in (
        TOK_DOUBLE,
        TOK_AROMATIC,
        TOK_TRIPLE,
        TOK_RING,
        TOK_DELOCALIZED,
        TOK_STEREO,
        TOK_CHARGE_POS,
        TOK_CHARGE_NEG,
    )


def _parse_one_atom(tokens, pos):
    """Parse a single atom with its bond prefix, ring marker, charge, etc.
    Returns (atom_symbol, bond_type, is_ring, charge, is_aromatic, new_pos)
    or (None, ..., new_pos) if no atom found at this position."""
    bond = ""
    is_ring = False
    charge = 0
    is_aromatic = False
    n = len(tokens)

    # Consume prefixes
    while pos < n and _is_bond_or_modifier(tokens[pos]):
        t = tokens[pos]
        if t.type == TOK_DOUBLE:
            bond = "="
        elif t.type == TOK_AROMATIC:
            is_aromatic = True
            bond = "aromatic"
        elif t.type == TOK_TRIPLE:
            bond = "#"
        elif t.type == TOK_RING:
            is_ring = True
        elif t.type == TOK_CHARGE_POS:
            charge = 1
        elif t.type == TOK_CHARGE_NEG:
            charge = -1
        elif t.type == TOK_DELOCALIZED:
            pass  # informational, skip
        elif t.type == TOK_STEREO:
            pass  # skip stereo markers
        pos += 1

    # Now expect an ATOM token
    if pos < n and tokens[pos].type == TOK_ATOM:
        symbol = tokens[pos].value
        pos += 1
        # Check for trailing charge after atom
        if pos < n and tokens[pos].type == TOK_CHARGE_POS:
            charge = 1
            pos += 1
        elif pos < n and tokens[pos].type == TOK_CHARGE_NEG:
            charge = -1
            pos += 1
        # Check for trailing aromatic marker (e.g., *C* pattern)
        if pos < n and tokens[pos].type == TOK_AROMATIC:
            is_aromatic = True
            if bond != "aromatic":
                bond = "aromatic"
            pos += 1
        # Check for trailing delocalization
        if pos < n and tokens[pos].type == TOK_DELOCALIZED:
            pos += 1
        return symbol, bond, is_ring, charge, is_aromatic, pos

    return None, bond, is_ring, charge, is_aromatic, pos


def _parse_sphere_block(tokens, pos, stop_types):
    """Parse sphere content until we hit a token type in stop_types or end of tokens.
    Returns (list_of_spheres, new_pos) where each sphere is a list of branches,
    and each branch is a list of (symbol, bond, is_ring, charge, is_aromatic) tuples."""
    spheres = []
    current_sphere = []
    current_branch = []
    n = len(tokens)

    while pos < n:
        if tokens[pos].type in stop_types:
            break

        if tokens[pos].type == TOK_SPHERE_SEP:
            current_sphere.append(current_branch)
            spheres.append(current_sphere)
            current_sphere = []
            current_branch = []
            pos += 1
        elif tokens[pos].type == TOK_BRANCH_SEP:
            current_sphere.append(current_branch)
            current_branch = []
            pos += 1
        else:
            sym, bond, is_ring, charge, is_arom, pos = _parse_one_atom(tokens, pos)
            if sym is not None:
                current_branch.append((sym, bond, is_ring, charge, is_arom))
            # If sym is None, we consumed modifiers but found no atom - skip

    # Don't forget the last branch/sphere
    if current_branch:
        current_sphere.append(current_branch)
    if current_sphere:
        spheres.append(current_sphere)

    return spheres, pos


def parse_hose(hose_code, central_atom="C"):
    """Parse a HOSE code into an atom tree rooted at the central atom.
    Returns the root TreeNode."""
    tokens = tokenize(hose_code)
    if not tokens:
        root = TreeNode(central_atom, sphere=0)
        return root

    root = TreeNode(central_atom, sphere=0)
    pos = 0
    n = len(tokens)

    # Phase 1: Parse sphere 0 atoms (before the first OPEN_PAREN)
    sphere0_atoms = []
    while pos < n and tokens[pos].type != TOK_OPEN:
        # Stop if we hit a sphere separator (shouldn't happen before '(' but be safe)
        if tokens[pos].type == TOK_SPHERE_SEP:
            pos += 1
            break
        sym, bond, is_ring, charge, is_arom, pos = _parse_one_atom(tokens, pos)
        if sym is not None:
            sphere0_atoms.append((sym, bond, is_ring, charge, is_arom))

    # Detect sphere-0 ring: if ANY sphere-0 atom has is_ring=True,
    # the non-H atoms form a chain (ring path) rather than siblings.
    has_ring = any(r for _, _, r, _, _ in sphere0_atoms)

    if has_ring:
        # Separate H atoms (branches) from non-H atoms (ring chain)
        h_atoms = [(s, b, r, c, a) for s, b, r, c, a in sphere0_atoms if s == "H"]
        chain_atoms = [(s, b, r, c, a) for s, b, r, c, a in sphere0_atoms if s != "H"]

        # H atoms are normal branches off root
        for sym, bond, is_ring, charge, is_arom in h_atoms:
            child = TreeNode(sym, bond=bond, sphere=1)
            child.is_ring_atom = False
            child.charge = charge
            child.is_aromatic = is_arom
            root.add_child(child)

        # Non-H atoms form a chain: root → first → second → ... → last
        # The last one closes back to root via ring closure
        if chain_atoms:
            prev = root
            for i, (sym, bond, is_ring, charge, is_arom) in enumerate(chain_atoms):
                child = TreeNode(sym, bond=bond, sphere=1)
                child.is_ring_atom = is_ring
                child.charge = charge
                child.is_aromatic = is_arom
                prev.add_child(child)
                prev = child
            # Mark the last chain atom for ring closure back to root
            prev._ring_close_to = root
    else:
        # Normal case: all sphere-0 atoms are siblings of root
        for sym, bond, is_ring, charge, is_arom in sphere0_atoms:
            child = TreeNode(sym, bond=bond, sphere=1)
            child.is_ring_atom = is_ring
            child.charge = charge
            child.is_aromatic = is_arom
            root.add_child(child)

    # Phase 2: Parse inner spheres (inside parentheses)
    inner_spheres = []
    if pos < n and tokens[pos].type == TOK_OPEN:
        pos += 1  # consume '('
        inner_spheres, pos = _parse_sphere_block(tokens, pos, {TOK_CLOSE})
        if pos < n and tokens[pos].type == TOK_CLOSE:
            pos += 1  # consume ')'

    # Attach inner spheres to the LAST sphere-0 child
    if inner_spheres and root.children:
        last_child = root.children[-1]
        _attach_spheres(last_child, inner_spheres)

    # Phase 3: Parse outer continuation (after ')')
    outer_spheres = []
    if pos < n:
        outer_spheres, pos = _parse_sphere_block(tokens, pos, set())

    # Attach outer spheres to the remaining sphere-0 children (all except last)
    if outer_spheres and len(root.children) > 1:
        remaining = root.children[:-1]
        _distribute_outer_spheres(remaining, outer_spheres)
    elif outer_spheres and len(root.children) == 1:
        # Only one sphere-0 child, attach outer as additional spheres
        _attach_spheres(root.children[0], outer_spheres)

    return root


def _attach_spheres(parent_node, spheres):
    """Attach parsed sphere data as children in the tree.
    spheres[0] = sphere 1 data (children of parent_node)
    spheres[1] = sphere 2 data (grandchildren) etc."""
    if not spheres:
        return

    # Sphere 1: create children of parent_node
    sphere1 = spheres[0]  # list of branches
    for branch in sphere1:
        for sym, bond, is_ring, charge, is_arom in branch:
            child = TreeNode(sym, bond=bond, sphere=parent_node.sphere + 1)
            child.is_ring_atom = is_ring
            child.charge = charge
            child.is_aromatic = is_arom
            parent_node.add_child(child)

    # For sphere 2+, distribute branches among expandable (non-H) children
    if len(spheres) > 1:
        _attach_deeper_spheres(parent_node, spheres[1:])


def _attach_deeper_spheres(parent_node, remaining_spheres):
    """Recursively attach deeper sphere data to expandable descendants."""
    if not remaining_spheres:
        return

    # Get expandable children (non-H atoms can have further neighbors)
    expandable = [c for c in parent_node.children if c.atom != "H"]

    sphere_data = remaining_spheres[0]  # list of branches
    # Each branch corresponds to one expandable child
    for i, branch in enumerate(sphere_data):
        if i >= len(expandable):
            break
        target = expandable[i]
        for sym, bond, is_ring, charge, is_arom in branch:
            child = TreeNode(sym, bond=bond, sphere=target.sphere + 1)
            child.is_ring_atom = is_ring
            child.charge = charge
            child.is_aromatic = is_arom
            target.add_child(child)

    # Continue to next sphere level
    if len(remaining_spheres) > 1:
        # Collect all expandable atoms at the current deepest level
        next_expandable = []
        for c in expandable:
            for gc in c.children:
                if gc.atom != "H":
                    next_expandable.append(gc)

        if next_expandable and remaining_spheres[1]:
            # Create a synthetic parent to pass to _attach_deeper_spheres
            # Actually we just distribute the next sphere's branches
            next_sphere_data = remaining_spheres[1]
            for i, branch in enumerate(next_sphere_data):
                if i >= len(next_expandable):
                    break
                target = next_expandable[i]
                for sym, bond, is_ring, charge, is_arom in branch:
                    child = TreeNode(
                        sym, bond=bond, sphere=target.sphere + 1
                    )
                    child.is_ring_atom = is_ring
                    child.charge = charge
                    child.is_aromatic = is_arom
                    target.add_child(child)

            # Handle spheres beyond 3 with the same pattern
            if len(remaining_spheres) > 2:
                deeper_expandable = []
                for ne in next_expandable:
                    for c in ne.children:
                        if c.atom != "H":
                            deeper_expandable.append(c)
                _distribute_flat(deeper_expandable, remaining_spheres[2:])


def _distribute_flat(expandable_nodes, remaining_spheres):
    """Distribute remaining sphere data across expandable nodes, one level at a time."""
    current_expandable = expandable_nodes
    for sphere_data in remaining_spheres:
        if not current_expandable or not sphere_data:
            break
        next_expandable = []
        for i, branch in enumerate(sphere_data):
            if i >= len(current_expandable):
                break
            target = current_expandable[i]
            for sym, bond, is_ring, charge, is_arom in branch:
                child = TreeNode(sym, bond=bond, sphere=target.sphere + 1)
                child.is_ring_atom = is_ring
                child.charge = charge
                child.is_aromatic = is_arom
                target.add_child(child)
                if child.atom != "H":
                    next_expandable.append(child)
        current_expandable = next_expandable


def _distribute_outer_spheres(remaining_children, spheres):
    """Distribute outer continuation sphere data to the remaining sphere-0 children."""
    if not spheres or not remaining_children:
        return

    # Filter to expandable (non-H)
    expandable = [c for c in remaining_children if c.atom != "H"]

    if not expandable:
        return

    # Sphere 1 of outer continuation = children of these remaining nodes
    sphere1 = spheres[0]
    for i, branch in enumerate(sphere1):
        if i >= len(expandable):
            break
        target = expandable[i]
        for sym, bond, is_ring, charge, is_arom in branch:
            child = TreeNode(sym, bond=bond, sphere=target.sphere + 1)
            child.is_ring_atom = is_ring
            child.charge = charge
            child.is_aromatic = is_arom
            target.add_child(child)

    # Deeper spheres
    if len(spheres) > 1:
        next_expandable = []
        for e in expandable:
            for c in e.children:
                if c.atom != "H":
                    next_expandable.append(c)
        _distribute_flat(next_expandable, spheres[1:])


def _format_atom_smiles(node):
    """Format a single atom for SMILES output."""
    atom = node.atom
    charge = node.charge

    # Aromatic atoms use lowercase
    if node.is_aromatic and atom in ("C", "N", "O", "S", "P"):
        atom = atom.lower()

    # Determine if we need bracket notation
    needs_bracket = (
        charge != 0
        or atom not in ORGANIC_SUBSET
        and atom not in ("c", "n", "o", "s", "p")
    )

    if needs_bracket:
        s = "[" + atom
        if charge > 0:
            s += "+"
        elif charge < 0:
            s += "-"
        s += "]"
        return s
    return atom


def _format_bond_smiles(bond):
    """Format a bond for SMILES output."""
    if bond == "=":
        return "="
    if bond == "#":
        return "#"
    # Aromatic and single bonds are implicit
    return ""


def _find_ring_pairs(root):
    """Find ring closure pairs from nodes marked with _ring_close_to.

    During tree construction, sphere-0 ring chains set _ring_close_to
    on the last chain atom to point back to root.  We collect these
    into (ancestor, closer, ring_id) triples for SMILES ring digits.
    """
    pairs = []
    ring_id = 1
    _collect_ring_close(root, pairs, ring_id)
    return pairs


def _collect_ring_close(node, pairs, ring_id):
    """Walk the tree collecting _ring_close_to markers."""
    if node._ring_close_to is not None:
        pairs.append((node._ring_close_to, node, len(pairs) + 1))
    for child in node.children:
        _collect_ring_close(child, pairs, ring_id)


def tree_to_smiles(root, ring_pairs=None):
    """Convert the atom tree to a SMILES string via DFS."""
    if ring_pairs is None:
        ring_pairs = []

    # Build lookup: which nodes open/close ring digits
    ring_opens = {}   # node -> list of ring_ids that OPEN here
    ring_closers = {} # node -> (ring_id, bond) that CLOSES here

    for ancestor, closer, rid in ring_pairs:
        ring_opens.setdefault(ancestor, []).append(rid)
        ring_closers[closer] = (rid, closer.bond)

    return _gen_smiles(root, ring_opens, ring_closers)


def _is_implicit_h(node):
    """Check if a hydrogen child can be suppressed (implicit in SMILES).
    Keep H explicit only if it has children or a charge."""
    if node.atom != "H":
        return False
    if node.charge != 0:
        return False
    if node.children:
        return False
    return True


def _gen_smiles(node, ring_opens, ring_closers):
    """Recursive SMILES generation with implicit hydrogen suppression."""
    parts = []

    # Write atom
    parts.append(_format_atom_smiles(node))

    # If this node is a ring closer, emit the ring digit and stop
    if node in ring_closers:
        rid, bond = ring_closers[node]
        parts.append(str(rid) if rid < 10 else "%" + str(rid))
        return "".join(parts)

    # Write ring-open digits at this node
    if node in ring_opens:
        for rid in ring_opens[node]:
            parts.append(str(rid) if rid < 10 else "%" + str(rid))

    # Collect children, suppressing implicit H
    real_children = [c for c in node.children if not _is_implicit_h(c)]

    if not real_children:
        return "".join(parts)

    if len(real_children) == 1:
        child = real_children[0]
        parts.append(_format_bond_smiles(child.bond))
        parts.append(_gen_smiles(child, ring_opens, ring_closers))
    else:
        # Pick the "heaviest" child as main chain (most descendants)
        main = max(real_children, key=lambda c: _subtree_size(c))
        branches = [c for c in real_children if c is not main]

        for bc in branches:
            parts.append("(")
            parts.append(_format_bond_smiles(bc.bond))
            parts.append(_gen_smiles(bc, ring_opens, ring_closers))
            parts.append(")")

        parts.append(_format_bond_smiles(main.bond))
        parts.append(_gen_smiles(main, ring_opens, ring_closers))

    return "".join(parts)


def _subtree_size(node):
    """Count total nodes in subtree."""
    count = 1
    for c in node.children:
        count += _subtree_size(c)
    return count


def hose_to_smiles(hose_code, central_atom="C"):
    """Convert a HOSE code to a SMILES fragment string.

    Args:
        hose_code: The HOSE code string (e.g., "C(=CC/HC,HHH/HHC)")
        central_atom: The central atom symbol (from the Nucleus field, e.g., "C")

    Returns:
        A SMILES fragment string, or None if conversion fails.
    """
    if not hose_code or not hose_code.strip():
        return None

    try:
        root = parse_hose(hose_code, central_atom)
        ring_pairs = _find_ring_pairs(root)
        smiles = tree_to_smiles(root, ring_pairs)
        return smiles if smiles else None
    except Exception:
        return None


def extract_central_atom(nucleus_str):
    """Extract the central atom symbol from a nucleus string like 'C-4', 'H-1', 'N-15'."""
    if not nucleus_str:
        return "C"
    # Take the first character(s) before '-' or digit
    parts = nucleus_str.split("-")
    if parts:
        atom = parts[0].strip()
        if atom:
            return atom
    return "C"


if __name__ == "__main__":
    # Quick test
    test_cases = [
        ("HHHC", "C"),
        ("C(=CC/HC,HHH/HHC)", "C"),
        ("H=CC(HC,HHC/H=C,@OCC/", "C"),
        ("=CCC(HC,HHH,HHH/HHC/", "C"),
        ("HHHC(HCC/H=C,HHH/", "C"),
        ("CC(HHC,H=C/=CC,CC/HC,HHH,HHH,HHH)H=C/HC/", "C"),
    ]

    print("HOSE -> SMILES Converter Test")
    print("=" * 70)
    for hose, central in test_cases:
        smiles = hose_to_smiles(hose, central)
        print(f"  HOSE:   {hose}")
        print(f"  SMILES: {smiles}")
        print()
