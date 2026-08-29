"""
IUPAC Graph Engine
==================
Module 1 of CS-Fundamentals-Lab.

Represents a saturated, acyclic hydrocarbon (alkane) as a graph:
    - nodes  = carbon atoms
    - edges  = single bonds (C-C skeleton; hydrogens are implicit)

Pipeline:
    1. Build the molecule graph from an edge list.
    2. DFS from every leaf (degree-1 atom) to every other leaf to find
       ALL longest paths (candidate "main chains") -- a tree can have
       several paths tied for longest.
    3. Apply IUPAC tie-breaking rules to pick the real main chain:
         a. most carbons
         b. most substituents attached to it
         c. lowest locants for substituents
    4. Number the chosen chain in the direction that gives the lowest
       locants to substituents.
    5. Identify each substituent's structure (methyl, ethyl, propyl,
       isopropyl, butyl, isobutyl, sec-butyl, tert-butyl) via DFS on
       the branch, and name it.
    6. Assemble the final name: locants + substituent names
       (alphabetical, with di/tri/tetra... multipliers) + parent chain
       name + "ane".

Scope / limitations (documented, not hidden):
    - Acyclic alkanes only (CnH2n+2). No rings, no double/triple bonds,
      no heteroatoms or functional groups.
    - Branch identification covers substituents up to 4 carbons
      (methyl -> tert-butyl). Larger substituents fall back to a
      generic "Cn-yl" label.
"""

from collections import defaultdict, deque
from itertools import combinations

# ---------------------------------------------------------------------------
# Naming tables
# ---------------------------------------------------------------------------

CHAIN_PREFIXES = {
    1: "Meth", 2: "Eth", 3: "Prop", 4: "But", 5: "Pent",
    6: "Hex", 7: "Hept", 8: "Oct", 9: "Non", 10: "Dec",
    11: "Undec", 12: "Dodec",
}

MULTIPLIER_PREFIXES = {
    1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta",
    6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca",
}


class Molecule:
    """Undirected graph of carbon atoms connected by single bonds."""

    def __init__(self, num_atoms, edges):
        self.num_atoms = num_atoms
        self.adj = defaultdict(list)
        for a, b in edges:
            self.adj[a].append(b)
            self.adj[b].append(a)
        self._validate_is_tree()

    def _validate_is_tree(self):
        # An alkane carbon skeleton must be a tree: exactly num_atoms - 1
        # edges, and connected (this engine doesn't support rings yet).
        edge_count = sum(len(v) for v in self.adj.values()) // 2
        if edge_count != self.num_atoms - 1:
            raise ValueError(
                "Graph is not a tree (ring or disconnected component "
                "detected). This engine only supports acyclic alkanes."
            )

    def leaves(self):
        return [n for n in range(self.num_atoms) if len(self.adj[n]) == 1]

    def degree(self, node):
        return len(self.adj[node])

    # -- DFS: every simple path from `start` to `end` in a tree is unique,
    #    so a straightforward DFS walk gives us the path directly. --------
    def dfs_path(self, start, end):
        visited = {start}
        path = [start]

        def _walk(node):
            if node == end:
                return True
            for neighbor in self.adj[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    path.append(neighbor)
                    if _walk(neighbor):
                        return True
                    path.pop()
            return False

        _walk(start)
        return path


# ---------------------------------------------------------------------------
# Step 2-3: find candidate main chains via DFS between every pair of leaves
# ---------------------------------------------------------------------------

def find_longest_chains(mol: Molecule):
    """DFS between every pair of leaf atoms; return all paths tied for
    the maximum length (a list of atom-index lists)."""
    leaves = mol.leaves()
    if len(leaves) < 2:
        # Single atom (methane) -- no leaf pairs to walk.
        return [[0]] if mol.num_atoms == 1 else []

    best_len = -1
    best_paths = []
    for a, b in combinations(leaves, 2):
        path = mol.dfs_path(a, b)
        if len(path) > best_len:
            best_len = len(path)
            best_paths = [path]
        elif len(path) == best_len:
            best_paths.append(path)
    return best_paths


def substituent_count(mol: Molecule, chain):
    """How many atoms hang off this chain (i.e. branch points)."""
    chain_set = set(chain)
    count = 0
    for atom in chain:
        for neighbor in mol.adj[atom]:
            if neighbor not in chain_set:
                count += 1
    return count


def choose_main_chain(mol: Molecule, candidates):
    """Apply IUPAC tie-break rules: longest, then most substituents."""
    if len(candidates) == 1:
        return candidates[0]
    max_subs = max(substituent_count(mol, c) for c in candidates)
    tied = [c for c in candidates if substituent_count(mol, c) == max_subs]
    # Still tied after this -> any is chemically equivalent for naming
    # purposes here; lowest-locant numbering (step 4) resolves the rest.
    return tied[0]


# ---------------------------------------------------------------------------
# Step 5: identify + name a substituent branch hanging off the main chain
# ---------------------------------------------------------------------------

def name_substituent(mol: Molecule, root, forbidden):
    """
    root: the first atom of the branch (attached to the main chain)
    forbidden: set of atoms NOT to walk into (the main chain itself)
    Returns (name, size_in_carbons).
    """
    # Collect every atom in this branch via DFS, since it's a tree.
    visited = {root}
    stack = [root]
    branch_atoms = [root]
    while stack:
        node = stack.pop()
        for neighbor in mol.adj[node]:
            if neighbor not in forbidden and neighbor not in visited:
                visited.add(neighbor)
                branch_atoms.append(neighbor)
                stack.append(neighbor)

    size = len(branch_atoms)

    if size == 1:
        return "methyl", size
    if size == 2:
        return "ethyl", size
    if size == 3:
        # propyl (straight) vs isopropyl (branched at attachment point)
        if mol.degree(root) - 1 == 2:  # two branches off the root itself
            return "isopropyl", size
        return "propyl", size
    if size == 4:
        # Classify by shape: straight butyl, isobutyl, sec-butyl, tert-butyl
        sub_edges = [(a, b) for a in branch_atoms for b in mol.adj[a]
                     if b in branch_atoms and a < b]
        sub_mol = Molecule(4, _relabel(branch_atoms, sub_edges))
        degrees = sorted(len(sub_mol.adj[n]) for n in range(4))
        root_relabel = branch_atoms.index(root)
        root_degree = len(sub_mol.adj[root_relabel])

        if degrees == [1, 1, 1, 3] and root_degree == 3:
            return "tert-butyl", size
        if degrees == [1, 1, 1, 3] and root_degree == 1:
            return "isobutyl", size
        if degrees == [1, 1, 2, 2] and root_degree == 2:
            return "sec-butyl", size
        return "butyl", size

    # Larger substituents: generic fallback name (documented limitation)
    return f"({size}-carbon substituent)", size


def _relabel(atoms, edges):
    index = {atom: i for i, atom in enumerate(atoms)}
    return [(index[a], index[b]) for a, b in edges]


# ---------------------------------------------------------------------------
# Step 4 + 6: number the chain both directions, pick lowest locants, name it
# ---------------------------------------------------------------------------

def name_molecule(mol: Molecule):
    if mol.num_atoms == 1:
        return "Methane"

    candidates = find_longest_chains(mol)
    main_chain = choose_main_chain(mol, candidates)
    chain_set = set(main_chain)

    # Gather (position_along_chain_from_start, atom) for every substituent
    def substituents_for_direction(chain):
        subs = []
        for pos, atom in enumerate(chain, start=1):
            for neighbor in mol.adj[atom]:
                if neighbor not in chain_set:
                    name, _ = name_substituent(mol, neighbor, chain_set | {atom} - {atom} | chain_set)
                    subs.append((pos, name))
        return subs

    forward = substituents_for_direction(main_chain)
    backward = substituents_for_direction(list(reversed(main_chain)))

    def locant_signature(subs):
        return sorted(pos for pos, _ in subs)

    # Lowest-locants rule: compare first point of difference
    chosen = forward if locant_signature(forward) <= locant_signature(backward) else backward

    # Group identical substituents together, sort alphabetically
    grouped = defaultdict(list)
    for pos, name in chosen:
        grouped[name].append(pos)

    name_parts = []
    for sub_name in sorted(grouped.keys(), key=lambda s: s.replace("tert-", "").replace("sec-", "")):
        positions = sorted(grouped[sub_name])
        count = len(positions)
        multiplier = MULTIPLIER_PREFIXES.get(count, f"{count}-")
        locants = ",".join(str(p) for p in positions)
        name_parts.append(f"{locants}-{multiplier}{sub_name}")

    chain_len = len(main_chain)
    parent = CHAIN_PREFIXES.get(chain_len, f"C{chain_len}-") + "ane"

    if name_parts:
        return "-".join(name_parts) + parent.lower()[0] + parent[1:]  # keep parent capitalized correctly
    return parent


# ---------------------------------------------------------------------------
# Demo / self-test with known compounds
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tests = []

    # Methane: C
    tests.append(("Methane", Molecule(1, [])))

    # Octane: straight chain of 8
    tests.append(("Octane", Molecule(8, [(i, i + 1) for i in range(7)])))

    # Isobutane (2-methylpropane): central C bonded to 3 CH3
    tests.append(("2-methylpropane (isobutane)", Molecule(4, [(0, 1), (0, 2), (0, 3)])))

    # 2,3-dimethylbutane: main chain 0-1-2-3, methyls on atoms 1 and 2
    tests.append((
        "2,3-dimethylbutane",
        Molecule(6, [(0, 1), (1, 2), (2, 3), (1, 4), (2, 5)])
    ))

    # 3-ethylpentane: main chain 0-1-2-3-4, ethyl (atoms 5-6) on atom 2
    tests.append((
        "3-ethylpentane",
        Molecule(7, [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5), (5, 6)])
    ))

    print(f"{'Expected':30s} | Engine output")
    print("-" * 60)
    for expected, mol in tests:
        result = name_molecule(mol)
        print(f"{expected:30s} | {result}")
