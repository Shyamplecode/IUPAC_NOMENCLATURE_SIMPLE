# IUPAC Graph Engine

A small script that looks at a molecule (just the carbon skeleton) and figures out its chemical name — the same way you'd do it by hand in chemistry class, but with code.

## Why I built this

Organic chemistry naming rules always felt like a graph problem to me. You're given a bunch of carbons connected to each other, and you need to:

1. Find the longest chain running through them
2. Notice anything branching off that chain
3. Number things so the branches get the smallest possible numbers
4. Stitch it all into a name

That's basically "find the longest path in a tree, then walk it and label what you find" — which is a nice, contained way to actually use DFS on something other than a textbook example.

## What it does

You give it a list of atoms and which ones are bonded to each other. It gives you back a name like `2,3-dimethylbutane`.

It only handles simple alkanes for now — plain carbon chains with single bonds, no rings, no double bonds, nothing exotic. That was a deliberate choice: get the core logic (longest chain + branch naming + numbering) working and correct before adding more chemistry on top.

## How it works, roughly

- **Build the graph.** Atoms are nodes, bonds are edges.
- **Find the main chain.** Every atom with only one connection is an "end" of the molecule. The script runs a DFS between every pair of ends and keeps whichever path turns out longest. If there's a tie, it picks whichever chain has more stuff branching off it — that's the actual tie-breaking rule chemists use.
- **Number it both ways.** Once it knows the main chain, it counts along it from each end and checks which direction gives the branches lower numbers. Lower numbers win.
- **Name the branches.** Each branch gets its own little DFS to figure out its shape — a single carbon is "methyl", a branch that itself forks is "isopropyl" or "tert-butyl", and so on.
- **Put it together.** Branches get sorted alphabetically, duplicate branches get grouped with "di-", "tri-", etc., and the whole thing gets glued onto the base name (`butane`, `pentane`, `hexane`...).

## Running it

```bash
python3 iupac_engine.py
```

This runs a handful of built-in test molecules — methane, octane, isobutane, 2,3-dimethylbutane, and 3-ethylpentane — and prints what the engine names each one, next to what it should be. Right now all five check out.

## Limitations

Being upfront about what this doesn't do yet:

- No rings (cyclohexane and friends are out)
- No double or triple bonds
- No oxygen, nitrogen, or anything besides carbon
- Branches bigger than 4 carbons just get labeled generically instead of properly named

