---
layout: default
title: element_pct.py
script: element_pct.py
description: 'Prints the atomic percentage of one element in an imago.skl file by counting its atoms in the fractional coordinate block. A small helper used by the DEAD-MD genetic algorithm and handy on its own.'
parent: Script Catalog
---

# element_pct

Prints the atomic percentage of a given element in an Imago skeleton
file. The script parses the `fractional` block of the skeleton, counts
how many atom lines begin with the requested element symbol, and
divides by the total atom count declared on the `fractional` line.
It is a small helper that the DEAD-MD genetic algorithm
(`deadmd.py`) calls to score candidate structures, and it is useful
on its own for a quick composition check.

## Usage

**Commandline:**

```
element_pct.py FILENAME ELEMENT
```

Both arguments are required and positional.

**Variable Table:**

| Variable   | Accepted Values               | Description |
|:----------:|:-----------------------------:|:------------|
|FILENAME    |String; file path              | Imago skeleton file. |
|ELEMENT     |Chemical symbol, e.g. h, c, si | Element to count. The comparison is case-insensitive. |

## Output

A single number on standard output, the percentage to two decimal
places, for example `33.33`. On a missing file, or a skeleton whose
`fractional` line does not declare an atom count, it prints an error
and exits with status 1.

## Importable Usage

```python
from element_pct import calculate_element_percentage
pct = calculate_element_percentage("imago.skl", "si")
```

The function raises `ValueError` if the atom count cannot be
determined and `FileNotFoundError` if the file does not exist.
