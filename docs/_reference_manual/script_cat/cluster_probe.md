---
layout: default
title: cluster_probe.py
script: cluster_probe.py
description: 'Writes a starter clusterrc.py cluster-dispatch settings file by asking the SLURM scheduler what it knows: partitions and their node sizes from sinfo, chargeable accounts from sacctmgr. Site conventions it cannot discover are left as marked FILL IN blanks.'
parent: Script Catalog
---

# cluster_probe

Discovers a site's cluster dispatch settings and writes a starter
`clusterrc.py`. That file is pure data: the parameters that the
dispatch generator behind
[build_initial_potentials.py](build_initial_potentials.html) and the
kaleidoscope flight runner reads when it submits work to a scheduler.
This tool creates a first copy of it by asking the scheduler. It runs
`sinfo` for the queues and the size of their nodes, and `sacctmgr`
for the accounts you may charge, and writes a `clusterrc.py` with
what it could learn already filled in and a brief plain-language note
on every setting.

What the scheduler cannot report is convention and policy: how a
worker brings up the Imago environment (`worker_init`) and which
account to charge are site-specific and human-chosen, so those are
left as clearly marked `FILL IN` blanks. Complete them before
dispatching.

**Two honesty rules.**

* **Only the scheduler is trusted, not the login node.** `sinfo`
  reports the *compute* nodes' core and memory counts. The tool does
  not read the login node's own CPU layout, which would describe the
  wrong machine, so no login-node facts ever reach the file.
* **A heterogeneous cluster is not guessed at.** When the cluster's
  nodes differ in core count or memory, the tool does not silently
  pick one value. It leaves that setting blank and lists the values
  it saw in a comment, so you choose deliberately.

The tool is self-contained: it carries its own copy of the settings
schema and only ever *writes* a `clusterrc.py`, never reads one. A
test keeps that schema identical to the one `clusterrc.py` itself
declares, so the two cannot drift apart.

## Usage

**Commandline:**

```
cluster_probe.py [-o PATH] [--force] [--help]
```

Run it on a login node of the cluster you intend to dispatch to.
The dispatch code looks for `clusterrc.py` first in the current
working directory, as a per-run override, and then in the directory
named by `$IMAGO_RC`, so copy the completed file to `$IMAGO_RC` to
make it the site-wide default.

**Option Table:**

| Option         | Defaults                              | Description |
|:--------------:|:-------------------------------------:|:------------|
|-o / --output   |clusterrc.py in the current directory  | Path for the starter file. |
|--force         |Refuse and keep the existing file      | Overwrite an existing output file. |
|--help          |N/A                                    | Print usage information and exit. |

## Output

A `clusterrc.py` whose `parameters_and_defaults()` returns the
dispatch settings: the partitions discovered with their node core
and memory counts, the accounts discovered, a default topology, and
the `worker_init` and account fields marked `FILL IN`. The closing
message names what was discovered and reminds you to complete the
blanks.
