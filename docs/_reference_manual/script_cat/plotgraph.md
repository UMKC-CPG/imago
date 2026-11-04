---
layout: default
title: plotgraph.py
script: plotgraph.py
description: 'General-purpose plotter for Imago output. Reads a column-oriented data file, writes a plotting script for matplotlib, plotly or veusz, and runs it. Has dedicated figure types for density of states, optical properties and symmetric band structures, with control over subplots, curve stacking, axis limits and curve styling.'
parent: Script Catalog
---

# plotgraph

The default plotting tool for Imago results. It reads a
column-oriented data file (a `.plot` file from an Imago calculation or
from a post-processing script), **writes a plotting script** for the
chosen display library, and then **runs that script**. Because the
generated script is an ordinary Python or veusz file left in the
working directory, you can open it, edit it, and re-run it by hand
when a figure needs a tweak the command line does not offer.

It supports four figure types and should be your first choice for
them: `general`, density of states (`dos`), optical properties
(`optc`), and symmetric band structure (`sybd`). If your data does
not fit those categories, or the script fails on it, look at the
other post-processing scripts in the catalog.

**Defaults.** Every option's default value is read from the resource
control file `plotgraphrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run. The curve style,
colour and marker palettes are read from `curve_styles.dat`,
`curve_colors.dat` and `curve_marks.dat` in `$IMAGO_DATA`.

## Usage

**Commandline:**

```
plotgraph.py [-d {mpl,veusz,plotly}] [-vl] [-i INFILE] [-o OUTFILE]
             [-t TITLE] [-pc] [-pl]
             [-po {portrait,landscape}] [-ff {letter,A4}]
             [-fh HEIGHT] [-fw WIDTH] [-ft {general,dos,optc,sybd}]
             [-l] [-mx] [-xc XCOL] [-yc YCOLS]
             [-xi XMIN] [-xa XMAX] [-yi YMIN] [-ya YMAX]
             [-xm XMULT] [-ym YMULT] [-xai XINC] [-xmi XMINOR]
             [-yma YMAJOR] [-ymi YMINOR] [-lx] [-ly]
             [-ss SUBSEP] [-sp N [N ...]] [-cs N [N ...]]
             [-cp SEP [SEP ...]]
             [-ct W [W ...]] [-cta START] [-cte STEP] [-cti SIZE]
             [-cy S [S ...]] [-cya START] [-cye STEP] [-cyi SIZE]
             [-cc C [C ...]] [-cca START] [-cce STEP] [-cci SIZE]
             [-cm M [M ...]] [-cma START] [-cme STEP] [-cmi SIZE]
             [--help]
```

**Examples:**

```
plotgraph.py -i gs_dos-fb.t.plot -ft dos            # total DOS
plotgraph.py -i PDOS.plot -ft dos -cs 1 -sp 4       # 4 stacked PDOS panels
plotgraph.py -i gs_sybd-fb.plot -ft sybd            # band structure
plotgraph.py -i gs_optc-eb.t.eps2.plot -ft optc -xa 40   # eps2 to 40 eV
plotgraph.py -i data.plot -yc "2 4..6" -cs 4 -pl    # columns 2,4,5,6 in one panel
```

**Variable Table:**

| Variable   | Accepted Values                          | Description |
|:----------:|:----------------------------------------:|:------------|
|INFILE      |String; file path                         | Column-oriented data file. The first line holds the column headers, except for `sybd` files where it is data. |
|OUTFILE     |String; file path                         | Name of the generated plotting script. |
|TITLE       |String                                    | Title for the figure set. |
|HEIGHT, WIDTH |float; percent of the page             | Figure height and width as a percentage of the page size. |
|XCOL        |integer; 1-based                          | Column holding the x values when every curve shares one x column. |
|YCOLS       |quoted string of integers and ranges      | Columns to plot as y values, e.g. `"2 4..6"` means columns 2, 4, 5 and 6. |
|XMIN, XMAX, YMIN, YMAX |float                          | Explicit axis limits. When not given the axis autoscales. |
|XMULT, YMULT |float                                    | Multiple to which autoscaled axis limits are rounded. |
|XINC        |float                                     | Spacing between x-axis major ticks. |
|XMINOR, YMINOR |integer                                | Number of minor tick intervals per major interval. |
|YMAJOR      |integer                                   | Number of y-axis major ticks. |
|SUBSEP      |integer                                   | Vertical separation between subplots. |
|N           |integer; may be a list                    | Count used by `-sp` and `-cs`. See the figure layout section. |
|SEP         |float; may be a list                      | Vertical offset between successive curves stacked in one subplot. |
|W, S, C, M  |explicit lists                            | Curve widths, style tokens, colour names and marker tokens. |
|START, STEP, SIZE |integer (float for width)           | Generate a list of widths, styles, colours or marks from the palette files instead of listing them explicitly. |

**Option Table:**

| Option                    | Defaults          | Description |
|:-------------------------:|:-----------------:|:------------|
|-d / --display             |mpl                | Display library to use: `mpl` (matplotlib), `veusz` or `plotly`. |
|-vl / --veuszlaunch        |Off                | Also open the generated `.vsz` document in the veusz GUI. Only meaningful with `-d veusz`. |
|-i / --infile              |gs_dos-fb.t.plot   | Input file name. |
|-o / --outfile             |Xscript.py         | Output script name. A leading `X` is replaced by the display name, so the default becomes `mplscript.py`, `plotlyscript.py` or `veuszscript.py`. |
|-t / --title               |Plotted by `$USER` | Title for the figure set. |
|-pc / --printcmd           |Off                | Print the command line used onto the figure. |
|-pl / --printlgd           |Off                | Print a legend on each subplot. |
|-po / --pageo              |portrait           | Page orientation: `portrait` or `landscape`. |
|-ff / --figformat          |letter             | Page size: `letter` or `A4`. |
|-fh / --figh               |40                 | Figure height, percent of page. |
|-fw / --figw               |80                 | Figure width, percent of page. |
|-ft / --figt               |general            | Figure type: `general`, `dos`, `optc` or `sybd`. See below. |
|-l / --link_subplots       |On                 | Link the subplots of a figure so they share axes. |
|-mx / --multix             |Off                | The input has a separate x column immediately to the left of every y column. Cannot be combined with `-xc`. |
|-xc / --xcol               |1                  | The single x column shared by all curves. |
|-yc / --ycol               |All columns except the x column | Which columns to plot, as a quoted list with `..` ranges. |
|-xi / -xa / -yi / -ya      |Autoscale          | Explicit x-min, x-max, y-min, y-max. |
|-xm / --xmult              |5                  | Round autoscaled x limits to a multiple of this. |
|-ym / --ymult              |1                  | Round autoscaled y limits to a multiple of this. |
|-xai / --xaxisinc          |5                  | X-axis major tick increment. |
|-xmi / --xminor            |4                  | X-axis minor ticks per major interval. |
|-yma / --ymajor            |10                 | Number of y-axis major ticks. |
|-ymi / --yminor            |4                  | Y-axis minor ticks per major interval. |
|-lx / --linkx              |Off                | Link the x axes of all subplots. |
|-ly / --linky              |Off                | Link the y axes of all subplots. |
|-ss / --subplotsep         |10                 | Separation between subplots. |
|-sp / --subplotsperpage    |1                  | Subplots per figure. A list cycles. |
|-cs / --curves             |1                  | Curves per subplot. A list cycles. |
|-cp / --curvesep           |10                 | Vertical offset between stacked curves in a subplot. A list cycles. |
|-ct / --curvewidth         |Generated: start 0.5, step 0.25, size 1 | Explicit list of curve widths, or generate one with `-cta`, `-cte`, `-cti`. |
|-cy / --curvestyle         |Generated: start 0, step 1, size 1 | Explicit list of style tokens, or generate one with `-cya`, `-cye`, `-cyi`. Note that `-cy` cannot take a list of tuple styles on the command line; use the generated form. |
|-cc / --curvecolor         |Generated: start 18, step 1, size 1 | Explicit list of colour names, or generate one with `-cca`, `-cce`, `-cci`. Names are xkcd colour names; quote any that contain spaces. |
|-cm / --curvemark          |Generated: start 0, step 1, size 1 | Explicit list of marker tokens, or generate one with `-cma`, `-cme`, `-cmi`. |
|--help                     |N/A                | Print usage information and exit. |

## Option Details

### Figure Types

> -ft / --figt

| Type      | Behaviour |
|:---------:|:----------|
|`general`  | Plain x-y curves. Column headers come from the first line of the file. |
|`dos`      | Density of states. Accepted so that the intent is recorded in the generated script, but at present it is drawn exactly like `general`. |
|`optc`     | Optical properties. As `dos`: currently identical to `general`. |
|`sybd`     | Symmetric band structure. The input file has no header line: column one is the k-path coordinate, the band energies follow, and the last two columns carry the high-symmetry k-point positions and names. All bands are drawn in one subplot, vertical lines are drawn at the high-symmetry points, and the x axis is relabelled with their names. Only the matplotlib backend supports this type; `plotly` and `veusz` refuse it with a message. |

### Choosing Columns

> -xc, -yc, -mx

Columns are numbered from 1. In the usual layout one column holds the
x values (`-xc`, column 1 by default) and every other column is a
curve. `-yc` selects a subset as a quoted, space-separated list in
which `a..b` expands to every integer from a through b, so
`-yc "2 4..6"` plots columns 2, 4, 5 and 6.

Some files instead carry a separate x column immediately to the left
of every y column. For those give `-mx`; the y columns then default to
the even-numbered columns and each is paired with the column before
it. `-mx` and `-xc` describe contradictory layouts and the script
refuses to run if both are given.

### Figure Layout

> -sp, -cs, -cp, -ss

The requested y columns are walked in order and dealt out into
curves, subplots and figures:

* `-cs` says how many curves go into a subplot before a new subplot
  is started.
* `-sp` says how many subplots go into a figure before a new figure
  is started.
* `-cp` is the vertical offset added to each successive curve within
  a subplot, so that several curves can be stacked for comparison
  rather than overlaid.

Each of these accepts a list. The list is consumed one entry per
subplot (or figure) and wraps around when exhausted, so
`-cs 1 3` alternates one-curve and three-curve subplots. A single
value applies to every subplot. For `sybd` the curves-per-subplot
value is overridden so that every band lands in one subplot.

### Curve Appearance

> -ct, -cy, -cc, -cm and their start/step/size companions

Each visual aspect of a curve, its width, line style, colour and
marker, can be set in one of two ways. Give an explicit list
(`-cc "xkcd:red" "xkcd:blue"`), or give a *start*, *step* and *size*
and let the script build the list by stepping through the palette
file in `$IMAGO_DATA`. The generated list is appended to the explicit
one, and both wrap around the palette, so a short list still covers
any number of curves. The three palette files each carry one token
per line for every backend, and only the column for the active
backend is used.

### Axis Limits

> -xi, -xa, -yi, -ya, -xm, -ym

Any limit not given on the command line is determined from the data
and rounded outward to a multiple of `-xm` (x) or `-ym` (y). Give an
explicit limit to override that for one axis end only; the other end
still autoscales.

## Output

The script writes `mplscript.py`, `plotlyscript.py` or a veusz
document, depending on the backend, and then executes it. For
matplotlib the figure window opens interactively with hover labels
on each curve. For veusz the `.vsz` document is written and, with
`-vl`, opened in the veusz GUI.
