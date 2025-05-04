# `classify_pt_atoms.py` Documentation

This script classifies Platinum (Pt) atoms within a molecular dynamics simulation based on their coordination numbers. It can either compute the coordination number cutoff automatically by analyzing the first minimum after the main peak in the Pt-Pt Radial Distribution Function (RDF), or use a user-specified cutoff value. The script then uses a KDTree to efficiently count neighbors within the determined cutoff distance and assigns coordination numbers to each Pt atom. Finally, atoms are classified into categories such as "bulk", "facet", "edge", and "vertex" based on their coordination numbers.

## Usage

```bash
python classify_pt_atoms.py \\
    --psf /path/to/system.psf \\
    --dcd /path/to/trajectory.dcd \\
    --output-dir /path/to/output_dir \\
    [--prefix system_prefix] \\
    [--frame-index INDEX] \\
    [--use-cutoff CUT_OFF] \\
    [--r-max 10.0] [--nbins 200] \\
    [--min-frame INDEX] [--max-frame INDEX] \\
    [--plot-rdf]
```

## Arguments

*   `--psf` (required): Path to the PSF file.
*   `--dcd` (required): Path to the DCD trajectory file.
*   `--output-dir` (required): Directory to save output files.
*   `--prefix` (optional, default: `system`): Prefix for output file names.
*   `--frame-index` (optional): Index of a single frame to analyze (default: middle frame).
*   `--min-frame` (optional): Start frame for RDF averaging.
*   `--max-frame` (optional): End frame for RDF averaging.
*   `--use-cutoff` (optional): Bypass RDF calculation and use this cutoff distance (in Å) directly.
*   `--r-max` (optional, default: `10.0`): Maximum distance (in Å) for RDF calculation.
*   `--nbins` (optional, default: `200`): Number of bins for RDF calculation.
*   `--plot-rdf` (optional, flag): Save a plot of the RDF curve with the chosen cutoff.

## Outputs

The script generates the following files in the specified output directory:

*   `{prefix}_coordination_numbers.csv`: A CSV file containing the atom index, coordination number, classification, and a boolean indicating if the atom is on the surface.
*   `{prefix}_coordination_numbers.txt`: A text file grouping atom indices by their classification and listing surface atoms.
*   `(optional) {prefix}_rdf_curve.png`: A PNG image visualizing the calculated RDF curve and the determined cutoff distance (generated if `--plot-rdf` is used).

## Classification Logic

Atoms are classified based on their coordination number (`cn`) using the following thresholds:

*   `bulk`: `cn >= 11`
*   `111 facet`: `cn == 9`
*   `100 facet`: `cn == 8`
*   `edge`: `6 <= cn < 8`
*   `vertex`: `cn < 6`