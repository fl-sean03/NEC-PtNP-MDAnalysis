# `simulation_info.py` Documentation

This script extracts metadata from a molecular dynamics simulation in PSF/DCD format. It calculates and reports information such as the number of total atoms, platinum (Pt) atoms, NEC fragments, number of frames, integration timestep, time between frames, total simulation time, box dimensions, and atom type counts.

The script outputs the extracted information into three formats: CSV, JSON, and Markdown.

## Usage

```bash
python simulation_info.py \\
    --psf /path/to/system.psf \\
    --dcd /path/to/trajectory.dcd \\
    --steps-per-frame 1000 \\
    --output-dir /path/to/output_dir
```

## Arguments

*   `--psf` (required): Path to the PSF file.
*   `--dcd` (required): Path to the DCD trajectory file.
*   `--steps-per-frame` (required): Number of MD integration steps between saved frames (e.g., 1000).
*   `--output-dir` (optional, default: `.`): Directory to save the output files.

## Outputs

The script generates the following files in the specified output directory:

*   `simulation_info.csv`: A CSV file containing a table of properties and their extracted values.
*   `simulation_info.json`: A JSON file with the structured metadata.
*   `simulation_info.md`: A human-readable Markdown summary of the simulation metadata.

## Details

*   Integration timestep is reported in femtoseconds (fs).
*   All other time metrics (time between frames, start time, end time, total time) are converted to nanoseconds (ns).
*   Box dimensions are reported in Ångströms (Å).