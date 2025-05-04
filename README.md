# Molecular Dynamics Analysis of Pt-NEC Interactions

This project provides a pipeline and scripts for analyzing molecular dynamics simulations of systems containing Platinum (Pt) nanoparticles and NEC fragments. The analysis focuses on understanding the interactions between the NEC fragments and the Pt surface, including atom classification, surface proximity analysis, residence event quantification, and binding metrics calculation.

## Motivation

Understanding the binding behavior and interactions of molecules with nanoparticle surfaces is crucial in various fields, including catalysis, materials science, and drug delivery. This project aims to provide a robust and automated workflow to extract key quantitative insights from molecular dynamics simulations of such systems, specifically focusing on Pt-NEC interactions.

## How it Works

The analysis is performed through a series of Python scripts located in the `src/` directory, orchestrated by a main `pipeline.py` script. The workflow follows these main steps:

1.  **Simulation Information:** Extracts basic metadata from the simulation trajectory and structure files.
2.  **Pt Atom Classification:** Classifies the Pt atoms on the nanoparticle surface based on their coordination environment (e.g., bulk, facet, edge, vertex). This step can use an automatically determined cutoff from the Radial Distribution Function (RDF) or a user-specified value.
3.  **Fragment Surface Analysis:** Analyzes the proximity of NEC fragments to the classified Pt nanoparticle surface.
4.  **Filter Surface Fragments:** Filters the surface analysis results based on a distance cutoff to identify fragments and atoms that are considered "on the surface".
5.  **Residence Event Analysis:** Quantifies how long NEC fragments reside on the Pt surface based on the filtered surface data and defined time criteria.
6.  **Binding Metrics:** Computes binding statistics such as mean residence time, dissociation constant (K_D), and free energy of binding (ΔG), aggregated per molecule and per Pt facet type.

The `pipeline.py` script automates the execution of these steps in the correct order, ensuring that the output of one step serves as the input for the next.

## Project Structure

*   `data/`: Contains example input simulation files (PSF and DCD).
*   `src/`: Contains the core Python scripts for each analysis step.
*   `output/`: Default directory for saving all generated analysis results (CSV, JSON, Markdown, plots).
*   `docs/scripts/`: Contains documentation files (in Markdown) for each script and the overall pipeline.
*   `sample_scripts/`: Contains additional example scripts, potentially for specific analysis tasks or visualizations.
*   `config.py`: Configuration file to set input/output paths and parameters for the pipeline.
*   `pipeline.py`: The main script to run the entire analysis workflow.
*   `README.md`: This file.

## Getting Started

1.  Ensure you have the necessary dependencies installed (e.g., MDAnalysis, pandas, numpy, matplotlib, scipy). A `requirements.txt` file could be added for easier environment setup.
2.  Place your simulation PSF and DCD files in the `data/` directory or update the paths in `config.py`.
3.  Modify the `config.py` file to adjust parameters for the analysis steps and specify output locations.
4.  Run the pipeline from the root directory of the project:

    ```bash
    python pipeline.py
    ```

The script will execute the analysis and save the results in the specified output directory.

## Configuration

The `config.py` file is the central place to configure the analysis pipeline. It uses a dictionary structure (`ANALYSIS_PARAMETERS`) to group parameters by the script they apply to. You can edit this file in any text editor to customize the input data, output locations, and analysis settings.

## Documentation

Detailed documentation for each individual script and the overall pipeline design can be found in the `docs/scripts/` directory.

*   [`pipeline.md`](docs/scripts/pipeline.md): Documentation for the overall analysis pipeline.
*   [`simulation_info.md`](docs/scripts/simulation_info.md): Documentation for the simulation information script.
*   [`classify_pt_atoms.md`](docs/scripts/classify_pt_atoms.md): Documentation for the Pt atom classification script.
*   [`fragment_surface_analysis.md`](docs/scripts/fragment_surface_analysis.md): Documentation for the fragment surface analysis script.
*   [`filter_surface_fragments.md`](docs/scripts/filter_surface_fragments.md): Documentation for the filter surface fragments script.
*   [`residence_event_analysis.md`](docs/scripts/residence_event_analysis.md): Documentation for the residence event analysis script.
*   [`binding_metrics.md`](docs/scripts/binding_metrics.md): Documentation for the binding metrics script.