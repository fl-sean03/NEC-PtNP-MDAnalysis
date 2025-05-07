# Molecular Dynamics Analysis of Pt-NEC Interactions

This project provides a pipeline and scripts for analyzing molecular dynamics simulations of systems containing Platinum (Pt) nanoparticles and NEC fragments. The analysis focuses on understanding the interactions between the NEC fragments and the Pt surface, including atom classification, surface proximity analysis, residence event quantification, and binding metrics calculation.

## Motivation

Understanding the binding behavior and interactions of molecules with nanoparticle surfaces is crucial in various fields, including catalysis, materials science, and drug delivery. This project aims to provide a robust and automated workflow to extract key quantitative insights from molecular dynamics simulations of such systems, specifically focusing on Pt-NEC interactions.

## How it Works

The analysis is performed through a series of Python scripts located in the `src/` directory, orchestrated by pipeline scripts in the workspace root. The workflow follows these main steps:

1.  **Simulation Information:** Extracts basic metadata from the simulation trajectory and structure files.
2.  **Pt Atom Classification:** Classifies the Pt atoms on the nanoparticle surface based on their coordination environment (e.g., bulk, facet, edge, vertex). This step can use an automatically determined cutoff from the Radial Distribution Function (RDF) or a user-specified value.
3.  **Fragment Surface Analysis:** Analyzes the proximity of NEC fragments to the classified Pt nanoparticle surface.
4.  **Filter Surface Fragments:** Filters the surface analysis results based on a distance cutoff to identify fragments and atoms that are considered "on the surface".
5.  **Residence Event Analysis:** Quantifies how long NEC fragments reside on the Pt surface based on the filtered surface data and defined time criteria.
6.  **Binding Metrics:** Computes binding statistics such as mean residence time, dissociation constant (K_D), and free energy of binding (ΔG), aggregated per molecule and per Pt facet type.
7.  **Iterative Binding Metrics (in `pipeline_iterative_binding.py`)**: Runs the binding metrics calculation iteratively over increasing time windows to assess convergence.
8.  **Plot Facet Metrics (in `pipeline_iterative_binding.py`)**: Generates plots of the iterative binding metrics.

The pipeline scripts automate the execution of these steps in the correct order, ensuring that the output of one step serves as the input for the next.

## Project Structure

*   `configs/`: Contains configuration files for the analysis pipelines.
*   `data/`: Contains example input simulation files (PSF and DCD).
*   `src/`: Contains the core Python scripts for each analysis step.
*   `output/`: Default directory for saving all generated analysis results (CSV, JSON, Markdown, plots).
*   `docs/`: Contains documentation files (in Markdown) for scripts and the overall pipeline.
*   `sample_scripts/`: Contains additional example scripts, potentially for specific analysis tasks or visualizations.
*   `pipeline.py`: The main script to run the standard analysis workflow.
*   `pipeline_iterative_binding.py`: A pipeline script to run the analysis workflow including iterative binding metrics and plotting.
*   `README.md`: This file.

## Getting Started

1.  Ensure you have the necessary dependencies installed (e.g., MDAnalysis, pandas, numpy, matplotlib, scipy). A `requirements.txt` file is provided for easier environment setup.
2.  Place your simulation PSF and DCD files in the `data/` directory or update the paths in the desired configuration file in the `configs/` directory.
3.  Modify the configuration file in the `configs/` directory to adjust parameters for the analysis steps and specify output locations.
4.  Run the desired pipeline script from the root directory of the project, specifying the configuration file:

    ```bash
    python pipeline.py --config configs/config.py
    # or for the iterative pipeline
    python pipeline_iterative_binding.py --config configs/config_pipeline_iterative.py
    ```

The script will execute the analysis and save the results in the specified output directory.

## Configuration

Configuration files are located in the `configs/` directory. Each configuration file uses a dictionary structure (`ANALYSIS_PARAMETERS`) to group parameters by the script they apply to. You can edit these files in any text editor to customize the input data, output locations, and analysis settings.

## Documentation

Detailed documentation for each individual script and the overall pipeline design can be found in the `docs/` directory.

*   [`docs/scripts/pipeline.md`](docs/scripts/pipeline.md): Documentation for the standard analysis pipeline.
*   [`docs/iterative_binding_metrics_pipeline.md`](docs/iterative_binding_metrics_pipeline.md): Documentation for the iterative binding metrics pipeline.
*   [`docs/scripts/simulation_info.md`](docs/scripts/simulation_info.md): Documentation for the simulation information script.
*   [`docs/scripts/classify_pt_atoms.md`](docs/scripts/classify_pt_atoms.md): Documentation for the Pt atom classification script.
*   [`docs/scripts/fragment_surface_analysis.md`](docs/scripts/fragment_surface_analysis.md): Documentation for the fragment surface analysis script.
*   [`docs/scripts/filter_surface_fragments.md`](docs/scripts/filter_surface_fragments.md): Documentation for the filter surface fragments script.
*   [`docs/scripts/residence_event_analysis.md`](docs/scripts/residence_event_analysis.md): Documentation for the residence event analysis script.
*   [`docs/scripts/binding_metrics.md`](docs/scripts/binding_metrics.md): Documentation for the binding metrics script.
*   [`docs/scripts/simulation_info.md`](docs/scripts/simulation_info.md): Documentation for the simulation information script.