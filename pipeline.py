import subprocess
import os
import sys
import argparse

# Assuming config.py is in the same directory
from config import AnalysisConfig

def run_command(command, cwd=None):
    """Runs a shell command and checks for errors."""
    print(f"Executing command: {' '.join(command)}")
    try:
        result = subprocess.run(command, check=True, cwd=cwd, capture_output=True, text=True)
        print("Command output:")
        print(result.stdout)
        if result.stderr:
            print("Command error output:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        print(f"Stderr: {e.stderr}")
        sys.exit(f"Pipeline failed at command: {' '.join(command)}")

def main():
    parser = argparse.ArgumentParser(description="Run the molecular dynamics analysis pipeline.")
    # For simplicity, initially hardcode using the default config.py
    # parser.add_argument("--config", default="config.py", help="Path to the configuration file.")
    args = parser.parse_args()

    # Load configuration
    config = AnalysisConfig()
    print("Loaded configuration.")

    # --- Step 1: Simulation Info ---
    print("\n--- Running Simulation Info ---")
    sim_info_cmd = [
        sys.executable, os.path.join("src", "simulation_info.py"),
        "--psf", config.psf_file,
        "--dcd", config.dcd_file,
        "--steps-per-frame", str(config.steps_per_frame),
        "--output-dir", config.sim_info_output_dir
    ]
    run_command(sim_info_cmd)

    # --- Step 2: Pt Classification ---
    print("\n--- Running Pt Classification ---")
    pt_classify_cmd = [
        sys.executable, os.path.join("src", "classify_pt_atoms.py"),
        "--psf", config.psf_file,
        "--dcd", config.dcd_file,
        "--output-dir", config.pt_classification_output_dir,
        "--prefix", config.pt_classification_prefix,
        "--r-max", str(config.pt_classification_r_max),
        "--nbins", str(config.pt_classification_nbins)
    ]
    if config.pt_classification_frame_index is not None:
        pt_classify_cmd.extend(["--frame-index", str(config.pt_classification_frame_index)])
    if config.pt_classification_min_frame is not None:
        pt_classify_cmd.extend(["--min-frame", str(config.pt_classification_min_frame)])
    if config.pt_classification_max_frame is not None:
        pt_classify_cmd.extend(["--max-frame", str(config.pt_classification_max_frame)])
    if config.pt_classification_use_cutoff is not None:
        pt_classify_cmd.extend(["--use-cutoff", str(config.pt_classification_use_cutoff)])
    if config.pt_classification_plot_rdf:
        pt_classify_cmd.append("--plot-rdf")

    run_command(pt_classify_cmd)

    # --- Step 3: Fragment Surface Analysis ---
    print("\n--- Running Fragment Surface Analysis ---")
    frag_surf_cmd = [
        sys.executable, os.path.join("src", "fragment_surface_analysis.py"),
        "--psf", config.psf_file,
        "--dcd", config.dcd_file,
        "--pt-csv", config.pt_classification_csv,
        "--output-dir", config.surface_analysis_output_dir
    ]
    if config.fragment_surface_analysis_fragment_ids is not None:
        frag_surf_cmd.extend(["--fragment-ids"] + [str(f) for f in config.fragment_surface_analysis_fragment_ids])
    if config.fragment_surface_analysis_frame_index is not None:
        frag_surf_cmd.extend(["--frame-index", str(config.fragment_surface_analysis_frame_index)])
    if config.fragment_surface_analysis_verbose:
        frag_surf_cmd.append("--verbose")
    if config.fragment_surface_analysis_all_atoms:
        frag_surf_cmd.append("--all-atoms")
    run_command(frag_surf_cmd)

    # --- Step 4: Filter Surface Fragments ---
    print("\n--- Running Filter Surface Fragments ---")
    filter_surf_cmd = [
        sys.executable, os.path.join("src", "filter_surface_fragments.py"),
        "--input-csv", config.surface_analysis_csv,
        "--output-dir", config.surface_filtered_output_dir,
        "--cutoff", str(config.filter_surface_fragments_cutoff)
    ]
    run_command(filter_surf_cmd)

    # --- Step 5: Residence Event Analysis ---
    print("\n--- Running Residence Event Analysis ---")
    res_event_cmd = [
        sys.executable, os.path.join("src", "residence_event_analysis.py"),
        "--filtered-csv", config.surface_filtered_csv,
        "--sim-info-json", config.sim_info_json,
        "--min-res-time", str(config.residence_event_analysis_min_res_time),
        "--max-off-time", str(config.residence_event_analysis_max_off_time),
        "--output-dir", config.residence_events_output_dir
    ]
    run_command(res_event_cmd)

    # --- Step 6: Binding Metrics ---
    print("\n--- Running Binding Metrics ---")
    binding_metrics_cmd = [
        sys.executable, os.path.join("src", "binding_metrics.py"),
        "--residence_events_csv", config.residence_events_csv,
        "--filtered_surface_csv", config.surface_filtered_csv, # Changed to use filtered surface CSV
        "--simulation_properties_csv", config.sim_info_csv,
        "--output_molecule_csv", config.binding_metrics_molecule_csv,
        "--output_facet_csv", config.binding_metrics_facet_csv
    ]
    run_command(binding_metrics_cmd)

    print("\n--- Pipeline finished successfully ---")

if __name__ == "__main__":
    main()