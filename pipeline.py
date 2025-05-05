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
    parser.add_argument("--config", default="config.py", help="Path to the configuration file.")
    args = parser.parse_args()

    # Load configuration
    # Dynamically import the configuration module
    config_path = args.config
    config_dir = os.path.dirname(config_path)
    config_module_name = os.path.basename(config_path).replace('.py', '')

    # Add the config directory to sys.path to allow importing
    sys.path.insert(0, config_dir)
    try:
        config_module = __import__(config_module_name)
        config = config_module.AnalysisConfig()
    finally:
        # Remove the config directory from sys.path
        sys.path.pop(0)
    print("Loaded configuration.")

    # Define the steps and their corresponding done files
    steps = [
        {
            "name": "Simulation Info",
            "script": "simulation_info.py",
            "done_file": "simulation_info.done",
            "command": [
                sys.executable, os.path.join("src", "simulation_info.py"),
                "--psf", config.psf_file,
                "--dcd", config.dcd_file,
                "--steps-per-frame", str(config.steps_per_frame),
                "--output-dir", config.sim_info_output_dir
            ],
           "optional_args": [
               ("start-frame", config.simulation_info_start_frame),
               ("end-frame", config.simulation_info_end_frame),
               ("timestep-fs", config.simulation_info_timestep_fs), # Add timestep-fs argument
           ]
        },
        {
            "name": "Pt Classification",
            "script": "classify_pt_atoms.py",
            "done_file": "pt_classification.done",
            "command": [
                sys.executable, os.path.join("src", "classify_pt_atoms.py"),
                "--psf", config.psf_file,
                "--dcd", config.dcd_file,
                "--output-dir", config.pt_classification_output_dir,
                "--prefix", config.pt_classification_prefix,
                "--r-max", str(config.pt_classification_r_max),
                "--nbins", str(config.pt_classification_nbins)
            ],
            "optional_args": [
                # Use the frame_index specified in the pt_classification config, defaulting to the middle frame if None
                ("frame-index", config.pt_classification_frame_index),
                ("min-frame", config.pt_classification_min_frame),
                ("max-frame", config.pt_classification_max_frame),
                ("use-cutoff", config.pt_classification_use_cutoff),
            ],
            "flag_args": [
                ("plot_rdf", config.pt_classification_plot_rdf),
            ]
        },
        {
            "name": "Fragment Surface Analysis",
            "script": "fragment_surface_analysis.py",
            "done_file": "fragment_surface_analysis.done",
            "command": [
                sys.executable, os.path.join("src", "fragment_surface_analysis.py"),
                "--psf", config.psf_file,
                "--dcd", config.dcd_file,
                "--pt-csv", config.pt_classification_csv,
                "--output-dir", config.surface_analysis_output_dir,
                "--resume"
            ],
            "optional_args": [
               ("frame_index", config.fragment_surface_analysis_frame_index),
               ("start-frame", config.fragment_surface_analysis_start_frame),
               ("end-frame", config.fragment_surface_analysis_end_frame),
            ],
            "list_args": [
               ("fragment-ids", config.fragment_surface_analysis_fragment_ids),
            ],
            "flag_args": [
               ("verbose", config.fragment_surface_analysis_verbose),
               ("all_atoms", config.fragment_surface_analysis_all_atoms),
            ]
        },
        {
            "name": "Filter Surface Fragments",
            "script": "filter_surface_fragments.py",
            "done_file": "filter_surface_fragments.done",
            "command": [
                sys.executable, os.path.join("src", "filter_surface_fragments.py"),
                "--input-csv", config.surface_analysis_csv,
                "--output-dir", config.surface_filtered_output_dir,
                "--cutoff", str(config.filter_surface_fragments_cutoff)
            ]
        },
        {
            "name": "Residence Event Analysis",
            "script": "residence_event_analysis.py",
            "done_file": "residence_event_analysis.done",
            "command": [
                sys.executable, os.path.join("src", "residence_event_analysis.py"),
                "--filtered-csv", config.surface_filtered_csv,
                "--sim-info-json", config.sim_info_json,
                "--min-res-time", str(config.residence_event_analysis_min_res_time),
                "--max-off-time", str(config.residence_event_analysis_max_off_time),
                "--output-dir", config.residence_events_output_dir
            ]
        },
        {
            "name": "Binding Metrics",
            "script": "binding_metrics.py",
            "done_file": "binding_metrics.done",
            "command": [
                sys.executable, os.path.join("src", "binding_metrics.py"),
                "--residence_events_csv", config.residence_events_csv,
                "--filtered_surface_csv", config.surface_filtered_csv, # Changed to use filtered surface CSV
                "--simulation_properties_csv", config.sim_info_csv,
                "--output_molecule_csv", config.binding_metrics_molecule_csv,
                "--output_facet_csv", config.binding_metrics_facet_csv,
                "--min-res-time", str(config.residence_event_analysis_min_res_time)
            ],
            "optional_args": [
                ("analysis-start-frame", config.binding_metrics_analysis_start_frame),
                ("analysis-end-frame", config.binding_metrics_analysis_end_frame),
            ]
        }
    ]

    # Ensure the base output directory and the .done directory exist
    done_dir = os.path.join(config.base_output_dir, ".done")
    os.makedirs(config.base_output_dir, exist_ok=True)
    os.makedirs(done_dir, exist_ok=True)

    # Run the pipeline steps with checkpointing
    for step in steps:
        step_name = step["name"]
        done_file_path = os.path.join(done_dir, step["done_file"])
        command = step["command"]

        print(f"\n--- Running {step_name} ---")

        if os.path.exists(done_file_path):
            print(f"Skipping {step_name}: Done file found at {done_file_path}")
            continue

        # Add optional arguments if they exist
        if "optional_args" in step:
            for arg_name, arg_value in step["optional_args"]:
                if arg_value is not None:
                    command.extend([f"--{arg_name}", str(arg_value)])

        # Add list arguments if they exist
        if "list_args" in step:
             for arg_name, arg_values in step["list_args"]:
                 if arg_values is not None:
                     command.extend([f"--{arg_name}"] + [str(v) for v in arg_values])

        # Add flag arguments if they exist and are True
        if "flag_args" in step:
            for arg_name, arg_value in step["flag_args"]:
                if arg_value:
                    command.append(f"--{arg_name}")

        run_command(command)

        # Create the done file after successful execution
        try:
            with open(done_file_path, 'w') as f:
                pass # Create an empty file
            print(f"Created done file: {done_file_path}")
        except IOError as e:
            print(f"Error creating done file {done_file_path}: {e}")
            # Decide how to handle this error - maybe exit or just warn?
            # For now, we'll just print a warning.
            pass


    print("\n--- Pipeline finished successfully ---")

if __name__ == "__main__":
    main()