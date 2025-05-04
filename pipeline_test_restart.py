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

    # Define the steps and their corresponding done files
    # Define the steps and their corresponding done files
    steps = [
        {
            "name": "Simulation Info",
            "script": "simulation_info.py",
            "done_file": "simulation_info.done",
        },
        {
            "name": "Pt Classification",
            "script": "classify_pt_atoms.py",
            "done_file": "pt_classification.done",
            "optional_args": [
                ("frame_index", config.pt_classification_frame_index),
                ("min_frame", config.pt_classification_min_frame),
                ("max_frame", config.pt_classification_max_frame),
                ("use_cutoff", config.pt_classification_use_cutoff),
            ],
            "flag_args": [
                ("plot_rdf", config.pt_classification_plot_rdf),
            ]
        },
        {
            "name": "Fragment Surface Analysis",
            "script": "fragment_surface_analysis_test_restart.py", # Use the test script
            "done_file": "fragment_surface_analysis.done",
             "optional_args": [
                ("frame_index", config.fragment_surface_analysis_frame_index),
            ],
            "list_args": [
                ("fragment_ids", config.fragment_surface_analysis_fragment_ids),
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
        },
        {
            "name": "Residence Event Analysis",
            "script": "residence_event_analysis.py",
            "done_file": "residence_event_analysis.done",
        },
        {
            "name": "Binding Metrics",
            "script": "binding_metrics.py",
            "done_file": "binding_metrics.done",
        }
    ]

    # Ensure the base output directory and the .done directory exist
    done_dir = os.path.join(config.base_output_dir, ".done")
    os.makedirs(config.base_output_dir, exist_ok=True)
    os.makedirs(done_dir, exist_ok=True)

    # Run the pipeline steps with checkpointing
    for step in steps:
        step_name = step["name"]
        script_name = step["script"]
        done_file_path = os.path.join(done_dir, step["done_file"])

        print(f"\n--- Running {step_name} ---")

        if os.path.exists(done_file_path):
            print(f"Skipping {step_name}: Done file found at {done_file_path}")
            continue

        # Construct the command dynamically
        command = [
            sys.executable, os.path.join("src", script_name)
        ]

        # Add arguments based on step name
        if step_name == "Simulation Info":
            command.extend([
                "--psf", config.psf_file,
                "--dcd", config.dcd_file,
                "--steps-per-frame", str(config.steps_per_frame),
                "--output-dir", config.sim_info_output_dir
            ])
        elif step_name == "Pt Classification":
             command.extend([
                "--psf", config.psf_file,
                "--dcd", config.dcd_file,
                "--output-dir", config.pt_classification_output_dir,
                "--prefix", config.pt_classification_prefix,
                "--r-max", str(config.pt_classification_r_max),
                "--nbins", str(config.pt_classification_nbins)
            ])
        elif step_name == "Fragment Surface Analysis":
            command.extend([
                "--psf", config.psf_file,
                "--dcd", config.dcd_file,
                "--pt-csv", config.pt_classification_csv,
                "--output-dir", config.surface_analysis_output_dir
            ])
        elif step_name == "Filter Surface Fragments":
            command.extend([
                "--input-csv", config.surface_analysis_csv,
                "--output-dir", config.surface_filtered_output_dir,
                "--cutoff", str(config.filter_surface_fragments_cutoff)
            ])
        elif step_name == "Residence Event Analysis":
            command.extend([
                "--filtered-csv", config.surface_filtered_csv,
                "--sim-info-json", config.sim_info_json,
                "--min-res-time", str(config.residence_event_analysis_min_res_time),
                "--max-off-time", str(config.residence_event_analysis_max_off_time),
                "--output-dir", config.residence_events_output_dir
            ])
        elif step_name == "Binding Metrics":
             command.extend([
                "--residence_events_csv", config.residence_events_csv,
                "--filtered_surface_csv", config.surface_filtered_csv,
                "--simulation_properties_csv", config.sim_info_csv,
                "--output_molecule_csv", config.binding_metrics_molecule_csv,
                "--output_facet_csv", config.binding_metrics_facet_csv
            ])

        # Add optional arguments if they exist
        if "optional_args" in step:
            for arg_name, arg_value in step["optional_args"]:
                if arg_value is not None:
                    command.extend([f"--{arg_name.replace('_', '-')}", str(arg_value)])

        # Add list arguments if they exist
        if "list_args" in step:
             for arg_name, arg_values in step["list_args"]:
                 if arg_values is not None:
                     command.extend([f"--{arg_name.replace('_', '-')}"] + [str(v) for v in arg_values])

        # Add flag arguments if they exist and are True
        if "flag_args" in step:
            for arg_name, arg_value in step["flag_args"]:
                if arg_value:
                    command.append(f"--{arg_name.replace('_', '-')}")

        run_command(command)

        # Create the done file after successful execution
        try:
            with open(done_file_path, 'w') as f:
                pass # Create an empty file
            print(f"Created done file: {done_file_path}")
        except IOError as e:
            print(f"Error creating done file {done_file_path}: {e}")
            pass # Just print a warning for now.


    print("\n--- Pipeline finished successfully ---")

if __name__ == "__main__":
    main()