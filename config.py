import os

# --- Analysis Parameters ---
# Configure the parameters for each step of the analysis pipeline here.
# These values will be used by the pipeline.py script.

ANALYSIS_PARAMETERS = {
    "general": {
        # Base directory for input files (e.g., PSF, DCD)
        "base_input_dir": "data",
        # Base directory for all output files
        "base_output_dir": "output/short", # Changed to output to output/short

        # Input file paths relative to base_input_dir
        "psf_file": "short/0HPt.psf", # Pointing to data/short/0HPt.psf
        "dcd_file": "short/Pt+0HNEC.dcd", # Pointing to data/short/Pt+0HNEC.dcd
    },
    "simulation_info": {
        # Number of MD integration steps between saved frames
        "steps_per_frame": 1000,
    },
    "pt_classification": {
        # Prefix for output file names
        "prefix": "0HPt",
        # Single frame to analyze (use None for middle frame)
        "frame_index": None,
        # Start frame for RDF averaging (use None to start from beginning)
        "min_frame": None,
        # End frame for RDF averaging (use None to go to end)
        "max_frame": None,
        # Bypass RDF: use this cutoff (Å) directly (use None to compute from RDF)
        "use_cutoff": 3.5, # Set manual cutoff to 3.5 Å as requested
        # Max distance (Å) for RDF calculation
        "r_max": 10.0,
        # Number of bins for RDF calculation
        "nbins": 200,
        # Save RDF curve plot with cutoff (True/False)
        "plot_rdf": False,
    },
    "fragment_surface_analysis": {
        # IDs of fragments to analyze (list of ints, or None for all)
        "fragment_ids": None, # e.g., [0, 2, 5]
        # Frame index to analyze (single int, or None for all frames)
        "frame_index": None, # e.g., 0
        # Enable verbose output to terminal (True/False)
        "verbose": False,
        # Save data for all NEC atoms; default is only the closest per fragment (True/False)
        "all_atoms": False,
    },
    "filter_surface_fragments": {
        # Distance cutoff (Å) to select fragment-frame pairs and atoms
        "cutoff": 3.5, # Increased cutoff to 3.5 Å
    },
    "residence_event_analysis": {
        # Minimum contiguous residence time (ns) to count as an event
        "min_res_time": 0.001, # Reduced minimum residence time
        # Maximum off-surface gap (ns) to merge adjacent events
        "max_off_time": 1.0, # Increased maximum off-surface time
    },
    "binding_metrics": {
        # Binding metrics parameters are derived from inputs; no specific params here yet.
        # Output paths are defined in the AnalysisConfig class.
    }
}


class AnalysisConfig:
    def __init__(self, parameters=ANALYSIS_PARAMETERS):
        self.params = parameters

        self.base_input_dir = self.params["general"]["base_input_dir"]
        self.base_output_dir = self.params["general"]["base_output_dir"]

        # Ensure output directories exist
        self.sim_info_output_dir = os.path.join(self.base_output_dir, "simulation_info")
        self.pt_classification_output_dir = os.path.join(self.base_output_dir, "pt-classification")
        self.surface_analysis_output_dir = os.path.join(self.base_output_dir, "surface_analysis")
        self.surface_filtered_output_dir = os.path.join(self.base_output_dir, "surface_filtered")
        self.residence_events_output_dir = os.path.join(self.base_output_dir, "residence_events")
        self.binding_metrics_output_dir = os.path.join(self.base_output_dir, "binding_metrics")

        os.makedirs(self.sim_info_output_dir, exist_ok=True)
        os.makedirs(self.pt_classification_output_dir, exist_ok=True)
        os.makedirs(self.surface_analysis_output_dir, exist_ok=True)
        os.makedirs(self.surface_filtered_output_dir, exist_ok=True)
        os.makedirs(self.residence_events_output_dir, exist_ok=True)
        os.makedirs(self.binding_metrics_output_dir, exist_ok=True)


        # --- Input Files ---
        self.psf_file = os.path.join(self.base_input_dir, self.params["general"]["psf_file"])
        self.dcd_file = os.path.join(self.base_input_dir, self.params["general"]["dcd_file"])

        # --- Simulation Info Parameters ---
        self.steps_per_frame = self.params["simulation_info"]["steps_per_frame"]

        # --- Pt Classification Parameters ---
        self.pt_classification_prefix = self.params["pt_classification"]["prefix"]
        self.pt_classification_frame_index = self.params["pt_classification"]["frame_index"]
        self.pt_classification_min_frame = self.params["pt_classification"]["min_frame"]
        self.pt_classification_max_frame = self.params["pt_classification"]["max_frame"]
        self.pt_classification_use_cutoff = self.params["pt_classification"]["use_cutoff"]
        self.pt_classification_r_max = self.params["pt_classification"]["r_max"]
        self.pt_classification_nbins = self.params["pt_classification"]["nbins"]
        self.pt_classification_plot_rdf = self.params["pt_classification"]["plot_rdf"]

        # --- Fragment Surface Analysis Parameters ---
        self.fragment_surface_analysis_fragment_ids = self.params["fragment_surface_analysis"]["fragment_ids"]
        self.fragment_surface_analysis_frame_index = self.params["fragment_surface_analysis"]["frame_index"]
        self.fragment_surface_analysis_verbose = self.params["fragment_surface_analysis"]["verbose"]
        self.fragment_surface_analysis_all_atoms = self.params["fragment_surface_analysis"]["all_atoms"]

        # --- Filter Surface Fragments Parameters ---
        self.filter_surface_fragments_cutoff = self.params["filter_surface_fragments"]["cutoff"]

        # --- Residence Event Analysis Parameters ---
        self.residence_event_analysis_min_res_time = self.params["residence_event_analysis"]["min_res_time"]
        self.residence_event_analysis_max_off_time = self.params["residence_event_analysis"]["max_off_time"]

        # --- Binding Metrics Parameters ---
        # No specific parameters defined in the dictionary yet, but can be added here if needed.

        # --- Intermediate File Paths (derived from outputs of previous steps) ---
        self.pt_classification_csv = os.path.join(
            self.pt_classification_output_dir,
            f"{self.pt_classification_prefix}_coordination_numbers.csv"
        )
        self.surface_analysis_csv = os.path.join(
            self.surface_analysis_output_dir,
            "surface_analysis.csv"
        )
        self.surface_filtered_csv = os.path.join(
            self.surface_filtered_output_dir,
            "surface_filtered.csv"
        )
        self.sim_info_json = os.path.join(
            self.sim_info_output_dir,
            "simulation_info.json"
        )
        self.sim_info_csv = os.path.join(
            self.sim_info_output_dir,
            "simulation_info.csv"
        )
        self.residence_events_csv = os.path.join(
            self.residence_events_output_dir,
            "residence_events.csv"
        )
        self.binding_metrics_molecule_csv = os.path.join(
            self.binding_metrics_output_dir,
            "molecule_results.csv"
        )
        self.binding_metrics_facet_csv = os.path.join(
            self.binding_metrics_output_dir,
            "facet_results.csv"
        )

# Example usage:
# config = AnalysisConfig()
# print(config.psf_file)
# print(config.filter_surface_fragments_cutoff)