# Pipeline-Level Checkpointing Implementation Plan (using done.txt files)

This document outlines the implementation plan for a pipeline-level checkpointing system using "done.txt" files. This system allows the pipeline to resume from the last successfully completed step in case of interruptions or crashes.

## Pipeline Execution Logic (`pipeline.py`)

The `pipeline.py` script will be modified to incorporate checkpointing logic for each analysis step using dedicated "done.txt" files.

**Implementation:**

1.  **Define Done Files:** For each analysis script in the pipeline sequence, define a specific "done file" name (e.g., `simulation_info.done`, `classify_pt_atoms.done`, `fragment_surface_analysis.done`, `residence_event_analysis.done`, `binding_metrics.done`). These files will be placed in a `.done` subdirectory within the main output directory specified in `config.py`.
2.  **Iterate Through Steps:** The `pipeline.py` script will iterate through the defined sequence of analysis steps in the correct order.
3.  **Check for Completion:** For each step, it will construct the full path to the expected "done file" based on the main output directory from `config.py`, the `.done` subdirectory, and the step's defined done file name.
4.  **Conditional Execution:**
    *   `pipeline.py` will check for the existence of the "done file" for the current step.
    *   If the "done file" exists, `pipeline.py` will assume the step was completed successfully in a previous run. It will print a message indicating that the step is being skipped and move directly to the next step in the sequence.
    *   If the "done file" is missing, `pipeline.py` will execute the script for that step using `subprocess`.
5.  **Marking Completion:** After a script for a step successfully completes (i.e., the `subprocess` call returns a zero exit code), `pipeline.py` will create an empty file at the designated "done file" path for that step. This marks the step as completed.

## Individual Analysis Scripts

The individual analysis scripts (e.g., `simulation_info.py`, `classify_pt_atoms.py`, `fragment_surface_analysis.py`, etc.) will run from start to finish whenever executed by the pipeline. They do not need to implement internal checkpointing for this pipeline-level system to work.

## Configuration (`config.py`)

The `fragment_surface_analysis.py` script will be modified to support resuming processing from a partial output file.

**Implementation:**

1.  **Output File Path Argument:** Ensure the script receives the full path to its output CSV file (`surface_analysis.csv`) as a required argument.
2.  **Check for Existing Output:** Upon script startup, check if the specified output CSV file exists and is not empty.
3.  **Determine Start Frame:**
    *   If the file exists and is not empty, read the last line of the CSV file. Parse this line to extract the `frame_index` of the last successfully processed frame. The script will then set its starting frame for trajectory processing to `last_frame_index + 1`.
    *   If the file does not exist or is empty, the script will set its starting frame to `0`.
4.  **File Handling:**
    *   Open the output CSV file in append mode (`'a'`) if resuming from an existing file.
    *   Open the output CSV file in write mode (`'w'`) if starting a new file.
    *   Implement logic to write the CSV header row only when opening in write mode (i.e., starting a new file).
5.  **Frame Processing and Flushing:**
    *   Iterate through the trajectory frames starting from the determined start frame.
    *   After processing each frame and writing the corresponding data row(s) to the output CSV file, immediately flush the write buffer to disk using the file object's `flush()` method. This minimizes data loss in case of a crash between frames.

## Configuration (`config.py`)

The `config.py` file needs to clearly define the output paths for all analysis steps.

**Implementation:**

1.  Ensure that the `ANALYSIS_PARAMETERS` dictionary (or equivalent structure) in `config.py` includes parameters specifying the output directory and filenames for each script.
2.  `pipeline.py` will read these output paths from `config.py` to perform the existence checks for pipeline-level checkpointing.

This comprehensive checkpointing system will significantly improve the robustness of the analysis pipeline against interruptions.