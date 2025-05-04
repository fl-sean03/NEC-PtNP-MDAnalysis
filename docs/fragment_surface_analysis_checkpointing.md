# `fragment_surface_analysis.py` Frame-Level Checkpointing Implementation Plan

This document outlines the implementation plan for adding frame-level checkpointing capabilities to the `fragment_surface_analysis.py` script. This will allow the script to resume processing a molecular dynamics trajectory from the last successfully processed frame in case of interruptions or crashes, without requiring manual intervention to specify a starting frame.

## Implementation Details

The `fragment_surface_analysis.py` script will be modified to support resuming processing from a partial output file.

1.  **Output File Path Argument:** Ensure the script receives the full path to its primary output CSV file (`surface_analysis.csv`) as a required argument. This path will be provided by the `pipeline.py` script based on the configuration in `config.py`.
2.  **Determine Start Frame:** Upon script startup, the script will determine the frame index from which to begin processing the trajectory:
    *   It will check if the specified output CSV file (`surface_analysis.csv`) exists and is not empty.
    *   If the file exists and is not empty, the script will read the last line of the CSV file. It will parse this line to extract the `frame_index` of the last successfully processed frame that was written to the file. The script will then set its starting frame for trajectory processing to `last_frame_index + 1`.
    *   If the file does not exist or is empty, the script will set its starting frame to `0` (the first frame of the trajectory).
3.  **File Handling:**
    *   Open the output CSV file in append mode (`'a'`) if resuming from an existing file (i.e., if the determined start frame is greater than 0).
    *   Open the output CSV file in write mode (`'w'`) if starting a new file (i.e., if the determined start frame is 0).
    *   Implement logic to write the CSV header row only when opening the file in write mode (i.e., when starting a new file). This prevents duplicate headers when resuming.
4.  **Frame Processing and Flushing:**
    *   Iterate through the trajectory frames starting from the determined start frame up to the end of the trajectory.
    *   For each frame:
        *   Load the necessary atomic coordinates and other data for the current frame into memory.
        *   Perform the surface analysis calculations for that frame (e.g., calculating distances to surface Pt atoms).
        *   Format the results for the current frame as a row (or multiple rows, depending on the output structure) for the output CSV file.
        *   Write the formatted data row(s) to the output CSV file.
        *   **Immediately flush the write buffer to disk** using the file object's `flush()` method. This is a critical step to ensure that the data for each processed frame is saved to disk as soon as it's available, minimizing the amount of data lost if the script crashes between processing frames.
5.  **Completion:** Once the script has successfully processed all frames from the starting frame to the end of the trajectory, it will exit with a zero exit code, signaling successful completion to the calling pipeline script.

This implementation ensures that `fragment_surface_analysis.py` can pick up where it left off, significantly reducing the time and resources required to restart the analysis after an interruption, especially for long simulations. By processing and flushing frame by frame, it also keeps memory usage manageable.