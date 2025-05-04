import os
import MDAnalysis as mda
import plotly.graph_objects as go

# 1. Define the base directory and select a system
data_dir = "/home/sf2/LabWork/Workspace/9-NAMDAnalysis/data/time-data-10ns"
system_name = "4HPt"  # Change as needed

# 2. Build file paths for the PSF and DCD files
psf_file = os.path.join(data_dir, system_name, "cbs.psf")
dcd_file = os.path.join(data_dir, system_name, f"{system_name}.dcd")

# 3. Load the simulation Universe using MDAnalysis
u = mda.Universe(psf_file, dcd_file)

# 4. Set the trajectory to the first frame
u.trajectory[0]

# 5. Select all NEC atoms using a selection that excludes Pt atoms
nec_atoms = u.select_atoms("not name PT*")

# 6. Obtain connected components (fragments) based on bond connectivity.
#    Since explicit molecule-level info may not be provided, we use the `fragments` attribute.
nec_fragments = nec_atoms.fragments

# 7. Print some information about the fragments
print("Number of NEC fragments (by connectivity):", len(nec_fragments))
for i, frag in enumerate(nec_fragments):
    print(f"NEC fragment {i+1} has {len(frag)} atoms")

# 8. Choose one fragment to visualize (e.g., the first one)
fragment = nec_fragments[0]
positions = fragment.positions

# 9. Create an interactive 3D Plotly scatter plot for the selected fragment
trace = go.Scatter3d(
    x=positions[:, 0],
    y=positions[:, 1],
    z=positions[:, 2],
    mode='markers',
    marker=dict(
        size=4,   # Marker size for atoms in the fragment
        color='blue'
    ),
    name="NEC Fragment 1"
)

fig = go.Figure(data=[trace])
fig.update_layout(
    title="Interactive 3D Plot of NEC Fragment 1",
    scene=dict(
        xaxis_title="X (Å)",
        yaxis_title="Y (Å)",
        zaxis_title="Z (Å)"
    )
)

fig.show()

