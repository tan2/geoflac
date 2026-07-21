# GeoFLAC User Guide: Post-Processing & Visualization with VisIt

This guide explains how to convert GeoFLAC's raw binary outputs into VTK formats and visualize the results using LLNL's **VisIt** (Virtual Instrument 3D, developed by Lawrence Livermore National Laboratory).

---

## 1. Converting Solver Binary Outputs to VTK

GeoFLAC outputs simulation data in raw binary formats (`*.0`, `*.rs`, etc.). To visualize these in VisIt, you must first convert them to XML-based VTK formats using the Python utilities located in the `util/` directory:

| Utility Script | Input Data | VTK Output Format | Description |
|:---------------|:-----------|:------------------|:------------|
| **`flac2vtk.py`** | `vel.0`, `temp.0`, `phase.0`, etc. | `.vts` (Structured Grid) | Converts cell and node-wise variables (e.g., velocity, temperature, viscosity, phases) onto the grid. |
| **`flacmarker2vtk.py`** | `_markers.0` or `marker1.*` | `.vtp` (PolyData) | Converts Lagrangian marker locations and properties (including thermochronology ages) into point datasets. |
| **`tcages2grid.py`** | `flacmarker.*.vtp` | `.vts` (Structured Grid) | Interpolates marker-based thermochronology ages onto the element centers of the grid for contoured grid rendering. |

### Running the Conversion:
Navigate to your simulation output directory (where the binary files reside) and run:
```bash
# 1. Convert standard grid variables (Required)
python3 ../../util/flac2vtk.py .

# 2. Convert Lagrangian marker points (Optional - only if markers are used)
#    (add -t option to compute offline thermochronology ages)
python3 ../../util/flacmarker2vtk.py -t .

# 3. Interpolate thermochronology ages to grid elements (Optional - only if using thermochronology)
python3 ../../util/tcages2grid.py .
```
This generates:
*   `flac.*.vts` (Grid files for each output step)
*   `flacmarker.*.vtp` (Marker point files, if Step 2 was run)

---

## 2. Setting Up VisIt

VisIt is a free, open-source, interactive parallel visualization and graphical analysis tool. 

1.  **Download and Install**: Download VisIt from the [LLNL VisIt Website](https://uci.llnl.gov/visit).
2.  **Launch VisIt**: Start VisIt via your terminal or GUI application.
3.  **Active Source Window**: You will see the main control window (GUI) on the left and the visualization window on the right.

---

## 3. Visualizing Grid Data (`.vts`)

### Step 1: Open the Grid Database
1.  Click **Open** in the Sources panel of the Main Window.
2.  Navigate to your output directory and select **`flac.*.vts database`**. VisIt automatically groups individual time step files into a single database time-series.

### Step 2: Render the Grid Mesh
1.  Click **Add** -> **Mesh** -> **mesh**.
2.  Click **Draw**. The undeformed or deformed mesh grid will appear in the visualization window.
3.  To style the mesh lines, double-click the **Mesh** plot entry to open the Mesh Plot Attributes (adjust line width, color, or opacity).

### Step 3: Render Temperature or Viscosity (Pseudocolor Plot)
1.  Click **Add** -> **Pseudocolor** -> **Temperature** (or `Viscosity`, `Pressure`, `Plastic strain`).
2.  Click **Draw**.
3.  Double-click the **Pseudocolor** plot entry to change the colormap (e.g., `hot` for temperature, `viridis` or `spectral` for strain).
4.  **Steady Colorbar Limits**: By default, VisIt dynamically scales the min/max color limits for each timestep based on the current frame's data range, which causes colors to shift erratically during animations. To lock the colorbar limits, double-click the **Pseudocolor** entry to open the Pseudocolor Plot Attributes. In the *Limits* section:
    *   Check **Min** and **Max**.
    *   Enter fixed constant values (e.g., Min=`0.0`, Max=`1200.0` for temperature; Min=`-15.0`, Max=`-12.5` for `Strain rate`; or Min=`0.0`, Max=`2.0` for `Plastic strain`).

### Step 4: Render Velocity Vector Fields
1.  Click **Add** -> **Vector** -> **Velocity**.
2.  Click **Draw**. Velocity vectors will display as arrows at node points.
3.  **Adjusting Arrow Density**: Double-click **Vector** to open the Vector Plot Attributes. In the *Sampling* tab, change the stride count (e.g., skip 2 or 3 nodes) to keep the viewport clean.
4.  **Steady Vector Scaling**: By default, vector arrows are auto-scaled, which causes arrow lengths to fluctuate confusingly between timesteps. To keep scaling consistent:
    *   Double-click the **Vector** entry to open the Vector Plot Attributes.
    *   Under the *Vector Geometry* tab, set the **Scale** options to a fixed constant value (rather than automatic scaling). This ensures velocity magnitudes can be visually compared accurately over time.

---

## 4. Visualizing Lagrangian Marker Data (`.vtp`)

Markers represent discrete rock packages moving through the mesh.

### Step 1: Open the Marker Database
1.  Click **Open** in the Sources panel and select **`flacmarker.*.vtp database`**.

### Step 2: Display Markers as Points
1.  Click **Add** -> **Pseudocolor** -> **Phase** (or `age_ZFT`, `age_AFT` to display thermochronology).
2.  Click **Draw**.
3.  Double-click **Pseudocolor** in the plot list. In the *Formatting* tab:
    *   Change **Point Type** to *Sphere* or *Box* for 3D/2D representation.
    *   Set **Point Size** (in pixels or absolute units) so that individual markers are visible without overlapping too heavily.

### Step 3: Synchronizing Grid and Marker Databases (Time Correlation)
When you load both the grid database (`.vts`) and the marker database (`.vtp`) in the same VisIt session:
1. VisIt automatically creates a **time correlation** between the two databases.
2. This links their cycles and timesteps together under a single active time slider. When you step through the animation using the VCR controls, both datasets will advance in sync, allowing you to see the Lagrangian markers moving exactly with the deforming grid.
3. If they do not synchronize automatically, you can manually establish a correlation by going to **Controls** -> **Time Relation** in the menu bar, selecting the two databases, and creating a new correlation group.

---

## 5. Advanced VisIt Operators & Filters

Operators are applied to plots to cut, filter, or modify data:

### A. Thresholding by Phase or Strain (The Threshold Operator)
To isolate specific rock phases (e.g., showing only the subducting slab) or to show only high-shear zones:
1.  In the Plot List, click on your plot (e.g., Pseudocolor of `Plastic strain`).
2.  Click **Operators** -> **Selection** -> **Threshold**.
3.  Expand the plot name in the list, double-click **Threshold** to open attributes:
    *   Select the variable to threshold by (e.g., `Phase` or `Plastic strain`).
    *   Set the **Lower bound** and **Upper bound** (e.g., set `Plastic strain` bounds to `0.5` to `2.5` to display only active shear zones).
4.  Click **Apply** and **Draw**.

### B. Making 2D Slices (The Slice Operator)
If you have 3D models or want to cut a plane:
1.  Select your plot, click **Operators** -> **Slicing** -> **Slice**.
2.  In Slice Attributes, configure the cutting plane normal (e.g., Z-intercept or arbitrary plane orientation).

### C. Elevating Topography (The Elevate Operator)
To render 1D topography profiles as 2D curves or to deform a surface plot vertically by elevation:
1.  Select your plot (e.g., Pseudocolor of `topo`).
2.  Click **Operators** -> **Transforms** -> **Elevate**.
3.  Double-click **Elevate** and choose `topo` as the elevating variable.

---

## 6. Animating and Exporting

1.  **Time Slider**: Use the VCR-style controls at the bottom of the Main Window to play through the time sequence.
2.  **Saving Images**: Click **File** -> **Save window** to export the current viewport as a `.png`, `.jpg`, or `.tiff` image.
3.  **Generating Movies**: Click **File** -> **Save movie** and select the format (e.g., MPEG, AVI) to render the entire transient evolutionary sequence.
