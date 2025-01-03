# Circular Track Reconstruction with ROOT

This project implements a simulation and visualization framework for reconstructing circular tracks in a 2D tracker using ROOT, a framework widely used in high-energy physics for data analysis and visualization. The program generates circular events, simulates tracker hits, and identifies the parameters of the circular trajectory (center and radius) using iterative cross-section analysis.

---

## Features

1. **Circular Event Generation**:
   - Generates synthetic tracker hits that lie on a circular trajectory.
   - Includes optional angular constraints for partial circles.

2. **Hit Visualization**:
   - Displays hits and reconstructed tracks in a 2D plot using ROOT's graphics capabilities.

3. **Track Reconstruction**:
   - Iteratively analyzes cross-sections of 3D cones formed by tracker hits.
   - Refines parameters (center and radius) using a zoom-in approach for increased precision.

---

## Dependencies

This program requires the following:
- **ROOT**: A data analysis framework (download from [ROOT CERN](https://root.cern/)).

---

## File Structure

- **Main Functions**:
  - `int CircTrack()`: Entry point of the program. Simulates an event, generates hits, reconstructs the track, and visualizes results.
  - `void CircularEventGenerator(...)`: Generates hits on a circular trajectory.
  - `track CreateConeCrossections(...)`: Reconstructs circular tracks by analyzing cross-sections.

- **Helper Functions**:
  - `void CreateCircle(...)`: Creates a circle for visualization and analysis.
  - Overloaded `CircularEventGenerator(...)`: Generates partial circular trajectories with angular constraints.

- **Structures**:
  - `struct hit`: Represents a tracker hit with coordinates and radius.
  - `struct track`: Represents the reconstructed circular track.

---
## Usage

### Compilation

Use ROOT's `cling` or your preferred C++ compiler to compile the code:
```bash
root-config --cxx --cflags --libs -o circular_track circular_track.cpp
```

### Execution
Run the compiled binary:
```bash 
root CircTrack.cpp 
```

### Output
The console will display:
```bash
number of hits: 25
Global Maximum Value: 12.3456
Global Maximum Bin: 123
At z value: 4.5678
```
This indicates the bin that indicates the data of the reconstructed circular track. 
A visualization of the circular track will be saved as: circular_track.png

## Parameters
The following parameters influence the simulation and reconstruction:
- n1, n2 (int): Dimensions of the tracker grid.
  Default: n1 = 10, n2 = 10.

- n_points (int): Number of points used for approximating circles.
  Default: 10000.

- delta (double): Bin size for histograms.
  Default: 0.5.

- zoom_scale (double): Zoom-in factor between iterations.
  Default: 10.

- num_of_iterations (int): Number of refinement iterations.
  Default: 4.

- x_min, x_max (double): Initial x-range of the tracker space.
  Default: x_min = -10, x_max = 20.

- y_min, y_max (double): Initial y-range of the tracker space.
  Default: y_min = -10, y_max = 20.

- z_min, z_max (double): Range of z-values for cross-section analysis.
  Default: z_min = 1, z_max = 10.

## Visualization
The generated circular_track.png contains:

- Small ellipses: Representing tracker hits.
- Large ellipse: Representing the reconstructed circular track.
