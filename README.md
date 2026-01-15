# Tree-point-cloud-processing-and-analysis

# Tree Point Cloud Analysis Tool

A Python tool for automated extraction of forestry metrics from 3D Point Clouds. This script processes batches of tree point clouds (e.g., PLY, PCD files) to calculate **Diameter at Breast Height (DBH)** using multiple fitting algorithms, **Tree Height**, **Crown Diameter**, and **Trunk Lean Angle**.

## 🚀 Features

* **Batch Processing:** Automatically iterates through species subfolders and processes all point cloud files.
* **Trunk Isolation:** Uses DBSCAN clustering to robustly separate the main trunk from noise and lower branches.
* **Multi-Method DBH Calculation:**
    * **Least Squares Circle:** Best fit circle.
    * **Minimum Enclosing Circle (MEC):** Welzl's algorithm.
    * **Cylindrical Fit:** Internal and External cylinder fitting.
    * **Elliptical Fit:** Major/minor axis and equivalent diameter.
* **Crown & Height:** Calculates full tree height and crown diameter.
* **Output:** Generates individual CSV reports per species and a global summary CSV for all trees.

## 🛠️ Installation

### Prerequisites
* Python 3.8+
* Anaconda (recommended)

### Dependencies
Install the required libraries using pip:

```bash
pip install numpy scipy open3d
