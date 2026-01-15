# Tree-point-cloud-trunk-segmentation-and-measurement

# Tree Point Cloud Trunk Analysis Tool

A Python tool for automated extraction of forestry metrics from 3D Point Clouds. This script processes batches of tree point clouds (e.g., PLY, PCD files) to calculate **Diameter at Breast Height (DBH)** using multiple fitting algorithms, **Tree Height**, **Crown Diameter**, and **Trunk Lean Angle**.

<img width="2087" height="996" alt="image" src="https://github.com/user-attachments/assets/319e8871-a55a-4308-b679-cd7ee938d5c8" />


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

## Configurations

* **Input Folder:** All the point cloud files should be saved in a folder, and dividing species folders is the recommended structure:
* /Input_Folder
    ├── Acer_platanoides/
    │   ├── tree_01.ply
    │   ├── tree_02.pcd
    ├── Betula_pendula/
    │   ├── tree_01.ply
    │   └── ...
* **Parameters:**
     * **Trunk Extraction:** This is the first height to cluster trunk out of other leaves, branches... The default setting is 0 to 1.5 meters.
     * **Trunk Slicing:** This is at which height to measure the trunk diameter. The default setting is from 1.2m to 1.4m.
     * **Clustering Parameter:** You can adjust the clustering density by epsilon distance. The default setting is 0.3m and 10 pts for a cluster.

## Remind
* All the point cloud trees should be cleaned as individual trees. Any attachments on trunks should be removed, as they will be counted as part of the trunk.
* The point cloud trees export can be operated in Cloud Compare. After segmentation and cleaning, run cloud_compare.py to export the results to the folder. Then run the delete.py to remove the noise clouds since a complete tree cloud is usually >100kb at least.

## Citiation
* Yao, C. and P. Fricker, Building Green Decarbonization for Urban Digital Twin – Estimating Carbon Sequestration of Urban Trees by Allometric Equations using Blend Types of Point Cloud. Proceedings of the XXVII International Conference of the Ibero-American Society of Digital Graphics, 2023(27): p. pp. 91–102.


## 🛠️ Installation

### Prerequisites
* Python 3.8+
* Anaconda (recommended)

### Dependencies
Install the required libraries using pip:

```bash
pip install numpy scipy open3d
