# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:16:46 2023

@author: Sleepyhead

Tree DBH metrics batch processing with crown diameter and global summary
"""

import os
import csv
import math
import numpy as np
import open3d as o3d
from scipy.spatial import ConvexHull

# =============================================================================
# --- CONFIGURATION -----------------------------------------------------------
# =============================================================================

# PATHS
INPUT_ROOT = r"F:/mobile/F/Tree"          # Folder containing species subfolders
OUTPUT_ROOT = r"F:/mobile/F/Tree_Results" # Folder where CSV results will be saved

# TRUNK SLICING PARAMETERS
DBH_CENTER_H = 1.30        # Height at which to measure DBH (meters)
SLICE_HALF_THICKNESS = 0.10 # Thickness of slice (+/- from center height)

# TRUNK EXTRACTION FILTERS
TRUNK_ZONE_LOW = 0.00      # Lower bound for trunk isolation (meters from min Z)
TRUNK_ZONE_HIGH = 1.50     # Upper bound for trunk isolation (meters from min Z)

# CLUSTERING PARAMETERS (DBSCAN)
DB_EPS = 0.3               # Epsilon distance for clustering (meters)
DB_MIN_PTS = 10            # Minimum points to form a cluster

# =============================================================================
# --- HELPER FUNCTIONS --------------------------------------------------------
# =============================================================================

def circle_ls(points_xy):
    """Least-squares circle fitting."""
    mean = np.mean(points_xy, axis=0)
    P = points_xy - mean
    x, y = P[:, 0], P[:, 1]
    A = np.c_[2*x, 2*y, np.ones(len(x))]
    b = x**2 + y**2
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    xc_rel, yc_rel, c2 = c[0], c[1], c[2]
    r = math.sqrt(max(0.0, c2 + xc_rel**2 + yc_rel**2))
    center = np.array([xc_rel, yc_rel]) + mean
    return center, r

def mec_minidisk(points_xy):
    """Minimum Enclosing Circle (Welzl's algorithm)."""
    import random
    P = points_xy.copy().tolist()
    random.shuffle(P)
    def circ_from_two(a, b):
        cx = 0.5*(a[0] + b[0])
        cy = 0.5*(a[1] + b[1])
        r = math.hypot(a[0]-cx, a[1]-cy)
        return (cx, cy), r
    def circ_from_three(a, b, c):
        ax, ay = a; bx, by = b; cx_, cy_ = c
        d = 2*(ax*(by-cy_) + bx*(cy_-ay) + cx_*(ay-by))
        if abs(d) < 1e-12:
            return None
        ux = ((ax*ax+ay*ay)*(by-cy_) + (bx*bx+by*by)*(cy_-ay) + (cx_*cx_+cy_*cy_)*(ay-by)) / d
        uy = ((ax*ax+ay*ay)*(cx_-bx) + (bx*bx+by*by)*(ax-cx_) + (cx_*cx_+cy_*cy_)*(bx-ax)) / d
        r = math.hypot(ux-ax, uy-ay)
        return (ux, uy), r
    def inside(pt, center, r):
        return math.hypot(pt[0]-center[0], pt[1]-center[1]) <= r + 1e-9
    def welzl(P):
        center = (0.0, 0.0); r = -1.0
        for i, p in enumerate(P):
            if r >= 0 and inside(p, center, r):
                continue
            center = (p[0], p[1]); r = 0.0
            for j, q in enumerate(P[:i]):
                if inside(q, center, r):
                    continue
                center, r = circ_from_two(p, q)
                for s in P[:j]:
                    if inside(s, center, r):
                        continue
                    res = circ_from_three(p, q, s)
                    if res is not None:
                        center, r = res
        return center, r
    c, r = welzl(P)
    return np.array(c), r

def fit_cylinder_internal(points):
    """Fits cylinder based on mean distance to axis."""
    mu = points.mean(axis=0)
    X = points - mu
    cov = np.cov(X.T)
    vals, vecs = np.linalg.eigh(cov)
    axis_dir = vecs[:, np.argmax(vals)]
    if axis_dir[2] < 0:
        axis_dir = -axis_dir
    t = (points - mu) @ axis_dir
    axis_points = mu + np.outer(t, axis_dir)
    dists = np.linalg.norm(points - axis_points, axis=1)
    radius = dists.mean()
    return mu, axis_dir, radius

def fit_cylinder_external(points):
    """Fits cylinder based on max distance to axis (enclosing)."""
    mu = points.mean(axis=0)
    X = points - mu
    cov = np.cov(X.T)
    vals, vecs = np.linalg.eigh(cov)
    axis_dir = vecs[:, np.argmax(vals)]
    if axis_dir[2] < 0:
        axis_dir = -axis_dir
    t = (points - mu) @ axis_dir
    axis_points = mu + np.outer(t, axis_dir)
    dists = np.linalg.norm(points - axis_points, axis=1)
    radius = dists.max()
    return mu, axis_dir, radius

def fit_external_elliptical_cylinder(points_xy):
    """Fits an ellipse to the 2D projection of points."""
    mu = points_xy.mean(axis=0)
    X = points_xy - mu
    cov = np.cov(X.T)
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vecs = vecs[:, order]
    Xp = X @ vecs
    x_min, x_max = Xp[:,0].min(), Xp[:,0].max()
    y_min, y_max = Xp[:,1].min(), Xp[:,1].max()
    major_d = (x_max - x_min)
    minor_d = (y_max - y_min)
    eq_d = math.sqrt(major_d * minor_d)
    return major_d, minor_d, eq_d

def trunk_orientation(trunk_pts):
    """Calculates the lean angle of the trunk in degrees."""
    mu = trunk_pts.mean(axis=0)
    X = trunk_pts - mu
    cov = np.cov(X.T)
    vals, vecs = np.linalg.eigh(cov)
    axis_dir = vecs[:, np.argmax(vals)]
    if axis_dir[2] < 0:
        axis_dir = -axis_dir
    angle = math.degrees(math.acos(abs(axis_dir[2])))
    return angle

def coverage_angle(points_xy):
    """Calculates coverage angle using Convex Hull."""
    if len(points_xy) < 3:
        return 0.0
    hull = ConvexHull(points_xy)
    hull_pts = points_xy[hull.vertices]
    cx, cy = hull_pts[:,0].mean(), hull_pts[:,1].mean()
    angles = np.degrees(np.arctan2(hull_pts[:,1]-cy, hull_pts[:,0]-cx))
    angles = np.sort((angles + 360) % 360)
    diffs = np.diff(np.r_[angles, angles[0]+360])
    coverage = 360 - diffs.max()
    return coverage

def bbox_xy(points_xy):
    """Calculates 2D Bounding Box dimensions."""
    xmin, ymin = points_xy.min(axis=0)
    xmax, ymax = points_xy.max(axis=0)
    return xmax-xmin, ymax-ymin

# =============================================================================
# --- MAIN PROCESSING LOOP ----------------------------------------------------
# =============================================================================

def main():
    # Create output folder
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    # Get species folders
    species_folders = [f for f in os.listdir(INPUT_ROOT) 
                       if os.path.isdir(os.path.join(INPUT_ROOT, f))]
    
    print(f"Found {len(species_folders)} species folders in {INPUT_ROOT}")

    global_results = []

    for species in species_folders:
        species_path = os.path.join(INPUT_ROOT, species)
        output_csv = os.path.join(OUTPUT_ROOT, f"{species}_dbh_metrics.csv")
        
        files = [f for f in os.listdir(species_path) 
                 if os.path.isfile(os.path.join(species_path, f))]
        
        results = []
        
        print(f"\n--- Processing Species: {species} ({len(files)} files) ---")

        for file in files:
            try:
                path = os.path.join(species_path, file)
                pcd = o3d.io.read_point_cloud(path)
                pts = np.asarray(pcd.points)
                if pts.size == 0: continue
                
                z = pts[:,2]
                z0 = z.min()
                tree_height = z.max() - z.min() + 1
                
                # Compute crown diameter from full tree BEFORE trunk slicing
                crown_XY = pts[:, :2]
                crown_bbox_x, crown_bbox_y = bbox_xy(crown_XY)
                crown_d = math.sqrt(crown_bbox_x * crown_bbox_y)
                
                # Trunk zone filter
                mask_zone = (z >= z0 + TRUNK_ZONE_LOW) & (z <= z0 + TRUNK_ZONE_HIGH)
                zone_pts = pts[mask_zone]
                if len(zone_pts) < 30: continue
                
                # DBSCAN to isolate trunk
                zone_pcd = o3d.geometry.PointCloud()
                zone_pcd.points = o3d.utility.Vector3dVector(zone_pts)
                labels = np.array(zone_pcd.cluster_dbscan(eps=DB_EPS, min_points=DB_MIN_PTS, print_progress=False))
                if labels.size == 0 or labels.max() < 0: continue
                
                tallest_lbl, tallest_h = None, -1
                for lbl in np.unique(labels):
                    if lbl == -1: continue
                    cp = zone_pts[labels == lbl]
                    h = cp[:,2].max() - cp[:,2].min()
                    if h > tallest_h:
                        tallest_h = h
                        tallest_lbl = lbl
                if tallest_lbl is None: continue
                trunk_pts = zone_pts[labels == tallest_lbl]
                
                angle = trunk_orientation(trunk_pts)
                
                # DBH slice extraction
                h1 = z0 + (DBH_CENTER_H - SLICE_HALF_THICKNESS)
                h2 = z0 + (DBH_CENTER_H + SLICE_HALF_THICKNESS)
                mask_slice = (trunk_pts[:,2] >= h1) & (trunk_pts[:,2] <= h2)
                slice_pts = trunk_pts[mask_slice]
                if len(slice_pts) < 10: continue
                
                XY = slice_pts[:, :2]
                
                # Metric Calculations
                ls_center, ls_r = circle_ls(XY); ls_d = 2*ls_r
                mec_center, mec_r = mec_minidisk(XY); mec_d = 2*mec_r
                cyl_mu, cyl_dir, cyl_r_int = fit_cylinder_internal(slice_pts); cyl_d_int = 2*cyl_r_int
                cyl_mu, cyl_dir, cyl_r_ext = fit_cylinder_external(slice_pts); cyl_d_ext = 2*cyl_r_ext
                ell_major_d, ell_minor_d, ell_eq_d = fit_external_elliptical_cylinder(XY)
                cover = coverage_angle(XY)
                bbox_x, bbox_y = bbox_xy(XY)
                
                row = {
                    "filename": file,
                    "slice_points": len(slice_pts),
                    "coverage_deg": round(float(cover),1),
                    "bbox_x_m": round(float(bbox_x),4),
                    "bbox_y_m": round(float(bbox_y),4),
                    "ls_d_m": round(float(ls_d),4),
                    "mec_d_m": round(float(mec_d),4),
                    "cyl_d_int_m": round(float(cyl_d_int),4),
                    "cyl_d_ext_m": round(float(cyl_d_ext),4),
                    "ell_major_d_m": round(float(ell_major_d),4),
                    "ell_minor_d_m": round(float(ell_minor_d),4),
                    "ell_eq_diameter_m": round(float(ell_eq_d),4),
                    "trunk_angle_deg": round(float(angle),2),
                    "tree_height_m": round(float(tree_height),2),
                    "crown_diameter_m": round(float(crown_d),4)
                }
                
                results.append(row)
                global_results.append({"species": species, **row})
                
                print(f"  > {file}: EQ={ell_eq_d:.3f} m, Crown={crown_d:.3f} m")
                
            except Exception as e:
                print(f"  ! ERROR {file}: {e}")
        
        # Save CSV for this species
        with open(output_csv, "w", newline="") as f:
            fieldnames = list(results[0].keys()) if results else []
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in results:
                w.writerow(r)
        
        print(f"  Saved metrics for {species} -> {output_csv}")

    # Save Global Summary
    global_csv = os.path.join(OUTPUT_ROOT, "all_species_summary.csv")
    if global_results:
        with open(global_csv, "w", newline="") as f:
            fieldnames = list(global_results[0].keys())
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in global_results:
                w.writerow(r)
        print(f"\nSaved global summary -> {global_csv}")
    else:
        print("\nNo results to save.")

if __name__ == "__main__":
    main()