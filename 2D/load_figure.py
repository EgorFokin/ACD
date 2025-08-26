import matplotlib.pyplot as plt
from svgpathtools import svg2paths
import numpy as np
import os
from scipy.spatial import cKDTree
import argparse

os.environ["QT_QPA_PLATFORM"] = "xcb"


def merge_duplicate_points(points, edges, tol=1e-5):
    """
    Merge points that are closer than `tol` and update edges.

    Parameters
    ----------
    points : np.ndarray, shape (N,2)
        Array of points.
    edges : list of tuple
        List of edges as index pairs.
    tol : float
        Distance threshold for merging points.

    Returns
    -------
    new_points : np.ndarray
        Merged points.
    new_edges : list of tuple
        Updated edges with merged point indices.
    """
    tree = cKDTree(points)
    groups = tree.query_ball_tree(tree, r=tol)

    # Map old indices to new indices
    index_map = {}
    visited = set()
    new_points = []
    for i, group in enumerate(groups):
        if i in visited:
            continue
        # Take the first point in the group as representative
        rep_idx = len(new_points)
        new_points.append(points[i])
        for j in group:
            index_map[j] = rep_idx
            visited.add(j)

    new_points = np.array(new_points)
    new_edges = [(index_map[i], index_map[j]) for i, j in edges if index_map[i] != index_map[j]]

    return new_points, new_edges

def normalize_points(points):
    """
    Normalize points to fit within a unit square centered at the origin.

    Parameters
    ----------
    points : np.ndarray
        Array of points to normalize, shape (N, 2).

    Returns
    -------
    normalized_points : np.ndarray
        Normalized points, shape (N, 2).
    """
    min_x, min_y = np.min(points, axis=0)
    max_x, max_y = np.max(points, axis=0)
    
    width = max_x - min_x
    height = max_y - min_y
    
    normalized_points = (points - [min_x, min_y]) / [width, height]
    normalized_points = normalized_points * 2 - 1  # Scale to [-1, 1]
    
    return normalized_points

def load_figure(svg_path, samples_per_segment=20, normalize=True):
    """
    Load an SVG file and convert it into a set of points and edges.

    Parameters
    ----------
    svg_path : str
        Path to the SVG file.
    samples_per_segment : int, optional
        Number of points to sample per segment (default=20).
    normalize : bool, optional
        Whether to normalize the points to fit within a unit square (default=True).

    Returns
    -------
    points : np.ndarray
        Array of sampled points, shape (N, 2).
    edges : list of tuple
        List of edges as (i, j) index pairs into points.
    """
    paths, attributes = svg2paths(svg_path)

    points = []
    edges = []

    for path in paths:
        for segment in path:
            segment_points = []
            for t in np.linspace(0, 1, samples_per_segment):
                point = segment.point(t)
                segment_points.append((point.real, -point.imag))  # Flip Y

            # Add edges between consecutive sampled points
            start_idx = len(points)
            points.extend(segment_points)
            for i in range(len(segment_points) - 1):
                edges.append((start_idx + i, start_idx + i + 1))

    points = np.array(points)   

    points, edges = merge_duplicate_points(points, edges)

    if normalize:
        points = normalize_points(points)

    return points, edges

def visualize_figure(points, edges):
    # Draw edges
    for i, j in edges:
        plt.plot([points[i, 0], points[j, 0]], 
                 [points[i, 1], points[j, 1]], c="black", lw=1)
    # Draw vertices
    plt.scatter(points[:, 0], points[:, 1], s=2, c="red")
    plt.axis("equal")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load and visualize an SVG figure.")
    parser.add_argument("--svg", type=str, default="Pi-symbol.svg", help="Path to the SVG file.")
    parser.add_argument("--samples", type=int, default=20, help="Number of samples per segment.")
    args = parser.parse_args()


    points, edges = load_figure(args.svg)
    
    plt.figure(figsize=(6, 6))
    visualize_figure(points, edges)

    plt.show()