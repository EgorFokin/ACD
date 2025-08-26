import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import os
from load_figure import load_figure, visualize_figure
from skimage import measure
import argparse

os.environ["QT_QPA_PLATFORM"] = "xcb"

#1.load figure
#2.compute SDF
#3.compute shell

def point_segment_distance(q, p1, p2):
    v = p2 - p1
    w = q - p1
    t = np.clip(np.dot(w, v) / np.dot(v, v), 0.0, 1.0)
    proj = p1 + t * v
    return np.linalg.norm(q - proj)

def extract_contour(sdf, xs, ys, level=0.0):
    # find iso-contour (zero level set)
    contours = measure.find_contours(sdf, level=level)

    remeshed = []
    for contour in contours:
        # convert pixel coords to real coords
        c = np.array([
            [xs[int(x)], ys[int(y)]]
            for y, x in contour  # note order swap
        ])
        remeshed.append(c)
    return remeshed

def contour_to_graph(contours):
    all_points = []
    all_edges = []
    offset = 0
    for c in contours:
        n = len(c)
        all_points.extend(c)
        all_edges.extend([(offset+i, offset+(i+1)%n) for i in range(n)])
        offset += n
    return np.array(all_points), np.array(all_edges)

def compute_sdf(vertices, edges, res=200):
    xs = np.linspace(vertices[:,0].min()-1, vertices[:,0].max()+1, res)
    ys = np.linspace(vertices[:,1].min()-1, vertices[:,1].max()+1, res)
    X, Y = np.meshgrid(xs, ys)

    # Build polygons from edges
    # Assume edges form closed loops, find all unique loops
    polygons = []

    # A naive approach: collect connected components from edges
    # Here we assume edges are ordered per polygon
    visited = set()
    for i, e in enumerate(edges):
        if i in visited:
            continue
        # trace a loop
        loop = [e[0], e[1]]
        visited.add(i)
        current = e[1]
        while loop[0] != current:
            for j, f in enumerate(edges):
                if j in visited:
                    continue
                if f[0] == current:
                    loop.append(f[1])
                    current = f[1]
                    visited.add(j)
                    break
                elif f[1] == current:
                    loop.append(f[0])
                    current = f[0]
                    visited.add(j)
                    break
        polygons.append(Polygon(vertices[loop]))

    sdf = np.zeros_like(X, dtype=float)
    for i in range(res):
        for j in range(res):
            q = np.array([X[i,j], Y[i,j]])
            # distance to edges
            dists = [point_segment_distance(q, vertices[e[0]], vertices[e[1]]) for e in edges]
            d = min(dists)
            # sign: negative if inside any polygon
            if any(p.contains(Point(q)) for p in polygons):
                d = -d
            sdf[i,j] = d
    return xs, ys, sdf

def compute_shell(sdf, xs, ys, level=0.1):
    contours = extract_contour(sdf, xs, ys, level=level)
    vertices, edges = contour_to_graph(contours)
    return vertices, edges

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute SDF from SVG figure.")
    parser.add_argument("--svg", type=str, default="Pi-symbol.svg", help="Path to SVG file.")
    parser.add_argument("--res", type=int, default=100, help="Resolution for SDF computation.")
    parser.add_argument("--cached", action="store_true", help="Use cached SDF if available.")
    parser.add_argument("--level", type=float, default=0.1, help="Contour level for extraction.")
    parser.add_argument("--samples", type=int, default=20, help="Number of samples per segment.")
    args = parser.parse_args()

    vertices, edges = load_figure(args.svg,samples_per_segment=args.samples)

    print(len(edges))

    if args.cached and os.path.exists(os.path.join('cache', 'sdf.npy')):
        sdf = np.load(os.path.join('cache', 'sdf.npy'))
        xs = np.linspace(min(vertices[:,0])-1, max(vertices[:,0])+1, args.res)
        ys = np.linspace(min(vertices[:,1])-1, max(vertices[:,1])+1, args.res)
    else:
        xs, ys, sdf = compute_sdf(vertices, edges, res=args.res)
        np.save(os.path.join('cache', 'sdf.npy'), sdf)


    plt.figure(figsize=(6,6))
    visualize_figure(vertices, edges)
    plt.imshow(sdf, extent=(xs.min(), xs.max(), ys.min(), ys.max()),
            origin='lower', cmap='RdBu_r')
    plt.colorbar(label='Signed Distance')
    plt.title('Signed Distance Field')
    plt.show()

    remeshed_vertices, remeshed_edges = compute_shell(sdf, xs, ys, level=args.level)

    plt.figure(figsize=(6, 6))
    visualize_figure(vertices, edges)
    visualize_figure(remeshed_vertices, remeshed_edges)

    plt.show()