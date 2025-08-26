import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from load_figure import load_figure, visualize_figure
from compute_shell import compute_sdf, compute_shell
from tqdm import tqdm
import PythonCDT as cdt
from collections import defaultdict
from shapely.geometry import Polygon
from utils import *
from difference import substract_hull

from scipy.spatial import ConvexHull


def intersect_shell(v1, v2, shell):
    """Check if segment v1-v2 intersects with any edge in the shell."""
    for e in shell['edges']:
        if edges_intersect(v1, v2, shell['vertices'][e[0]], shell['vertices'][e[1]]):
            return True
    return False

def build_adjacency(edges):
    """Build adjacency list from edge list."""
    adj = defaultdict(set)
    for e in edges:
        adj[e[0]].add(e[1])
        adj[e[1]].add(e[0])
    return adj


def decompose_iteration(used, adj, verts, shell):
    """Perform one iteration of the decomposition."""
    # Only check neighbors of already used vertices
    candidates = set()
    
    for v in used:
        candidates.update(adj[v])

    candidates -= set(used)

    
    for vi in list(candidates):
        flag = True
        for vj in used:
            if intersect_shell(verts[vi], verts[vj], shell):
                flag = False
                break
        if flag:
            used.append(vi)
            return True
    return False


def decompose(mesh, shell, sdf):
    

    # build adjacency from triangles

    n_vertices = len(mesh['vertices'])
    unused = set(range(n_vertices))
    hulls = []

    
    while unused:
        segment = [unused.pop()]
        adj = build_adjacency(mesh['edges'])
        
        while True:
            flag = decompose_iteration(segment, adj, mesh['vertices'], shell)
            if not flag:
                break
            print(segment)

        unused -= set(segment)

        if len(segment) < 3:
            continue
            

        
        mapping = {new: orig for new, orig in enumerate(segment)}
        hull = ConvexHull(mesh['vertices'][segment], qhull_options='QJ')
        hull_verts = np.asarray([mapping[v] for v in hull.vertices])

        hull_edges = set()
        for i in range(len(hull_verts)):
            hull_edges.add((min(hull_verts[i], hull_verts[(i+1)%len(hull_verts)]), max(hull_verts[i], hull_verts[(i+1)%len(hull_verts)])))

        new_verts, new_edges = substract_hull(mesh['vertices'], mesh['edges'], hull_edges, hull, sdf)
        
        mesh['vertices'] = new_verts
        mesh['edges'] = new_edges

        hull_points = mesh['vertices'][hull_verts]

        visualize_figure(mesh['vertices'], list(mesh['edges']))
        plt.fill(hull_points[:,0], hull_points[:,1], alpha=0.3)
        plt.show()

        hulls.append(hull_points)

        

        


        # if len(segments) == 2:
        #     break

    return hulls



def triangulate_polygon(points, edges):
    t = cdt.Triangulation(cdt.VertexInsertionOrder.AS_PROVIDED, cdt.IntersectingConstraintEdges.TRY_RESOLVE, 0.0)

    cdt_vertices = []
    cdt_edges = []
    for v in points:
        cdt_vertices.append(cdt.V2d(v[0], v[1]))

    for e in edges:
        cdt_edges.append(cdt.Edge(e[0],e[1]))

    t.insert_vertices(cdt_vertices)
    t.insert_edges(cdt_edges)

    t.erase_outer_triangles_and_holes()

    triangles = [list(tri.vertices) for tri in  t.triangles]
    return triangles


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute shell from SVG figure.")
    parser.add_argument("--svg", type=str, default="Pi-symbol.svg", help="Path to the SVG file.")
    parser.add_argument("--res", type=int, default=100, help="Resolution for SDF computation.")
    parser.add_argument("--level", type=float, default=0.1, help="Level for contour extraction.")
    parser.add_argument("--cached", action='store_true', help="Use cached SDF if available.")
    parser.add_argument("--samples", type=int, default=20, help="Number of samples per segment.")   
    
    args = parser.parse_args()

    # Load figure
    vertices, edges = load_figure(args.svg, samples_per_segment=args.samples)

    # Compute SDF
    if args.cached and os.path.exists(os.path.join('cache', 'sdf.npy')):
        sdf = np.load(os.path.join('cache', 'sdf.npy'))
        xs = np.linspace(min(vertices[:,0])-1, max(vertices[:,0])+1, args.res)
        ys = np.linspace(min(vertices[:,1])-1, max(vertices[:,1])+1, args.res)
    else:
        xs, ys, sdf = compute_sdf(vertices, edges, res=args.res)
        np.save(os.path.join('cache', 'sdf.npy'), sdf)

    shell_vertices, shell_edges = compute_shell(sdf, xs, ys, level=args.level)

    # visualize_figure(shell_vertices, shell_edges)

    # print(vertices)
    vertices = np.array(vertices)
    edges = np.array(edges)

    

    triangles = triangulate_polygon(vertices, edges)

   
    # plt.show()
    
    mesh_edges = set()
    for tri in triangles:
        mesh_edges.add((min(tri[0], tri[1]), max(tri[0], tri[1])))
        mesh_edges.add((min(tri[1], tri[2]), max(tri[1], tri[2])))
        mesh_edges.add((min(tri[2], tri[0]), max(tri[2], tri[0])))




    mesh = {'vertices': vertices, 'edges': edges, 'edges': mesh_edges}
    shell = {'vertices': shell_vertices, 'edges': shell_edges}


    # if args.cached and os.path.exists(os.path.join('cache', 'graph.npy')):
    #     g = np.load(os.path.join('cache', 'graph.npy'))
    # else:
    #     g = compute_graph(mesh, shell)
    #     np.save(os.path.join('cache', 'graph.npy'), g)

    hulls = decompose(mesh, shell, sdf)
    print(len(hulls))

    plt.figure(figsize=(6, 6))
    plt.xlim((-2,2))
    plt.ylim((-2,2))

    # visualize_figure(mesh['vertices'], list(mesh['edges']))
    for h in hulls:
        plt.fill(h[:,0], h[:,1], alpha=0.3)
    

    # x = vertices[:, 0]
    # y = vertices[:, 1]
    # plt.triplot(x, y, triangles, color='gray', alpha=0.5)
    # colors = np.asarray([1,0,0]*len(vertices)).reshape(-1,3)
    # colors[list(border_verts)] = [0,0,1]
    # plt.scatter(x, y, color=colors, s=10)


    # segment = segments[1]
    # np.save('segment.npy', segment)
    # mapping = {new: orig for new, orig in enumerate(segment)}
    # hull = ConvexHull(mesh['vertices'][segment], qhull_options='QJ')
    # hull_verts = np.asarray([mapping[v] for v in hull.vertices])

    
    # hull_points = mesh['vertices'][hull_verts]
    # plt.fill(hull_points[:,0], hull_points[:,1], alpha=0.3)

    # np.save('hull.npy', hull_verts)
    # np.save('vertices.npy', mesh['vertices'])
    # np.save('triangles.npy', mesh['triangles'])
    
    plt.show()

    # colors = np.zeros(len(mesh['vertices']), dtype=int)
    # for i, seg in enumerate(segments):
    #     for v in seg:
    #         colors[v] = i % 10
    


    
    # for i, j in edges:
    #     plt.plot([vertices[i, 0], vertices[j, 0]], 
    #              [vertices[i, 1], vertices[j, 1]], c="black", lw=1)
    # # Draw vertices
    # plt.scatter(vertices[:, 0], vertices[:, 1], 
    #         s=50, c=colors, cmap="tab10", edgecolors="k")
    
    # visualize_figure(shell_vertices, shell_edges)

    




