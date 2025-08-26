import matplotlib.pyplot as plt
import numpy as np
import os
from shapely.geometry import Point, Polygon
from collections import defaultdict

from scipy.spatial import ConvexHull

from utils import *

os.environ["QT_QPA_PLATFORM"] = "xcb"

def point_in_hull(point, hull, tol=1e-4):
    """
    Check if point is inside convex hull.
    point: (x,y)
    hull: scipy.spatial.ConvexHull
    """
    return np.all(np.dot(hull.equations[:,:-1], point) + hull.equations[:,-1] <= tol)


def point_in_mesh(point, sdf):
    """
    Check if point is inside mesh using SDF.
    point: (x,y)
    sdf: 2D array of signed distance values
    """
    x, y = point
    h, w = sdf.shape
    sdf_xmin, sdf_xmax = -2, 2  
    sdf_ymin, sdf_ymax = -2, 2  
    ix = int((x - sdf_xmin) / (sdf_xmax - sdf_xmin) * w)
    iy = int((y - sdf_ymin) / (sdf_ymax - sdf_ymin) * h)
    if ix < 0 or ix >= w or iy < 0 or iy >= h:
        return False
    return sdf[iy, ix] < -0.01


def substract_hull(vertices, mesh_edges, hull_edges, hull, sdf):

    hull_verts = set()
    for e in hull_edges:
        hull_verts.add(e[0])
        hull_verts.add(e[1])

    new_verts = vertices[:]
    new_edges = mesh_edges.copy()
    remove_edges = set()

    for e in mesh_edges:
        if (min(e[0], e[1]), max(e[0], e[1])) in hull_edges:
            center_point = (vertices[e[0]] + vertices[e[1]]) / 2
            if point_in_mesh(center_point, sdf):
                continue
        for he in hull_edges:
            if point_in_hull(vertices[e[0]], hull) and point_in_hull(vertices[e[1]], hull):
                remove_edges.add(e)
                break
            if edges_intersect(vertices[e[0]], vertices[e[1]], new_verts[he[0]], new_verts[he[1]]):
                #add intersection to new_verts

                intersect_point = intersection_point(vertices[e[0]], vertices[e[1]], new_verts[he[0]], new_verts[he[1]])
                if point_distance(intersect_point, vertices[e[0]]) < 1e-5 or point_distance(intersect_point, vertices[e[1]]) < 1e-5:
                    continue

                hull_edges.remove(he)
                hull_edges.add((he[0], len(new_verts)))
                hull_edges.add((he[1], len(new_verts)))


                center_point1 = (new_verts[he[0]] + intersect_point) / 2
                center_point2 = (new_verts[he[1]] + intersect_point) / 2
                if point_in_mesh(center_point1, sdf):
                    new_edges.add((he[0], len(new_verts)))
                if point_in_mesh(center_point2, sdf):
                    new_edges.add((he[1], len(new_verts)))

                new_verts = np.vstack([new_verts, intersect_point])

                

                
                if point_in_hull(vertices[e[0]], hull):
                    # if 20 in e and 114 in e:
                    #     print(point_in_hull(vertices[e[0]], hull))
                    #     print(point_in_hull(vertices[e[1]], hull))
                    new_edges.add((e[1], len(new_verts)-1))
                else:
                    new_edges.add((e[0], len(new_verts)-1))
                remove_edges.add(e)
                break
                
            


    # for e in list(mesh_edges):
    #     if point_in_hull(vertices[e[0]], hull) or point_in_hull(vertices[e[1]], hull):
    #         remove_edges.add(e)
                
    
    new_edges -= remove_edges

    return new_verts, new_edges 


def visualize_inside_points(vertices, hull, hull_edges):

    hull_verts = set()
    for e in hull_edges:
        hull_verts.add(e[0])
        hull_verts.add(e[1])

    inside_points = []
    outside_points = []
    for i, v in enumerate(vertices):
        if point_in_hull(v, hull):
            inside_points.append(v)
        else:
            outside_points.append(v)
    
    inside_points = np.array(inside_points)
    outside_points = np.array(outside_points)

    plt.scatter(inside_points[:,0], inside_points[:,1], c='red', s=2)
    plt.scatter(outside_points[:,0], outside_points[:,1], c='blue', s=2)
    plt.show()


if __name__ == "__main__":
    vertices = np.load("vertices.npy")
    triangles = np.load("triangles.npy")
    segment = np.load("segment.npy")

    mapping = {new: orig for new, orig in enumerate(segment)}
    hull = ConvexHull(vertices[segment], qhull_options='QJ')
    hull_verts = np.asarray([mapping[v] for v in hull.vertices])

    plt.figure(figsize=(6, 6))
    plt.xlim((-2,2))
    plt.ylim((-2,2))

    


    # plt.triplot(vertices[:,0], vertices[:,1], triangles, color='lightgray')


    
    hull_points = vertices[hull_verts]
    # plt.fill(hull_points[:,0], hull_points[:,1], alpha=0.3)
    # plt.show()

    mesh_edges = set()
    for tri in triangles:
        mesh_edges.add((min(tri[0], tri[1]), max(tri[0], tri[1])))
        mesh_edges.add((min(tri[1], tri[2]), max(tri[1], tri[2])))
        mesh_edges.add((min(tri[2], tri[0]), max(tri[2], tri[0])))
    
    hull_edges = set()
    for i in range(len(hull_verts)):
        hull_edges.add((min(hull_verts[i], hull_verts[(i+1)%len(hull_verts)]), max(hull_verts[i], hull_verts[(i+1)%len(hull_verts)])))

    # visualize_inside_points(vertices, hull, hull_edges)

    print("before:", len(mesh_edges))

    


    new_verts, new_edges = substract_hull(vertices, mesh_edges, hull_edges, hull)

    print("after:", len(new_edges))

    #plot edges
    for e in new_edges:
        plt.plot([new_verts[e[0],0], new_verts[e[1],0]], [new_verts[e[0],1], new_verts[e[1],1]], c='black', lw=1)

    colors = np.asarray([1,0,0]*len(new_verts)).reshape(-1,3)
    #color 20 and 114 in blue
    colors[20] = [0,0,1]
    plt.scatter(new_verts[:,0], new_verts[:,1], s=2, c=colors)
    plt.show()





