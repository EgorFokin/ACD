from collections import defaultdict

def orientation(p, q, r):
    """Return orientation of triplet (p,q,r)."""
    val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
    if val > 0:
        return 1  # clockwise
    elif val < 0:
        return 2  # counterclockwise
    else:
        return 0  # collinear

def on_segment(p, q, r, eps=1e-5):
    """Check if point q lies on segment pr."""
    return (min(p[0], r[0]) - eps <= q[0] <= max(p[0], r[0]) + eps and
            min(p[1], r[1]) - eps <= q[1] <= max(p[1], r[1]) + eps)

def edges_intersect(a,b,c,d):
    """Check if segment ab intersects cd."""


    o1 = orientation(a, b, c)
    o2 = orientation(a, b, d)
    o3 = orientation(c, d, a)
    o4 = orientation(c, d, b)

    x1, y1 = a
    x2, y2 = b
    x3, y3 = c
    x4, y4 = d

    denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    # General case
    if o1 != o2 and o3 != o4 and abs(denom) > 1e-5:
        return True


    return False

def point_distance(a, b):
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2)**0.5

def intersection_point(a, b, c, d):
    """
    Get intersection point of two segments ab and cd.
    Returns:
      (x,y) if they intersect at a single point,
      (p,q) if they overlap in a segment (collinear),
    """

    x1, y1 = a
    x2, y2 = b
    x3, y3 = c
    x4, y4 = d

    denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)


    if abs(denom) > 1e-12:
        # Proper intersection point
        px = ((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4)) / denom
        py = ((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4)) / denom
        return (px, py)

    return False


        


