# SPDX-License-Identifier: LGPL-2.1-or-later

__title__ = "Surface tangency solver"
__author__ = "Andrew Shkolik"
__license__ = "LGPL 2.1+"

from collections import defaultdict
import math
from time import time
import FreeCAD as App
import Part
import math

### Debug visualization utilities ###

def show_object(shape, name="DebugShape", color=(0.0,1.0,0.0)):
    """
    Show a shape in FreeCAD for debugging.

    Args:
        shape: Part.Shape
        name: Name for the object
        color: Color (r,g,b)
    Returns:
        Part.Feature or None: The created object or None if not shown
    """
    obj = App.ActiveDocument.addObject("Part::Feature", name)
    obj.Shape = shape
    obj.ViewObject.ShapeColor = color
    obj.ViewObject.LineColor = color
    obj.ViewObject.PointColor = color
    return obj

def show_surface(bs, name="DebugSurface", color=None):
    """
    Show a BSpline surface in FreeCAD for debugging.

    Args:
        bs: Part.BSplineSurface
        name: Name for the object
        color: Color (r,g,b)
    Returns:
        Part.Feature or None: The created object or None if not shown
    """
    obj = App.ActiveDocument.addObject("Part::Feature", name)
    obj.Shape = Part.Face(bs.copy())
    if color:
        obj.ViewObject.ShapeColor = color
    return obj

def draw_point(point, color=(1.0,0,0), size=0.2, name="Pt", show=True):    
    """
    Method to draw points as spheres in FreeCAD for visualization.

    Args:
        point: point (App.Vector)
        color: Sphere color (r,g,b)
        size: Sphere radius
        name: Name for the object
        show: Whether to create the object in FreeCAD document
    Returns:
        Part.Shape: The created sphere
    """    
    p = Part.makeSphere(size, point)
    
    if show:
        show_object(p, name=name, color=color)
    return p

def draw_vector(point, vector, color=(1.0,0,0), scale=10.0, name="Vec", show=True):
    """
    Method to draw a vector as a line in FreeCAD for visualization.
    
    Args:
        point: Starting point (App.Vector)
        vector: Vector to draw (App.Vector)
        color: Line color (r,g,b)
        scale: Scale factor for vector length if vector is unit length
        name_prefix: Name prefix for the object
        show: Whether to create the object in FreeCAD document
    Returns:
        Part.Edge or None: The created edge or None if vector is invalid
    """
    if vector is None or vector.Length < 1e-9:
        return None
    
    end = App.Vector(vector)
    if end.Length <= 1:
        # end = vector.normalize().multiply(scale)
        end.normalize()
        end.multiply(scale)

    p1 = point
    p2 = point.add(end)
    edge = Part.makeLine(p1, p2)
    p = draw_point(p2, color=color, show=False)
    comp = Part.makeCompound([edge, p])
    if show:
        show_object(comp, name, color=color)
    return comp

def draw_vectors(points, vectors, color=(1.0,0,0), scale=10.0, name="Vectors", show=True):    
    """
    Method to draw multiple vectors as lines in FreeCAD for visualization.

    Args:
        points: List of starting points (App.Vector)
        vectors: List of vectors to draw (App.Vector)
        color: Line color (r,g,b)
        scale: Scale factor for vector length if vector is unit length
        name: Name for the object
        show: Whether to create the object in FreeCAD document
    Returns:
        Part.Compound or None: The Compound of vectors
    """
    vecs = []
    for point, vector in zip(points, vectors):
        vec = draw_vector(point, vector, color=color, scale=scale, name=name, show=False)
        if vec:
            vecs.append(vec)
    comp = Part.makeCompound(vecs)
    if show:
        show_object(comp, name=name, color=color)
    return comp

def draw_points(points, color=(1.0,0,0), size=2.0, name="Pt", show=True):    
    """
    Method to draw points as spheres in FreeCAD for visualization.

    Args:
        points: List of points (App.Vector)
        color: Sphere color (r,g,b)
        size: Sphere radius
        name: Name for the object
        show: Whether to create the object in FreeCAD document
    Returns:
        Part.Compound or None: The Compound of spheres
    """    
    compound = Part.makeCompound([Part.makeSphere(size, point) for point in points])

    if show:
        obj = App.ActiveDocument.addObject("Part::Feature", name)
        obj.Shape = compound
        obj.ViewObject.ShapeColor = color
        obj.ViewObject.LineColor = color
        obj.ViewObject.PointColor = color
    return compound
### End debug visualization utilities ###


#### Core functionality ###

class FairingDriver:
    def __init__(self, axis, side, max_dist, strength, decay):
        self.axis = axis
        self.side = side
        self.max_dist = max_dist
        self.strength = strength
        self.decay = decay

def insert_uniform_knots(bs, along_v, target_poles=8, tol=1e-9):
    """
    Insert uniform knots into a BSplineSurface to reach at least target pole count.
    Shape is preserved exactly.

    Args:
        bs (Part.BSplineSurface): surface to refine (MODIFIED IN PLACE)
        along_v (bool): True → boundary is U=const (refine V)
        target_poles (int): minimum number of poles in that direction
        tol (float): knot insertion tolerance
    Returns:
        bool: True if the surface was modified, False otherwise
    """
    modified = False

    if along_v:
        get_nb_poles = lambda: bs.NbVPoles
        get_knots    = lambda: bs.getVKnots()
        get_mults    = lambda: bs.getVMultiplicities()
        insert_knot  = lambda v: bs.insertVKnot(v, 1, tol)
        param_range  = bs.bounds()[2:4]
    else:
        get_nb_poles = lambda: bs.NbUPoles
        get_knots    = lambda: bs.getUKnots()
        get_mults    = lambda: bs.getUMultiplicities()
        insert_knot  = lambda u: bs.insertUKnot(u, 1, tol)
        param_range  = bs.bounds()[0:2]

    cur_poles = get_nb_poles()
    if cur_poles >= target_poles:
        return modified  # nothing to do

    umin, umax = param_range
    existing_knots = get_knots()
    existing_mults = get_mults()

    # ---- compute how many new poles we need ----
    needed = target_poles - cur_poles

    # ---- generate candidate uniform parameters ----
    total_segments = cur_poles + needed - 1
    step = (umax - umin) / total_segments

    candidates = [umin + i * step for i in range(1, total_segments)]

    # ---- filter out knots already present ----
    def is_existing(u):
        for k, m in zip(existing_knots, existing_mults):
            if abs(k - u) < tol and m >= 1:
                return True
        return False

    new_knots = [u for u in candidates if not is_existing(u)]

    # ---- insert knots until target reached ----
    for u in new_knots:
        insert_knot(u)
        modified = True
        if get_nb_poles() >= target_poles:
            break

    return modified

def insert_boundary_parallel_knots(bs, along_v, edge, layers, tol=1e-9):
    """
    Insert knots parallel to boundary direction to increase tangential DOF.
    """
    if layers <= 0:
        return False

    # Determine knot insertion direction
    if along_v:
        insert = lambda t: bs.insertUKnot(t, 1, tol)
        param_min, param_max = bs.bounds()[0:2]
    else:
        insert = lambda t: bs.insertVKnot(t, 1, tol)
        param_min, param_max = bs.bounds()[2:4]

    u0, u1 = edge.ParameterRange
    eps = (param_max - param_min) * 1e-6

    modified = False

    # Uniform subdivision of boundary parameter
    for k in range(1, layers + 1):
        s = k / (layers + 1)
        t = u0 + s * (u1 - u0)

        # Clamp safely
        t = max(param_min + eps, min(param_max - eps, t))
        insert(t)
        modified = True

    return modified

def insert_boundary_biased_knots(bs, along_v, boundary_index, layers=3, bias=3.0, tol=1e-9):
    """
    Insert knots clustered near a boundary to increase DOF
    without deforming the surface.

    Args:
        bs (Part.BSplineSurface)
        along_v (bool): True → boundary is U=const (refine V)
        boundary_index (int): 0 or max index of boundary row/column
        layers (int): number of layers to insert
        bias (float): >1 clusters knots near boundary
        tol (float) : knot insertion tolerance
    Returns:
        bool: True if the surface was modified, False otherwise
    """
    modified = False

    # --- Direction setup ---
    if along_v:
        # Boundary is U=const → refine V
        get_knots   = lambda: bs.getVKnots()
        get_mults   = lambda: bs.getVMultiplicities()
        insert_knot = lambda u: bs.insertVKnot(u, 1, tol)
        param_min, param_max = bs.bounds()[2:4]
    else:
        # Boundary is V=const → refine U
        get_knots   = lambda: bs.getUKnots()
        get_mults   = lambda: bs.getUMultiplicities()
        insert_knot = lambda u: bs.insertUKnot(u, 1, tol)
        param_min, param_max = bs.bounds()[0:2]
    

    # --- Boundary parameter ---
    at_start = (boundary_index == 0)
    boundary_param = param_min if at_start else param_max

    # --- Existing knots ---
    existing_knots = get_knots()
    existing_mults = get_mults()

    def has_knot(u):
        for k, m in zip(existing_knots, existing_mults):
            if abs(k - u) < tol and m >= 1:
                return True
        return False

    # --- Generate biased offsets ---
    span = param_max - param_min
    eps = span * 1e-6  # avoid exact boundary

    for i in range(1, layers + 1):
        s = (i / (layers + 1)) ** bias
        d = s * span * 0.5

        if at_start:
            u = boundary_param + d
        else:
            u = boundary_param - d

        # Clamp safely inside domain
        u = max(param_min + eps, min(param_max - eps, u))

        if not has_knot(u):
            insert_knot(u)
            modified = True

    return modified

def weighted_refine_params(params, weights, layers):
    inserts = []
    n = len(params)

    for layer in range(layers):
        acc = 0.0
        target = (layer + 1) / (layers + 1)

        for i in range(n - 1):
            acc += weights[i]
            if acc >= target:
                t = 0.5 * (params[i] + params[i + 1])
                inserts.append(t)
                break

    return inserts

def compute_mismatch_weights(driver_dirs, target_dirs, eps=1e-6):
    weights = []
    for td, ts in zip(driver_dirs, target_dirs):
        w = App.Vector(td).sub(ts).Length
        weights.append(max(w, eps))
    return weights

def normalize_weights(w):
    s = sum(w)
    if s < 1e-12:
        return [1.0 / len(w)] * len(w)
    return [x / s for x in w]

def insert_mismatch_weighted_knots(bs, along_v, boundary_index, edge,
        driver_dirs, target_dirs, layers=3, tol=1e-9):
    """
    Insert transverse knots near boundary weighted by tangent mismatch.

    Args:
        bs (Part.BSplineSurface): target surface
        along_v (bool): True → boundary is U=const (refine V)
        boundary_index (int): 0 or max index of boundary row/column
        edge (Part.Edge): edge corresponding to boundary
        driver_dirs (list of App.Vector): driver transverse tangents
        target_dirs (list of App.Vector): target transverse tangents
        layers (int): number of layers to insert
        tol (float) : knot insertion tolerance
    Returns:
        bool: True if the surface was modified, False otherwise
    """

    if not driver_dirs or not target_dirs:
        return False

    samples = min(len(driver_dirs), len(target_dirs))
    driver_dirs = driver_dirs[:samples]
    target_dirs = target_dirs[:samples]

    # Direction setup
    if along_v:
        insert = lambda t: bs.insertVKnot(t, 1, tol)
        param_min, param_max = bs.bounds()[2:4]
    else:
        insert = lambda t: bs.insertUKnot(t, 1, tol)
        param_min, param_max = bs.bounds()[0:2]

    u0, u1 = edge.ParameterRange
    edge_params = [u0 + (u1 - u0) * i / (samples - 1)
                   for i in range(samples)]

    # Compute mismatch weights
    weights = compute_mismatch_weights(driver_dirs, target_dirs)
    weights = normalize_weights(weights)

    # Generate insertion parameters
    refine_params = weighted_refine_params(
        edge_params,
        weights,
        layers
    )

    # Map edge → surface parameter
    eps = (param_max - param_min) * 1e-6
    at_start = (boundary_index == 0)
    boundary_param = param_min if at_start else param_max

    modified = False
    for t in refine_params:
        s = (t - u0) / (u1 - u0)
        d = s * (param_max - param_min) * 0.5

        u = boundary_param + d if at_start else boundary_param - d
        u = max(param_min + eps, min(param_max - eps, u))

        insert(u)
        modified = True

    return modified

def estimate_max_tangent_angle(driver_dirs, target_dirs):
    """
    Estimate maximum angular mismatch between driver and target boundary tangents.

    Args:
        driver_dirs: list[App.Vector or None] - normalized directions
        target_dirs: list[App.Vector or None] - normalized directions

    Returns:
        float - maximum angle in radians
    """
    if not driver_dirs or not target_dirs:
        return 0.0

    max_angle = 0.0

    for d,t in zip(driver_dirs, target_dirs):

        if d is None or t is None:
            continue

        if d.Length < 1e-9 or t.Length < 1e-9:
            continue

        # Ensure unit length (cheap safety)
        d = App.Vector(d)
        t = App.Vector(t)
        d.normalize()
        t.normalize()

        # Ignore sign (tangent orientation)
        dot = abs(d.dot(t))

        # Numerical clamp
        dot = max(-1.0, min(1.0, dot))

        angle = math.acos(dot)
        max_angle = max(angle, max_angle)

    return max_angle

def estimate_boundary_angular_capacity(bs, along_v, boundary_index):
    """
    Estimate how much tangent angle the surface can represent along the boundary
    without adding parallel knots.
    Args:
        bs (Part.BSplineSurface): target surface
        along_v (bool): True → boundary is V=const (measure along U)
        boundary_index (int): index of boundary row/column
    Returns:
        float - maximum angle in radians
    """
    if along_v:
        # boundary is U = const → measure along U
        row = [bs.getPole(i+1, boundary_index + 1) for i in range(bs.NbUPoles)]        
    else:
        # boundary is V = const → measure along V
        row = [bs.getPole(boundary_index + 1, j+1) for j in range(bs.NbVPoles)]

    angles = []
    for i in range(len(row) - 2):
        v1 = row[i+1].sub(row[i])
        v2 = row[i+2].sub(row[i+1])
        if v1.Length < 1e-9 or v2.Length < 1e-9:
            continue
        v1.normalize()
        v2.normalize()
        angles.append(math.acos(max(-1.0, min(1.0, v1.dot(v2)))))

    return max(angles) if angles else 0.0

def estimate_parallel_refinement_need(
    sample_params,
    angle_samples_rad,
    hot_angle_rad=0.17,   # ~10 degrees
    min_ratio=2.0
):
    """
    Decide whether localized parallel refinement is needed.

    Args:
        sample_params: list[float] – boundary params (monotone)
        angle_samples_rad: list[float] – tangent mismatch angles
        hot_angle_rad: float – threshold for "hot" region
        min_ratio: float – stiffness threshold

    Returns:
        dict or None:
            {
                "u0": float,
                "u1": float,
                "ratio": float,
                "severity": float
            }
    """
    assert len(sample_params) == len(angle_samples_rad)

    # Detect hot region
    hot_idx = [i for i,a in enumerate(angle_samples_rad) if a >= hot_angle_rad]
    if not hot_idx:
        return None

    i0, i1 = hot_idx[0], hot_idx[-1]
    u0, u1 = sample_params[i0], sample_params[i1]
    hot_width = max(u1 - u0, 1e-12)

    # Effective tangential stiffness proxy:
    # assume first interior row "sees" one knot span
    # we normalize by boundary param length = 1.0
    stiffness_ratio = 1.0 / hot_width

    if stiffness_ratio < min_ratio:
        return None

    severity = max(angle_samples_rad[i0:i1+1]) / hot_angle_rad

    return {
        "u0": u0,
        "u1": u1,
        "ratio": stiffness_ratio,
        "severity": severity
    }

def estimate_parallel_layers_from_region(region, max_layers=8):
    """
    Estimate number of parallel knots needed from region severity.
    """
    if region is None:
        return 0

    # severity ≈ how many "curvature quanta"
    layers = int(math.ceil(region["severity"]))

    return min(max(1, layers), max_layers)

def estimate_layers_from_region_width(width, min_trans, min_parallel):
    if width <= 1e-6:
        return min_trans, min_parallel

    transverse = max(min_trans, int(math.ceil(min_trans / width)))
    parallel   = max(min_parallel, int(math.ceil(min_parallel / width)))

    return transverse, parallel

def insert_localized_parallel_knots(bs, along_v, region, layers=8, tol=1e-9):
    """
    Insert parallel knots localized to hot curvature region.

    Args:
        bs: Part.BSplineSurface
        along_v: bool - True if boundary is V=const (so refine V)
        region: dict from estimate_parallel_refinement_need()
        layers: int - number of layers to insert
        tol: float - knot insertion tolerance
    Returns:
        bool - True if surface was modified
    """
    if region is None or layers <= 0:
        return False
    
    u0, u1 = region["u0"], region["u1"]
    if u1 <= u0:
        return False

    # Determine knot insertion direction
    if not along_v:
        insert = lambda t: bs.insertUKnot(t, 1, tol)
    else:
        insert = lambda t: bs.insertVKnot(t, 1, tol)

    # Gaussian-like distribution around center
    mid = 0.5 * (u0 + u1)
    span = 0.5 * (u1 - u0)

    modified = False
    for i in range(layers):
        # Bias toward center
        w = (i + 1) / (layers + 1)
        t = mid + span * (2*w - 1) * 0.6

        insert(t)
        modified = True

    return modified

def is_straight_line(edge, tol=1e-7):
    """
    Check if an edge is (almost) a straight line.
    Compare edge length with distance between its endpoints.
    Args:
        edge: Part.Edge
        tol: Tolerance for straightness
    Returns:
        True if edge is straight within tolerance
    """
    length = edge.Length
    start_pt = edge.Vertexes[0].Point
    end_pt = edge.Vertexes[-1].Point
    chord = (end_pt.sub(start_pt)).Length
    return abs(length - chord) <= tol

def points_are_close(p1, p2, tol=1e-7):
    """
    Check if two points are close within tolerance.
    
    Args:
        p1: First point (App.Vector)
        p2: Second point (App.Vector)
        tol: Tolerance for closeness
    """
    return (p1.sub(p2)).Length <= tol

def edges_are_close(e1, e2, tol=1e-4, samples=50):
    """
    Check if two edges are close by sampling points along their curves.
    
    Args:
        e1: First edge (Part.Edge)
        e2: Second edge (Part.Edge)
        tol: Tolerance for closeness
        samples: Number of sample points along each edge
    Returns:
        True if edges are close within tolerance
    """
    pts1 = e1.discretize(samples)
    pts2 = e2.discretize(samples)

    if pts1[0].sub(pts2[0]).Length > tol:
        pts2 = list(reversed(pts2))

    for pt1, pt2 in zip(pts1, pts2):
        if (pt1.sub(pt2)).Length > tol:
            return False
    return True

def find_shared_edge(face1, face2, tol=1e-3):
    """
    Detect shared edge between two faces.
    
    Args:
        face1: First face
        face2: Second face
        tol: Tolerance for edge comparison (detailed sampling)
    Returns:
        Shared edge (Part.Edge) or None
    """
    # Quick vertex-based filter before expensive sampling
    for e1 in face1.Edges:
        v1_start, v1_end = e1.Vertexes[0].Point, e1.Vertexes[-1].Point
        for e2 in face2.Edges:
            v2_start, v2_end = e2.Vertexes[0].Point, e2.Vertexes[-1].Point

            # Check if vertices roughly match (either orientation)
            if (points_are_close(v1_start, v2_start, tol=1e-7) and points_are_close(v1_end, v2_end, tol=1e-7)) or \
               (points_are_close(v1_start, v2_end, tol=1e-7) and points_are_close(v1_end, v2_start, tol=1e-7)):

                # Quick straight-line check
                if is_straight_line(e1, tol=1e-5) and is_straight_line(e2, tol=1e-5):
                    return e1
                
                # Then check detailed sampling closeness
                if edges_are_close(e1, e2, tol, samples=50):
                    return e1

    return None

def lerp(v1, v2, t):
    """
    Linear Interpolation
    
    Args:
        v1: First vector
        v2: Second vector
        t: Interpolation factor (0.0 - 1.0)
    Returns:
        Interpolated unit vector
    """
    a = App.Vector(v1)
    b = App.Vector(v2)

    if a.Length < 1e-9:
        return b
    if b.Length < 1e-9:
        return a

    v = a.multiply(1 - t).add(b.multiply(t))
    if v.Length < 1e-9:
        return a

    v.normalize()
    return v
              
def v_avg(vecs):
        v = App.Vector(0,0,0)
        for vec in vecs:
            v = v.add(vec)
        return v.multiply(1.0 / len(vecs))
    
def w_avg(ws):
    d = 0.0
    for w in ws:
        d += w
    return d / len(ws)

def isRational(bs):
    """
    Check if BSpline surface is rational in either U or V direction.

    Args:
        bs: BSpline surface
    Returns:
        bool: True if rational in U or V
    """
    return bs.isURational() or bs.isVRational()

def detect_edge_direction_and_boundary(bs, edge):
    """
    Detect if edge runs along U or V direction of the BSpline surface,
    and determine which boundary row/column it corresponds to.

    Args:
        bs: BSpline surface
        edge: Edge to analyze
    Returns:
        along_v: True if edge runs along V (U=const), False if along U (V=const)
        boundary_index: index of boundary row/column (0-based)
    """
    us = []
    vs = []

    u0, u1 = edge.ParameterRange

    for i in range(5):
        t = u0 + (u1 - u0) * i / 4.0
        p = edge.Curve.value(t)
        u, v = bs.parameter(p)
        us.append(u)
        vs.append(v)

    du = max(us) - min(us)
    dv = max(vs) - min(vs)

    u_min, u_max, v_min, v_max = bs.bounds()
    nu = bs.NbUPoles
    nv = bs.NbVPoles

    if du < dv:
        # Edge runs along V, boundary is U = const
        along_v = True

        if abs(us[0] - u_min) < abs(us[0] - u_max):
            boundary_index = 0      # first U row
        else:
            boundary_index = nu - 1     # last U row
    else:
        # Edge runs along U, boundary is V = const
        along_v = False

        if abs(vs[0] - v_min) < abs(vs[0] - v_max):
            boundary_index = 0      # first V row
        else:
            boundary_index = nv - 1     # last V row

    return along_v, boundary_index

def transverse_tangent_iso(surface, point, edge_tangent):
    """
    Compute surface tangent transverse to the edge using isoparametric curves.

    Args:
        surface: Part.Surface
        point: App.Vector
        edge_tangent: App.Vector (unit)
    Returns:
        transverse tangent App.Vector (unit) or None
    """

    u, v = surface.parameter(point)

    Tu = None
    Tv = None

    # --- U isoparametric tangent
    try:
        uiso = surface.uIso(u)
        Tu = uiso.tangent(v)[0]
    except:
        pass

    # --- V isoparametric tangent
    try:
        viso = surface.vIso(v)
        Tv = viso.tangent(u)[0]
    except:
        pass

    candidates = []

    if Tu:
        Tu_proj = Tu.sub(App.Vector(edge_tangent).multiply(Tu.dot(edge_tangent)))
        if Tu_proj.Length > 1e-9:
            candidates.append(Tu_proj)

    if Tv:
        Tv_proj = Tv.sub(App.Vector(edge_tangent).multiply(Tv.dot(edge_tangent)))
        if Tv_proj.Length > 1e-9:
            candidates.append(Tv_proj)

    if not candidates:
        return None

    # pick strongest transverse direction
    Tt = max(candidates, key=lambda x: x.Length)
    Tt.normalize()
    return Tt

def transverse_tangent_derivative(surface, point, edge_tangent):
    """
    Compute surface tangent transverse to the edge using surface normal.
        
    Args:
        surface: Part.Surface
        point: App.Vector
        edge_tangent: App.Vector (unit)
    Returns:
        transverse tangent App.Vector (unit) or None
    """
    try:
        u, v = surface.parameter(point)
        du = surface.getDN(u, v, 1, 0)
        dv = surface.getDN(u, v, 0, 1)
    except:
        return None

    # Surface normal
    N = du.cross(dv)
    if N.Length < 1e-9:
        return None
    N.normalize()

    # Transverse direction
    Tt = N.cross(edge_tangent)
    if Tt.Length < 1e-9:
        return None

    Tt.normalize()
    return Tt

def transverse_derivative_basis(surface, point, edge_tangent):
    """
    Return orthonormal basis vectors spanning the surface's transverse
    derivative space at this boundary point.
    Args:
        surface: (Part.Surface) driver surface
        point: (App.Vector) point on surface
        edge_tangent: (App.Vector) unit vector tangent to edge
    """
    u, v = surface.parameter(point)

    try:
        du = surface.getDN(u, v, 1, 0)
        dv = surface.getDN(u, v, 0, 1)
    except:
        return []

    basis = []

    for d in (du, dv):    
        # Remove edge tangent component (Gram–Schmidt)    
        d = d.sub(App.Vector(edge_tangent).multiply(d.dot(edge_tangent)))        
        if d.Length > 1e-9:
            d.normalize()
            basis.append(d)

    # Orthonormalize (Gram–Schmidt)
    ortho = []
    for b in basis:
        v = App.Vector(b)
        for o in ortho:
            v = v.sub(App.Vector(o).multiply(v.dot(o)))
        if v.Length > 1e-9:
            v.normalize()
            ortho.append(v)

    return ortho

def best_possible_transverse(driver_face, point, edge_tangent, driver_dir):
    """
    Returns the best achievable transverse direction for G1 continuity.
    Args:
        driver_face: Part.Face
        point: App.Vector
        edge_tangent: App.Vector (unit)
        driver_dir: App.Vector (unit)
    Returns:
        transverse tangent App.Vector (unit) or None
    """    
    surface = driver_face.Surface
    
    # Try transverse derivative basis projection
    basis = transverse_derivative_basis(surface, point, edge_tangent)
    
    if basis:
        Td = App.Vector(0, 0, 0)
        for b in basis:
            Td = Td.add(b.multiply(driver_dir.dot(b)))
        
        if Td.Length > 1e-9:
            Td.normalize()
            return Td
    else:
        App.Console.PrintWarning("No transverse derivative basis found\n")
        
    # Fallback: try isoparametric curves
    Td = transverse_tangent_iso(surface, point, edge_tangent)
    if Td:
        Td.normalize()
        return Td
    else:
        App.Console.PrintWarning("Failed to compute transverse tangent using isoparametric curves\n")

    # Fallback: try surface derivative
    Td = transverse_tangent_derivative(surface, point, edge_tangent)
    if Td:
        Td.normalize()
        return Td
    else:
        App.Console.PrintWarning("Failed to compute transverse tangent using surface derivative\n")
    
    # Final fallback: return original driver_dir
    return driver_dir

def find_transverse_tangent(driver_face, point, edge_tangent):
    """    
    Compute surface tangent transverse to the edge using projection onto basis.

    Args:
        driver_face: Part.Face
        point: App.Vector
        edge_tangent: App.Vector (unit)

    Returns:
        transverse tangent App.Vector (unit) or None
    """
    u, v = driver_face.Surface.parameter(point)
    N = driver_face.normalAt(u, v)

    # Guard against zero normal (rare but possible near degeneracies)
    if N.Length < 1e-9:
        App.Console.PrintWarning(f"Zero normal vector at point {point}\n")
        return None
    
    N.normalize()

    driver_dir = App.Vector(edge_tangent).cross(N)
    if driver_dir.Length < 1e-9:
        App.Console.PrintWarning(f"Zero driver direction vector at point {point}\n")
        return None

    driver_dir.normalize()

    d_dir = best_possible_transverse(driver_face, point, edge_tangent, driver_dir)

    return d_dir

def boundary_pole_params(bs, along_v, boundary_index):
    """
    Return surface parameters corresponding to boundary poles.

    Args:
        bs (Part.BSplineSurface): target surface
        along_v (bool):
            True  → boundary is U = const, poles vary along V
            False → boundary is V = const, poles vary along U
        boundary_index (int):
            0 → first boundary
            >0 → last boundary

    Returns:
        list[float]: parameter values along boundary direction
                     (U or V depending on along_v)
    """

    nu = bs.NbUPoles
    nv = bs.NbVPoles

    params = []

    if along_v:
        # U = const boundary, iterate along V
        i = 1 if boundary_index == 0 else nu
        for j in range(1, nv + 1):
            P = bs.getPole(i, j)
            u, v = bs.parameter(P)
            params.append(v)
    else:
        # V = const boundary, iterate along U
        j = 1 if boundary_index == 0 else nv
        for i in range(1, nu + 1):
            P = bs.getPole(i, j)
            u, v = bs.parameter(P)
            params.append(u)

    return params

def boundary_sample_params(bs, along_v, samples):
    """
    Sample surface parameters along boundary edge.

    Args:
        bs (Part.BSplineSurface): target surface
        along_v (bool): True if boundary/sampling direction is along V, False if along U
        samples (int): number of samples along edge
    Returns:
        list of float: sampled surface parameters along boundary direction
    """
    params = []

    if along_v:
        # U = const boundary, iterate along V
        v_min, v_max = bs.bounds()[2:4]

        for k in range(samples):
            s = k / (samples - 1)
            v = v_min + s * (v_max - v_min)
            params.append(v)
    else:
        # V = const boundary, iterate along U
        u_min, u_max = bs.bounds()[0:2]

        for k in range(samples):
            s = k / (samples - 1)
            u = u_min + s * (u_max - u_min)
            params.append(u)

    return params

def get_surface_boundary_tangents(face, edge, samples = 30, pole_params = None):
    """
    Get transverse tangents at boundary poles.

    Args:
        face: Part.Face
        edge: Part.Edge
        samples: int    (number of samples along edge)    
        pole_params: list of float (surface parameters u or v of boundary poles)
    Returns:
        list of App.Vector (unit) per boundary pole
    """
    u0, u1 = edge.ParameterRange
    eps = (u1 - u0) * 1e-9

    # Sample driver tangents along edge
    sample_params = []
    sample_dirs = []

    for i in range(samples):
        if i == 0:
            t = u0 + eps
        elif i == samples - 1:
            t = u1 - eps
        else:
            t = u0 + (u1 - u0) * i / (samples - 1)

        p = edge.Curve.value(t)

        # Edge tangent
        Te = edge.Curve.tangent(t)[0]
        if Te.Length < 1e-9:
            continue
        Te.normalize()

        # Transverse tangent on surface
        Td = find_transverse_tangent(face, p, Te)
        if Td is None or Td.Length < 1e-9:
            continue

        Td.normalize()

        sample_params.append(t)
        sample_dirs.append(Td)

    if len(sample_dirs) < 2:
        raise RuntimeError("Insufficient valid driver tangents")

    # Ensure directional consistency
    for i in range(1, len(sample_dirs)):
        if sample_dirs[i - 1].dot(sample_dirs[i]) < 0:
            sample_dirs[i] = sample_dirs[i].negative()

    if not pole_params:
        return sample_dirs

    # Map pole parameters to edge parameters
    # Normalize pole params to edge parameter space
    pole_dirs = []
    pole_locs = []
    for p in pole_params:
        # Clamp to sample params range
        pe = max(min(sample_params), min(max(sample_params), p))

        # Find bracketing samples
        for k in range(len(sample_params) - 1):
            t0 = sample_params[k]
            t1 = sample_params[k + 1]

            if pe >= t0 and pe <= t1:                
                w = (pe - t0) / (t1 - t0)
                d0 = App.Vector(sample_dirs[k])
                d1 = App.Vector(sample_dirs[k + 1])

                # Ensure same orientation - redundant but safer
                if d0.dot(d1) < 0:
                    d1 = d1.negative()

                Td = d0.multiply(1 - w).add(d1.multiply(w))
                if Td.Length < 1e-9:
                    Td = App.Vector(sample_dirs[k])

                Td.normalize()
                pole_locs.append(edge.Curve.value(pe))
                pole_dirs.append(Td)
                break
        else:
            # Outside range (numerical edge case)
            if pe < sample_params[0]:
                pole_dirs.append(App.Vector(sample_dirs[0]))
            else:
                pole_dirs.append(App.Vector(sample_dirs[-1]))
                
    # draw_vectors(
    #     [edge.Curve.value(t) for t in sample_params],
    #     sample_dirs,
    #     color=(1.0,0.5,0.0),
    #     scale=10.0,
    #     name="Driver_Tangents_Orange",
    #     show=True
    # )
    # draw_vectors(
    #     pole_locs,
    #     pole_dirs,
    #     color=(0.0,0.0,1.0),
    #     scale=10.0,
    #     name="Pole_Tangents_Blue",
    #     show=True
    # )
    return pole_dirs

def compute_g1_row(orig_poles, driver_tangents_dir, along_v, boundary_index, beta):
    """
    Compute new pole row enforcing G1 continuity.

    Args:
        orig_poles: 2D list [nu][nv] of App.Vector
        driver_tangents_dir: list of App.Vector or None (per pole)
        along_v: True if boundary is U=const, False if V=const
        boundary_index: index of boundary row/column (0-based)
        beta: magnitude blending factor [0..1]
    Returns:
        g1_row: list of App.Vector (new adjacent pole row)
    """
    g1_row = []

    row_len = len(driver_tangents_dir)
    step = 1 if boundary_index == 0 else -1
    adj_row = boundary_index + step
    
    for j in range(row_len):
        if along_v:
            P1 = orig_poles[boundary_index][j]
            P2 = orig_poles[adj_row][j]
        else:
            P1 = orig_poles[j][boundary_index]
            P2 = orig_poles[j][adj_row]

        T = driver_tangents_dir[j]
        if T is None or T.Length < 1e-9:
            g1_row.append(P2)
            continue

        Dorig = P2.sub(P1)
        d_orig = Dorig.Length
        if d_orig < 1e-9:
            g1_row.append(P2)
            continue

        # Ensure same orientation
        if Dorig.dot(T) < 0:
            T = T.negative()
            
        Norig = App.Vector(Dorig)
        Norig.normalize()
        
        # Direction blending
        Tblend = Norig.multiply(1.0 - beta).add(T.multiply(beta))
        if Tblend.Length < 1e-6:
            Tblend = T

        Tblend.normalize()

        Tfinal = Tblend.multiply(d_orig)
        # draw_vector(P1, Tfinal, (0.0,0.0,1.0), 10.0, "Tfinal_Blue")

        g1_row.append(P1.add(Tfinal))

    return g1_row

def get_spread_rows(g1_idx, spread_rows, step):
    """
    Find interior rows to apply G1 spreading.

    Args:
        g1_idx: index of boundary g1 row (0-based)
        spread_rows: requested number of interior rows
        step: direction step (1 or -1)

    Returns:
        list of row indices (0-based), ordered from edge inward
    """

    r1 = g1_idx + step
    rows = []
    for _ in range(spread_rows):
        rows.append(r1)
        r1 += step        

    return rows

def get_deweighted_pole(bs, i, j):
    """
    Get de-weighted pole from BSpline surface.
    Args:
        bs: (Part.BSplineSurface) target rational surface
        i: pole index in U (0-based)
        j: pole index in V (0-based)
    """    
    P = bs.getPole(i+1, j+1)
    w = bs.getWeight(i+1, j+1)
    if abs(w) < 1e-12:
        return P
    return P.multiply(1.0 / w)

def pole_spacing(bs, along_v, boundary_index):
    """
    Gets pole spacing near boundary
    
    Args:
    bs: (Part.BSplineSurface) target surface
    along_v: direction of boundary
    boundary_index: index of boundary row/column (0-based)
    """
    spacings = []
    i0 = boundary_index
    i1 = boundary_index + (1 if boundary_index == 0 else -1)   
    if along_v:
        for j in range(bs.NbVPoles):
            P0 = get_deweighted_pole(bs, i0, j)
            P1 = get_deweighted_pole(bs, i1, j)
            spacings.append(P1.sub(P0).Length)
    else:
        for i in range(bs.NbUPoles):
            P0 = get_deweighted_pole(bs, i, i0)
            P1 = get_deweighted_pole(bs, i, i1)
            spacings.append(P1.sub(P0).Length)

    return sum(spacings) / max(len(spacings), 1)

def estimate_refinement_layers(driver_dirs, target_dirs, beta, degree, energy_tol=0.15, max_layers=8):
    """
    Estimate transverse refinement layers from predicted G1 energy
    Args:
        driver_dirs: list of App.Vector (unit) from driver surface
        target_dirs: list of App.Vector (unit) from target surface
        beta: blending factor
        degree: surface degree in refinement direction
        energy_tol: tolerance for acceptable energy
        max_layers: maximum allowed refinement layers
    """

    if degree <= 0 or not driver_dirs:
        return 0

    max_mismatch = 0.0

    for td, ts in zip(driver_dirs, target_dirs):
        d = App.Vector(td).sub(ts)
        max_mismatch = max(max_mismatch, d.Length)

    # predicted dimensionless energy
    E = beta * max_mismatch / degree

    if E < energy_tol:
        return 0
    
    if E > 6.0:
        App.Console.PrintWarning(f"High predicted G1 energy={E:.2f}. G1 may not be possible.\n")
        return max_layers

    layers = int(math.ceil(math.sqrt(E / energy_tol)))
    return max(0, min(layers, max_layers))

def resolve_constraints(desired_dirs):
    """
    Build per-pole motion constraints from desired directions.

    Args:
        desired_dirs: dict[(i,j)] -> list of entries:
            {
                "dir": App.Vector (unit),
                "mag": float,
                "along_v": bool,
                "boundary": int,
                "driver_idx": int,
                "priority": int
            }

    Returns:
        pole_constraints: dict[(i,j)] -> constraint dict
    """

    pole_constraints = {}

    for key, entries in desired_dirs.items():

        if not entries:
            continue

        # Priority filtering
        hard = [e for e in entries if e["priority"] == 0]

        if hard:
            active = hard
        else:
            active = entries

        if not active:
            continue


        # sanitize
        dirs = []
        mags = []
        priorities = []

        for e in active:
            d = App.Vector(e["dir"])
            if d.Length < 1e-9:
                continue
            d.normalize()
            dirs.append(d)
            mags.append(abs(e["mag"]))
            priorities.append(e["priority"])

        if not dirs:
            continue

        # Single driver
        if len(dirs) == 1:
            pole_constraints[key] = {
                "mode": "single",
                "dir": dirs[0],
                "max_mag": mags[0],
                "priority": priorities[0]
            }
            continue

        # Cluster directions by alignment
        clusters = []  # list of {"dir": vec, "mags": [], "dirs": [], "priorities": []}

        for d, m, p in zip(dirs, mags, priorities):
            placed = False
            for c in clusters:
                if abs(c["dir"].dot(d)) > 1.0 - 1e-6:
                    c["mags"].append(m)
                    c["dirs"].append(d)
                    c["priorities"].append(p)
                    placed = True
                    break
            if not placed:
                clusters.append({
                    "dir": App.Vector(d),
                    "dirs": [d],
                    "mags": [m],
                    "priorities": [p]
                })

        def cluster_max_mag_priority(c):
            max_mag = max(c["mags"])            
            for m, p in zip(c["mags"], c["priorities"]):
                if m == max_mag:
                    return p
            return min(c["priorities"])
        
        # --- one dominant cluster ---
        if len(clusters) == 1:
            c = clusters[0]
            pole_constraints[key] = {
                "mode": "single",
                "dir": c["dir"],
                "max_mag": max(c["mags"]),
                "priority": min(c["priorities"])
            }
            continue

        # --- try to form an orthogonal plane (corner case) ---
        best_pair = None
        best_dot = 1.0

        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d1 = clusters[i]["dir"]
                d2 = clusters[j]["dir"]
                dot = abs(d1.dot(d2))
                if dot < best_dot:
                    best_dot = dot
                    best_pair = (clusters[i], clusters[j])

        # sufficiently orthogonal?
        if best_pair and best_dot < 0.3:  # ~72 degrees
            c1, c2 = best_pair

            # Orthonormalize
            u = App.Vector(c1["dir"])
            u.normalize()
            v = App.Vector(c2["dir"])
            v = v.sub(App.Vector(u).multiply(v.dot(u)))
            if v.Length > 1e-9:
                v.normalize()

                pole_constraints[key] = {
                    "mode": "plane",
                    "basis": [u, v],
                    "max_mag": min(max(c1["mags"]), max(c2["mags"])),
                    "priority": min(min(c1["priorities"]), min(c2["priorities"]))
                }
                continue

        # --- fallback: clamp ---
        # Multiple conflicting directions, low DOF
        # Allow movement but cap magnitude and priority
        avg_dir = v_avg(dirs)

        if avg_dir.Length > 1e-9:
            avg_dir.normalize()

        pole_constraints[key] = {
            "mode": "clamp",
            "dir": avg_dir,
            "max_mag": min(mags),
            "priority": min(priorities)
        }

    return pole_constraints

def effective_max_mag(c):
    if c["priority"] <= 0:
        return c["max_mag"]
    return c["max_mag"] * (1.0 / (0.3 + c["priority"]))

def apply_constraint(delta, c):
    if delta.Length < 1e-9:
        return delta

    max_mag = effective_max_mag(c)

    if c["mode"] == "single":
        d = App.Vector(c["dir"])
        mag = delta.dot(d)
        mag = max(-max_mag, min(max_mag, mag))
        return d.multiply(mag)

    if c["mode"] == "plane":
        proj = App.Vector(0, 0, 0)
        for b in c["basis"]:
            proj = proj.add(App.Vector(b).multiply(delta.dot(b)))

        if proj.Length > max_mag:
            proj.normalize()
            proj = proj.multiply(max_mag)
        return proj

    if c["mode"] == "clamp":
        if delta.Length > max_mag:
            delta = App.Vector(delta)
            delta.normalize()
            delta = delta.multiply(max_mag)
        return delta

    return delta

def build_pole_basis(bs):
    nu, nv = bs.NbUPoles, bs.NbVPoles
    basis = [[None for _ in range(nv)] for _ in range(nu)]

    poles = [[bs.getPole(i+1, j+1) for j in range(nv)] for i in range(nu)]

    for i in range(nu):
        for j in range(nv):
            b = []

            # U direction
            if 0 < i < nu - 1:
                du = poles[i+1][j].sub(poles[i-1][j])
                if du.Length > 1e-9:
                    du.normalize()
                    b.append(du)

            # V direction
            if 0 < j < nv - 1:
                dv = poles[i][j+1].sub(poles[i][j-1])
                if dv.Length > 1e-9:
                    dv.normalize()
                    b.append(dv)

            # Orthonormalize
            ortho = []
            for v in b:
                w = App.Vector(v)
                for o in ortho:
                    w = w.sub(o.multiply(w.dot(o)))
                if w.Length > 1e-9:
                    w.normalize()
                    ortho.append(w)

            basis[i][j] = ortho  # [] | [u] | [u,v]

    return basis

def linear_decay(d, max_d):
    return max(0.0, 1.0 - d / max_d)

def smooth_decay(d, max_d):
    t = min(1.0, d / max_d)
    return 1.0 - (3*t*t - 2*t*t*t)

def gaussian_decay(d, sigma):
    return math.exp(-(d*d) / (2*sigma*sigma))

def needs_degree_elevation(bs, along_v, boundary_index, transverse_layers, parallel_layers):
    """
    Check if current degrees are insufficient for requested refinement.
    Args:
        bs: (Part.BSplineSurface) target surface
        along_v: (bool) direction of boundary
        boundary_index: (int) index of boundary row/column (0-based)
        transverse_layers: (int) requested transverse refinement layers
        parallel_layers: (int) requested parallel refinement layers
    Returns:
        (uDeg, vDeg): tuple of bool indicating if degree elevation is needed in U and V
    """
    udeg, vdeg = bs.UDegree, bs.VDegree

    # refinement directions
    trans_deg = vdeg if not along_v else udeg
    para_deg  = udeg if not along_v else vdeg

    (uDeg, vDeg) = (False, False)
    
    # heuristic safety rule
    if transverse_layers > 2 * trans_deg:
        vDeg |= True if not along_v else False
        uDeg |= True if along_v else False
    if parallel_layers > 2 * para_deg:
        uDeg |= True if not along_v else False
        vDeg |= True if along_v else False
        
    return (uDeg, vDeg)

def count_poles_in_hot_region(bs, along_v, hot_region):
    """
    Count existing poles supporting curvature inside hot region.
    Returns:
        (parallel_poles, transverse_rows)
    """
    u0, u1 = hot_region["u0"], hot_region["u1"]

    u_knots = bs.getUKnots()
    v_knots = bs.getVKnots()

    if along_v:
        # boundary is V=const → parallel is U
        parallel_knots = u_knots
        transverse_knots = v_knots
    else:
        parallel_knots = v_knots
        transverse_knots = u_knots

    # Count knot spans overlapping hot region
    parallel_spans = sum(
        1 for i in range(len(parallel_knots) - 1)
        if not (parallel_knots[i+1] <= u0 or parallel_knots[i] >= u1)
    )

    transverse_spans = len(transverse_knots) - 1

    return parallel_spans, transverse_spans

def relax_layers_by_existing_dof(
    parallel_layers,
    transverse_layers,
    existing_parallel,
    existing_transverse,
    degree_parallel,
    degree_transverse
):
    """
    Reduce requested layers based on existing DOF.
    """

    parallel_capacity = existing_parallel * degree_parallel
    transverse_capacity = existing_transverse * degree_transverse

    parallel_layers = max(0, parallel_layers - parallel_capacity // 2)
    transverse_layers = max(0, transverse_layers - transverse_capacity // 2)

    return transverse_layers, parallel_layers

def estimate_surface_refinement(bs, driver_dirs, target_dirs, along_v, boundary_index, beta=1.0, energy_tol=0.15, hot_angle_rad=0.2):
    """
    Estimates required transverse + parallel refinement layers and degree elevation
    to achieve G1 continuity along shared edge between driver_face and target_face.
    Args:
        bs: (Part.BSplineSurface) target surface
        driver_dirs: (list of App.Vector) tangent directions on driver face
        target_dirs: (list of App.Vector) tangent directions on target face
        along_v: (bool) True if boundary is V=const, False if U=const
        boundary_index: (int) index of boundary row/column (0-based)
        energy_tol: (float) tolerance for acceptable G1 energy
        hot_angle_rad: (float) angle threshold for hot region detection (radians)
    Returns:
        dict with keys:
            "transverse_layers": int
            "parallel_layers": int
            "hot_region": dict or None
            "elevate_degree": (uDeg, vDeg) tuple of bool
    """

    # Transverse refinement from tangent mismatch
    degree = bs.VDegree if not along_v else bs.UDegree

    transverse_layers = estimate_refinement_layers(driver_dirs, target_dirs, beta, degree, energy_tol)
    parallel_layers = 0

    if transverse_layers <= 0:
        return {
            "transverse_layers": 0,
            "parallel_layers": parallel_layers,
            "hot_region": None,
            "elevate_degree": (False, False)
        }

    # Angular differences
    angles = [
        math.acos(max(-1.0, min(1.0, abs(d.dot(t)))))
        for d, t in zip(driver_dirs, target_dirs)
    ]

    # Hot region detection
    samples = len(angles)
    params = boundary_sample_params(bs, along_v, samples)

    hot_region = estimate_parallel_refinement_need(params, angles, hot_angle_rad)

    if hot_region:
        parallel_layers = estimate_parallel_layers_from_region(hot_region)

        print(f"Pre-estimated {transverse_layers} transverse layers and {parallel_layers} parallel layers")
        print(f"Detected hot region between {hot_region['u0']:.4f} and {hot_region['u1']:.4f}"
              f" with severity {hot_region['severity']:.4f} and stiffness {hot_region['ratio']:.4f}")
        
        # Coupled correction
        hot_width = hot_region["u1"] - hot_region["u0"]
        transverse_layers, parallel_layers = estimate_layers_from_region_width(hot_width, transverse_layers, parallel_layers)

        print(f"Region width corrected {transverse_layers} transverse layers and {parallel_layers} parallel layers")

        # relax layers based on existing DOF
        existing_parallel, existing_transverse = count_poles_in_hot_region(bs, along_v, hot_region)

        print(f"Existing poles in hot region: {existing_transverse} Transverse and {existing_parallel} Parallel")

        degree_parallel = bs.UDegree if not along_v else bs.VDegree
        degree_transverse = bs.VDegree if not along_v else bs.UDegree

        transverse_layers, parallel_layers = relax_layers_by_existing_dof(
            parallel_layers=parallel_layers,
            transverse_layers=transverse_layers,
            existing_parallel=existing_parallel,
            existing_transverse=existing_transverse,
            degree_parallel=degree_parallel,
            degree_transverse=degree_transverse
        )

        # ensure to add at least parallel + 1 transverse layers
        if parallel_layers > 0:
            transverse_layers = max(transverse_layers, max(1, parallel_layers + 1))

        print(f"Relaxed to {transverse_layers} transverse layers and {parallel_layers} parallel layers")


    (uDeg, vDeg) = needs_degree_elevation(bs, along_v, boundary_index, transverse_layers, parallel_layers)
    
    return {
        "transverse_layers": transverse_layers,
        "parallel_layers": parallel_layers,
        "hot_region": hot_region,
        "elevate_degree": (uDeg, vDeg)
    }

def enforce_G1_multiple(target_face, drivers, collision_mode='average', refinement=None):
    """
    Enforce G1 continuity on target_face along multiple driver surfaces.

    Args:
        target_face: Part.Face to modify
        drivers: list of dicts, each with keys:
            'driver_face': Part.Face,
            'samples': int,
            'beta': float,
            'spread_rows': int,
            'falloff': float
        collision_mode: 'average' or 'directional_clamp'
        refinement: dict with optional keys

    Returns:
        True if all edges processed, False if some failed.
    """
    # Get BSpline surface from target face   
    if not isinstance(target_face.Surface, Part.BSplineSurface):
        raise RuntimeError("Target face is not a BSpline surface")

    # Record the start time
    start_time = time()

    refine_surface = refinement is not None and refinement.get("method", None) is not None
    bs = target_face.Surface.copy()
    
    # track if surface modified during preprocessing - for debugging
    refined = False 

    # for uniform refinement we do it independent of drivers
    if refine_surface and refinement.get("method", None) == "uniform":
        min_u = refinement.get("min_poles_u", 8)
        min_v = refinement.get("min_poles_v", 8)
        tol = refinement.get("tol", 1e-9)
        degreeU, degreeV = bs.UDegree, bs.VDegree
        nu, nv = bs.NbUPoles, bs.NbVPoles
        # Ensure minimum degree 3 for better refinement
        if degreeU < 3:
            bs.increaseDegree(3, degreeV)  
            refined = True          
        if degreeV < 3:
            bs.increaseDegree(degreeU, 3)
            refined = True

        if nu < min_u:
            refined |= insert_uniform_knots(bs, along_v=False, target_poles=min_u, tol=tol)
        if nv < min_v:
            refined |= insert_uniform_knots(bs, along_v=True, target_poles=min_v, tol=tol)

    g1_poles = set()
    # Preprocess drivers
    for driver_idx, drv in enumerate(drivers):
        driver_face = drv['driver_face']
        samples = drv.get('samples', 30)

        # Find shared edge
        edge = find_shared_edge(target_face, driver_face)
        if edge is None:
            App.Console.PrintWarning("Failed to find shared edge with driver face\n")
            continue

        # Detect direction and boundary index 
        along_v, boundary_index = detect_edge_direction_and_boundary(bs, edge)
        side = 1 if boundary_index == 0 else -1

        if refine_surface:
            # check degree in refinement direction and if less than 3, elevate
            # along_v indicates boundary direction, so refinement is perpendicular
            degree = bs.VDegree if not along_v else bs.UDegree
            if degree < 3:
                if along_v:
                    bs.increaseDegree(3, bs.VDegree)
                else:
                    bs.increaseDegree(bs.UDegree, 3)
                refined = True

            refinement_method = refinement.get("method", "uniform")                     
            tol = refinement.get("tol", 1e-9)
            beta = drv.get('beta', 1.0)
           
            driver_dirs = get_surface_boundary_tangents(driver_face, edge, samples, None)
            target_dirs = get_surface_boundary_tangents(target_face, edge, samples, None)

            print(f"Estimating refinement for driver {driver_idx}: \n")

            while True:
                estimate = estimate_surface_refinement(
                    bs, driver_dirs, target_dirs, along_v, boundary_index,
                    beta=beta, energy_tol=refinement.get("tolerance", 0.15),
                    hot_angle_rad=refinement.get("hot_angle_rad", 0.2)
                )
                layers = estimate["transverse_layers"]
                parallel_layers = estimate["parallel_layers"]
                hot_region = estimate["hot_region"]
                (elevate_Udeg, elevate_Vdeg) = estimate["elevate_degree"]
                
                print(f"Current surface degree: UDegree={bs.UDegree}, VDegree={bs.VDegree}\n")

                if elevate_Udeg or elevate_Vdeg:
                    udeg = bs.UDegree
                    vdeg = bs.VDegree
                    if elevate_Udeg:
                        udeg += 1
                    if elevate_Vdeg:
                        vdeg += 1
                    bs.increaseDegree(udeg, vdeg)

                    boundary_index = 0 if side == 1 else (bs.NbUPoles - 1 if along_v else bs.NbVPoles - 1)
                    refined = True
                else:
                    break

            if parallel_layers > 0 and hot_region:
                refined |= insert_localized_parallel_knots(bs, along_v, hot_region, parallel_layers, tol)
            
            if refinement_method == "weighted":
                # compute weights based on driver tangents
                refined |= insert_mismatch_weighted_knots(bs, not along_v, boundary_index, edge,
                            driver_dirs, target_dirs, layers, tol)
            else:
                boundary_bias = refinement.get("boundary_bias", 1.2)                    
                refined |= insert_boundary_biased_knots(bs, not along_v, boundary_index, 
                    layers, boundary_bias, tol)

        # After refinement (if any) boundary index may have changed
        boundary_index = 0 if side == 1 else (bs.NbUPoles - 1 if along_v else bs.NbVPoles - 1)
        
        g1_index = boundary_index + side

        if along_v:
            for j in range(1, bs.NbVPoles-1):  # skip boundaries
                g1_poles.add((g1_index, j))
        else:
            for i in range(1, bs.NbUPoles-1):  # skip boundaries
                g1_poles.add((i, g1_index))

        fairing_distance = drv.get('fairing_distance', 0.6)
        spread = max(0, int(fairing_distance * (bs.NbUPoles if along_v else bs.NbVPoles)) - 2) # exclude boundaries and g1 row        
        spread_rows = get_spread_rows(boundary_index+side, max(0, spread), side)
        
        decay_type = drv.get('decay_type', 'smooth')
        if decay_type == 'linear':
            decay = lambda d, md=spread: linear_decay(d, md)
        elif decay_type == 'gaussian':
            sigma = drv.get('gaussian_sigma', spread / 3.0)
            decay = lambda d, s=sigma: gaussian_decay(d, s)
        else:
            decay = lambda d, md=spread: smooth_decay(d, md)

        drv['edge'] = edge                              # shared edge
        drv['along_v'] = along_v                        # boundary direction
        drv['side'] = side                              # +1 if boundary_index == 0 else -1 
        drv['spread_rows'] = spread_rows                # list of rows to spread G1 influence
        drv['boundary_index'] = boundary_index          # updated boundary index after refinement
        drv['decay'] = decay                            # decay function
    
    # Debug: show modified surface after refinement
    # if refined:
    # show_surface(bs, "Refined_Surface")

    is_rational = isRational(bs)

    # Save original poles and weights
    orig_poles = bs.getPoles()    
    orig_weights = None
    if is_rational:
        orig_weights = bs.getWeights()

    # Accumulators for pole deltas
    # (i, j): {
    #     "parallel": None,          # strongest projected delta
    #     "orthogonal": [],          # list of residual deltas
    #     "fairing": [],             # fairing-only deltas
    #     "raw": [],                 # raw deltas for debugging or averaging
    #     "weights": []              # rational only
    # }
    delta_accum = {}
    # Desired transverse directions per boundary pole (before spreading)
    # key = (u,v) - pole indices
    desired_dirs = defaultdict(list)
    
    # Process each driver
    for driver_idx, drv in enumerate(drivers):
        driver_face = drv['driver_face']
        samples = drv.get('samples', 30)
        beta = drv.get('beta', 1.0)
        edge = drv['edge']
                
        if edge is None:
            # Find shared edge
            edge = find_shared_edge(target_face, driver_face)

        if edge is None:
            App.Console.PrintWarning(f"Failed to find shared edge with driver face. Driver {driver_idx} will be ignored.\n")
            continue

        boundary_index = drv['boundary_index']
        along_v = drv['along_v']
        spread_rows = drv.get('spread_rows', [])
        side = drv['side']
        decay = drv['decay']
        
        # Get exact boundary pole parameters
        pole_params = boundary_pole_params(bs, along_v, boundary_index)

        # Sample driver tangents aligned to those poles
        driver_tangents_dir = get_surface_boundary_tangents(driver_face, edge, samples, pole_params)

        row_len = bs.NbVPoles if along_v else bs.NbUPoles

        if row_len != len(driver_tangents_dir):
            App.Console.PrintWarning(f"Mismatch row length. tangents length: {len(driver_tangents_dir)}, row length: {row_len}\n")
            continue

        # Compute new G1 row
        g1_row = compute_g1_row(orig_poles, driver_tangents_dir, along_v, boundary_index, beta)

        # Next to boundary row        
        adj_row = boundary_index + side

        # Capture desired transverse directions and intents to move poles
        for j in range(1, row_len-1):  # skip boundaries
            # Get adjacent row pole
            if along_v:
                i, jj = adj_row, j
            else:
                i, jj = j, adj_row

            P_adj = orig_poles[i][jj]

            # G1 target point
            P_g1 = App.Vector(g1_row[j])

            raw_delta = P_g1.sub(P_adj)
            mag = raw_delta.Length

            if mag < 1e-12:
                continue

            delta = App.Vector(raw_delta)
            dir = raw_delta.normalize()

            key = (i, jj)
            g1_poles.add(key)

            desired_dirs[key].append({
                "dir": dir,
                "mag": mag,
                "along_v": along_v,
                "boundary": boundary_index,
                "driver_idx": driver_idx,
                "priority": 0           # high priority for adjacent row
            })

            # Spread to interior rows
            for r in spread_rows:
                if decay is None:
                    continue
                
                d = abs(r - adj_row)
                d_mag = mag*decay(d)

                if d_mag < 1e-6:
                    continue

                if along_v:
                    key = (r, j)
                else:
                    key = (j, r)
                
                if key in g1_poles:                    
                    continue  # do not override g1 poles

                desired_dirs[key].append({
                    "dir": dir,
                    "mag": d_mag,
                    "along_v": along_v,
                    "boundary": boundary_index,
                    "driver_idx": driver_idx,
                    "priority": d           # further from g1 rows get lower priority based on distance
                })


    pole_constraints = resolve_constraints(desired_dirs)
    pole_basis = build_pole_basis(bs)
    
    def add_delta(i, j, delta, weight=None):
        if delta.Length < 1e-9:
            return
        
        # apply pole constraint
        constraint = pole_constraints.get((i, j))
        if constraint:
            delta = apply_constraint(delta, constraint)
            if delta.Length < 1e-9:
                return

        key = (i, j)
        if key not in delta_accum:
            delta_accum[key] = {
                "parallel": None,
                "orthogonal": [],
                "fairing": [],
                "raw": [],
                "weights": []
            }

        acc = delta_accum[key]
        acc["raw"].append(App.Vector(delta))

        if is_rational and weight is not None:
            acc["weights"].append(weight)

        basis = pole_basis[i][j]

        # no basis: pure accumulation
        if not basis:
            acc["orthogonal"].append(delta)
            return
        
        # Delta decomposition

        # project into basis
        d_proj = App.Vector(0, 0, 0)
        for b in basis:
            d_proj = d_proj.add(App.Vector(b).multiply(delta.dot(b)))

        d_orth = delta.sub(d_proj)

        # parallel: keep strongest projection
        cur = acc["parallel"]
        if cur is None or d_proj.Length > cur.Length:
            acc["parallel"] = d_proj

        # orthogonal: store residuals
        if d_orth.Length > 1e-9:
            acc["orthogonal"].append(d_orth)

    # Accumulate deltas from all intents and store boundary poles
    for key, items in desired_dirs.items():
        i, j = key        
        for item in items:
            dir = item["dir"]
            mag = item["mag"]
            delta = App.Vector(dir).multiply(mag)
            if is_rational:
                add_delta(i, j, delta, orig_weights[i][j])
            else:
                add_delta(i, j, delta)

    # Apply accumulated deltas to row poles solving conflicts

    # Compute final poles by averaging deltas
    if collision_mode == 'average':
        for (i, j), acc in delta_accum.items():
            new_pole = orig_poles[i][j].add(v_avg(acc["raw"]))
            bs.setPole(i + 1, j + 1, new_pole)

            if is_rational:
                avg_weight = orig_weights[i][j]
                if acc["weights"] and len(acc["weights"]) > 0:
                    avg_weight = w_avg(acc["weights"])
                bs.setWeight(i + 1, j + 1, avg_weight)
    # Compute final poles by derectionaly averaging deltas
    elif collision_mode == 'directional_clamp':
        for (i, j), acc in delta_accum.items():
            d = App.Vector(0, 0, 0)

            if acc["parallel"]:
                d = d.add(acc["parallel"])

            if acc["orthogonal"] and len(acc["orthogonal"]) > 0:
                d = d.add(v_avg(acc["orthogonal"]))

            if d.Length < 1e-9:
                continue

            bs.setPole(i + 1, j + 1, orig_poles[i][j].add(d))

            if is_rational:
                avg_weight = orig_weights[i][j]
                if acc["weights"] and len(acc["weights"]) > 0:
                    avg_weight = w_avg(acc["weights"])
                bs.setWeight(i + 1, j + 1, avg_weight)
    
    # draw_points(deltas, (1.0, 0.0, 1.0), 0.2, "Faired_Poles")

    # Record the end time
    end_time = time()

    # Calculate and print the elapsed time
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")

    return bs

def getSelectedFaces():
    sel = App.Gui.Selection.getSelectionEx()
    faces = []
    for s in sel:
        for fi in s.SubElementNames:
            if fi.startswith("Face"):
                face_index = int(fi[4:])  # "Face1" -> 1
                face = s.Object.Shape.Faces[face_index - 1]
                faces.append(face)
    return faces

# --------------------------------
# USER INPUT
# --------------------------------
selection_faces = getSelectedFaces()

if len(selection_faces) < 2:
    raise RuntimeError("Select at least two faces: target face first, driver face(s) second")
target_face = selection_faces[0]
driver_faces = selection_faces[1:]

drivers = []
for f in driver_faces:
    drivers.append({
        'driver_face': f,       # tool face
        'samples': 50,          # (optional, default 30) number of samples along edge
        'beta': 1.0,            # g1 factor (0.0 - 1.0). Higher = stronger driver influence.
        'decay_type': 'smooth', # 'linear', 'smooth', 'gaussian' (optional, default 'smooth')
        'fairing_distance': 0.6 # (optional, default 0.8) fairing distance as fraction of pole count along fairing direction
    })

# Optional: refinement before G1 enforcement
g1_refinement = {
    "method": "weighted",   # "uniform", "biased", "weighted" (None to disable)
    "min_poles_u": 8,       # (optional, default 8) for "uniform" method
    "min_poles_v": 8,       # (optional, default 8) for "uniform" method
    "boundary_bias": 1.2,   # for "biased" method
    "layers": 3,            # (optional, will be estimated) number of knot insertion layers for "biased" and "weighted" methods (overrides automatic estimation)
    "tol": 1e-9             # (optional, default 1e-9) knot insertion tolerance
}

#bs = enforce_G1_multiple(target_face, drivers, collision_mode='average') # more gentle corners blending
bs = enforce_G1_multiple(target_face, drivers, collision_mode='directional_clamp', refinement=g1_refinement) # better G1 in tight corners

new_face = Part.Face(bs)

doc = App.ActiveDocument
obj = doc.addObject("Part::Feature", "G1_ConstrainedSurface")
obj.Shape = new_face
doc.recompute()