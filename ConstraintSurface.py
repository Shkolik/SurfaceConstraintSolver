# SPDX-License-Identifier: LGPL-2.1-or-later

__title__ = "Surface tangency solver"
__author__ = "Andrew Shkolik"
__license__ = "LGPL 2.1+"

from collections import defaultdict
import math
import FreeCAD as App
import Part

### Debug visualization utilities ###

def draw_vector(point, vector, color=(1.0,0,0), scale=10.0, name_prefix="Vec", show=True):
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

    if show:
        obj = App.ActiveDocument.addObject("Part::Feature", name_prefix)
        obj.Shape = edge
        obj.ViewObject.ShapeColor = color
        obj.ViewObject.LineColor = color
        obj.ViewObject.PointColor = color
    return edge

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
        Part.Compound or None: The Compound of spheres or None if not shown
    """    
    compound = Part.makeCompound([Part.makeSphere(size, point) for point in points])

    if show:
        obj = App.ActiveDocument.addObject("Part::Feature", name)
        obj.Shape = compound
        obj.ViewObject.ShapeColor = color
        obj.ViewObject.LineColor = color
        obj.ViewObject.PointColor = color
    return compound

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

def transverse_curvature_proxy(surface, point, transverse_dir, tol=1e-9):
    """
    Approximate geometric curvature magnitude in transverse direction.
    Args:
        surface: (Part.Surface) driver surface
        point: (App.Vector) origin point
        transverse_dir: (App.Vector) direction (unit)
        tol: (float) tolerance for first derivative
    """    
    u, v = surface.parameter(point)

    # First derivatives (world space)
    Du = surface.getDN(u, v, 1, 0)
    Dv = surface.getDN(u, v, 0, 1)

    # Second derivatives
    Duu = surface.getDN(u, v, 2, 0)
    Dvv = surface.getDN(u, v, 0, 2)
    Duv = surface.getDN(u, v, 1, 1)

    # Solve Td ≈ a*Du + b*Dv (least squares)

    # Normal equation coefficients
    suu = Du.dot(Du)
    svv = Dv.dot(Dv)
    suv = Du.dot(Dv)

    rhs_u = transverse_dir.dot(Du)
    rhs_v = transverse_dir.dot(Dv)

    det = suu * svv - suv * suv
    if abs(det) < 1e-12:
        # Fallback: Project onto dominant direction
        if Du.Length > Dv.Length:
            d1 = Du
            d2 = Duu
        else:
            d1 = Dv
            d2 = Dvv

        if d1.Length < tol:
            return 0.0

        return d2.Length / (d1.Length * d1.Length)

    inv = 1.0 / det
    a = inv * ( rhs_u * svv - rhs_v * suv )
    b = inv * ( rhs_v * suu - rhs_u * suv )

    # First derivative along transverse_dir
    d1 = Du.multiply(a).add(Dv.multiply(b))
    d1_len = d1.Length
    if d1_len < tol:
        return 0.0

    # Second derivative along transverse_dir
    d2 = (
        Duu.multiply(a * a)
        .add(Dvv.multiply(b * b))
        .add(Duv.multiply(2.0 * a * b))
    )

    # Geometric curvature magnitude
    return d2.Length / (d1_len * d1_len)

def sample_boundary_curvature(surface, edge, transverse_dirs):
    """
    Returns curvature values per boundary sample.
    Args:
        surface: (Part.Surface) driver surface
        edge: (Part.Edge) shared edge
        transverse_dirs: list of (App.Vector) (unit) or None
    Returns:    
        curvatures: list of float curvature values
    """
    u0, u1 = edge.ParameterRange
    samples = len(transverse_dirs)
    curvatures = []

    for i in range(samples):
        t = u0 + (u1 - u0) * i / (samples - 1)
        p = edge.Curve.value(t)

        Td = transverse_dirs[i]
        if Td is None:
            curvatures.append(0.0)
            continue

        k = transverse_curvature_proxy(surface, p, Td)
        curvatures.append(k)

    return curvatures

def curvature_weighted_params(params, curvatures, layers=3, min_weight=0.1):
    """
    Generate parameter locations biased by curvature.
    Args:
        params: list of (float) parameter locations
        curvatures: list of (float) curvature values
        layers: (int) number of layers to insert
        min_weight: (float) minimum weight to avoid zero insertions
    Returns:
        list of (float) new parameter locations
    """
    if not params or not curvatures:
        return []

    # Normalize curvature
    max_k = max(curvatures)
    if max_k < 1e-9:
        return []

    weights = [max(k / max_k, min_weight) for k in curvatures]

    selected = []
    used = set()

    for _ in range(layers):
        # Pick strongest unused location
        best_i = None
        best_k = -1

        for i, k in enumerate(weights):
            if i in used:
                continue
            if k > best_k:
                best_k = k
                best_i = i

        if best_i is None:
            break

        selected.append(params[best_i])
        used.add(best_i)

    return selected

def insert_curvature_weighted_knots(bs, along_v, boundary_index, 
                                    edge, curvatures, layers=3, tol=1e-9):
    """
    Insert knots near boundary weighted by transverse curvature.
    Args:
        bs: (Part.BSplineSurface) surface to refine (MODIFIED IN PLACE)
        along_v: (bool) True → boundary is U=const (refine V)
        boundary_index: (int) 0 or max index of boundary row/column
        edge: (Part.Edge) shared edge
        curvatures: list of (float) geometric curvature values per sample
        layers: (int) number of layers to insert
        tol: (float) knot insertion tolerance
    """
    if curvatures is None or len(curvatures) == 0:
        return False
    
    samples = len(curvatures)

    # Direction setup
    if along_v:
        insert = lambda t: bs.insertVKnot(t, 1, tol)
        param_min, param_max = bs.bounds()[2:4]
    else:
        insert = lambda t: bs.insertUKnot(t, 1, tol)
        param_min, param_max = bs.bounds()[0:2]

    u0, u1 = edge.ParameterRange
    params = [u0 + (u1 - u0) * i / (samples - 1) for i in range(samples)]

    # Compute insert locations
    refine_params = curvature_weighted_params(
        params,
        curvatures,
        layers=layers
    )

    # Map edge parameter → surface parameter
    eps = (param_max - param_min) * 1e-6
    at_start = (boundary_index == 0)
    boundary_param = param_min if at_start else param_max

    modified = False
    for t in refine_params:
        s = (t - u0) / (u1 - u0)
        d = s * (param_max - param_min) * 0.5

        if at_start:
            u = boundary_param + d
        else:
            u = boundary_param - d

        u = max(param_min + eps, min(param_max - eps, u))
        insert(u)
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

def sample_driver_tangents(driver_face, edge, samples, row_len):
    """
    Sample driver surface tangents and interpolate to match target row poles.
    
    Args:
        driver_face: Driver surface (Part.Face)
        edge: Shared edge (Part.Edge)
        samples: Number of samples along the edge
        row_len: Number of poles in the target row. If None, no interpolation will be done.
    Returns:
        list of App.Vector (unit) or None
    """
    u0, u1 = edge.ParameterRange

    tangents = [] # tangent directions at samples

    for i in range(samples):
        if i == samples - 1:
            t = u1 - 1e-9  # small offset before the end
        elif i == 0:
            t = u0 + 1e-9  # small offset after the start
        else:
            t = u0 + (u1 - u0) * i / (samples - 1)
        
        p = edge.Curve.value(t)

        Te = edge.Curve.tangent(t)[0]
        if Te.Length < 1e-9:
            tangents.append(None)
            continue
        Te.normalize()
        
        d_dir = find_transverse_tangent(driver_face, p, Te)

        if d_dir is None:
            tangents.append(None)
            continue

        tangents.append(d_dir)

    if not row_len:
        return tangents
    
    # Interpolate to match row_len
    tangent_at_pole_dir = []

    for j in range(row_len):
        s = j * (samples - 1) / (row_len - 1)
        i0 = int(s)
        i1 = min(i0 + 1, samples - 1)
        w = s - i0

        T0 = tangents[i0]
        T1 = tangents[i1]

        if T0 is None and T1 is None:
            raise RuntimeError("No valid driver tangents")

        if T0 is None:
            T = App.Vector(T1)
        elif T1 is None:
            T = App.Vector(T0)
        else:
            if T0.dot(T1) < 0:
                T1 = T1.negative()
            T = lerp(T0, T1, w)

        if T.Length < 1e-9:
            raise RuntimeError("Failed to compute valid tangent")

        T.normalize()
        tangent_at_pole_dir.append(T)
        
    return tangent_at_pole_dir

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

def get_spread_rows(boundary_index, spread_rows, nu, nv, along_v):
    """
    Find interior rows to apply G1 spreading.

    Args:
        boundary_index: index of boundary row/column (0-based)
        spread_rows: requested number of interior rows
        nu, nv: number of poles in U and V
        along_v: True if boundary is U=const, False if V=const

    Returns:
        list of row indices (0-based), ordered from edge inward
    """

    num_rows = nu if along_v else nv

    max_safe = max(1, num_rows // 2)
    spread = max(1, min(spread_rows, max_safe))

    # Determine direction
    if boundary_index == 0:
        step = 1
        r = 1
    else:
        step = -1
        r = boundary_index - 1

    # Collect rows
    rows = []

    for _ in range(spread):
        if not (0 <= r < num_rows):
            break
        rows.append(r)
        r += step

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

def estimate_refinement_layers(h, kappa, tol=0.25):
    """
    Estimate the number of refinement layers needed based on pole spacing and curvature.
    
    Args:
        h: (float) Current pole spacing near boundary
        kappa: (float) Curvature magnitude
        tol: (float) Tolerance for acceptable deviation
    Returns: number of layers to insert
    """
    if kappa < 1e-12:  # avoid division by zero, or very small curvature
        return 0
    h_target = (tol / kappa) ** 0.5
    layers = max(0, math.ceil(math.log2(h / h_target)))
    return layers

def can_apply_weighted_refinement(surface, alongV):
    """
    Check if surface can be refined using curvature-weighted method.
    
    Args:
        surface: target BSpline surface
        alongV: True if refining along V, False if U
    Returns:
        True if refinement can be applied, False otherwise
    """
    degree = surface.VDegree if alongV else surface.UDegree
    poles_count = surface.NbVPoles if alongV else surface.NbUPoles
    # Minimum degree and poles count for meaningful weighted refinement
    if degree < 3 or poles_count < 4:
        return False
    return True

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
                "driver_idx": int
            }

    Returns:
        pole_constraints: dict[(i,j)] -> constraint dict
    """

    pole_constraints = {}

    for key, entries in desired_dirs.items():

        # --- sanitize ---
        dirs = []
        mags = []

        for e in entries:
            d = App.Vector(e["dir"])
            if d.Length < 1e-9:
                continue
            d.normalize()
            dirs.append(d)
            mags.append(abs(e["mag"]))

        if not dirs:
            continue

        # --- single driver: trivial ---
        if len(dirs) == 1:
            pole_constraints[key] = {
                "mode": "single",
                "dir": dirs[0],
                "max_mag": mags[0]
            }
            continue

        # --- cluster directions by alignment ---
        clusters = []  # list of {"dir": vec, "mags": [], "dirs": []}

        for d, m in zip(dirs, mags):
            placed = False
            for c in clusters:
                if abs(c["dir"].dot(d)) > 1.0 - 1e-6:
                    c["mags"].append(m)
                    c["dirs"].append(d)
                    placed = True
                    break
            if not placed:
                clusters.append({
                    "dir": App.Vector(d),
                    "dirs": [d],
                    "mags": [m]
                })

        # --- one dominant cluster ---
        if len(clusters) == 1:
            c = clusters[0]
            pole_constraints[key] = {
                "mode": "single",
                "dir": c["dir"],
                "max_mag": max(c["mags"])
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
            v = v.sub(u.multiply(v.dot(u)))
            if v.Length > 1e-9:
                v.normalize()

                pole_constraints[key] = {
                    "mode": "plane",
                    "basis": [u, v],
                    "max_mag": min(max(c1["mags"]), max(c2["mags"]))
                }
                continue

        # --- fallback: clamp ---
        # Multiple conflicting directions, low DOF
        # Allow movement but cap magnitude
        avg_dir = App.Vector(0, 0, 0)
        for d in dirs:
            avg_dir = avg_dir.add(d)

        if avg_dir.Length > 1e-9:
            avg_dir.normalize()

        pole_constraints[key] = {
            "mode": "clamp",
            "dir": avg_dir,
            "max_mag": min(mags)
        }

    return pole_constraints

def apply_constraint(delta, c):
    if delta.Length < 1e-9:
        return delta

    if c["mode"] == "single":
        d = App.Vector(c["dir"])
        mag = delta.dot(d)
        mag = max(-c["max_mag"], min(c["max_mag"], mag))
        return d.multiply(mag)

    if c["mode"] == "plane":
        proj = App.Vector(0, 0, 0)
        for b in c["basis"]:
            proj = proj.add(b.multiply(delta.dot(b)))

        if proj.Length > c["max_mag"]:
            proj.normalize()
            proj = proj.multiply(c["max_mag"])
        return proj

    if c["mode"] == "clamp":
        if delta.Length > c["max_mag"]:
            delta = App.Vector(delta)
            delta.normalize()
            delta = delta.multiply(c["max_mag"])
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

def compute_distance_to_fixed(nu, nv, fixed):
    dist = [[1e9]*nv for _ in range(nu)]

    from collections import deque
    q = deque()

    for (i,j) in fixed:
        dist[i][j] = 0
        q.append((i,j))

    while q:
        i,j = q.popleft()
        for di,dj in ((1,0),(-1,0),(0,1),(0,-1)):
            ni,nj = i+di, j+dj
            if 0 <= ni < nu and 0 <= nj < nv:
                if dist[ni][nj] > dist[i][j] + 1:
                    dist[ni][nj] = dist[i][j] + 1
                    q.append((ni,nj))
    return dist

def extract_g1_indices(locked_poles, axis):
    idx = set()
    for (i, j) in locked_poles.keys():
        idx.add(i if axis == 'U' else j)
    return idx

def fair_surface(poles, drivers, fixed_poles, g1_poles, dist, add_delta):
    nu = len(poles)
    nv = len(poles[0])

    for driver in drivers:
        axis = driver.axis        # 'U' or 'V'
        side = driver.side        # +1 or -1
        max_dist = driver.max_dist
        strength = driver.strength
        decay = driver.decay

        for i in range(nu):
            for j in range(nv):

                idx = (i, j)

                # do not modify fixed or G1 poles
                if idx in fixed_poles:
                    continue
                if idx in g1_poles:
                    continue
                
                # do not modify poles beyond max distance
                d = dist[i][j]
                if d <= 0 or d > max_dist:
                    continue

                # directional neighbors
                if axis == 'U':
                    # driver along U → fair in V
                    jp = j + side
                    jm = j - side
                    if jp < 0 or jp >= nv or jm < 0 or jm >= nv:
                        continue

                    P  = poles[i][j]
                    Pp = poles[i][jp]
                    Pm = poles[i][jm]
                else:
                    # driver along V → fair in U
                    ip = i + side
                    im = i - side
                    if ip < 0 or ip >= nu or im < 0 or im >= nu:
                        continue

                    P  = poles[i][j]
                    Pp = poles[ip][j]
                    Pm = poles[im][j]

                # 1D Laplacian
                lap = (Pp + Pm) * 0.5 - P

                # distance-weighted decay
                w = decay(d)
                if w <= 0.0:
                    continue

                delta = lap * (strength * w)

                add_delta(i, j, delta)

def propagate_g1_fairing(g1_intents, drivers, nu, nv, add_delta):
    for drv in drivers:
        axis = drv.axis
        max_dist = drv.max_dist
        decay = drv.decay
        side = drv.side

        for intent in g1_intents:
            i0, j0 = intent['pole']
            delta0 = intent['delta']

            for d in range(1, max_dist + 1):
                w = decay(d)
                if w < 1e-6:
                    continue

                if axis == 'U':
                    i = i0 + d * side
                    j = j0
                else:
                    i = i0
                    j = j0 + d * side

                # stay within interior poles only
                if not (1 <= i < nu-1 and 1 <= j < nv-1):
                    break

                add_delta(i, j, delta0, w, d)

def linear_decay(d, max_d):
    return max(0.0, 1.0 - d / max_d)

def smooth_decay(d, max_d):
    t = min(1.0, d / max_d)
    return 1.0 - (3*t*t - 2*t*t*t)

def gaussian_decay(d, sigma):
    return math.exp(-(d*d) / (2*sigma*sigma))
              
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

    # Preprocess drivers
    for drv in drivers:
        driver_face = drv['driver_face']
        samples = drv.get('samples', 30)

        # Find shared edge
        edge = find_shared_edge(target_face, driver_face)
        if edge is None:
            App.Console.PrintWarning("Failed to find shared edge with driver face\n")
            continue

        # Detect direction and boundary index 
        along_v, boundary_index = detect_edge_direction_and_boundary(bs, edge)

        drv['edge'] = edge
        drv['along_v'] = along_v
        drv['side'] = 1 if boundary_index == 0 else -1
        drv['axis'] = 'U' if along_v else 'V'

        if refine_surface:
            needs_refinement = False # even if user requested, doesnt mean we need it

            # check if weighted refinement is possible before degree elevation
            # weighted_refinement inserts knots based on curvature, and if initial degree is too low, 
            # it may not work well even after elevation -> curvature will stay the same
            weighted_possible = can_apply_weighted_refinement(bs, along_v)

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

            # Sample tangents from driver face along edge, do not interpolate yet
            driver_tangents_dir = sample_driver_tangents(driver_face, edge, samples, None)

            # check curvature vs pole spacing
            h = pole_spacing(bs, along_v, boundary_index)
            curvatures = sample_boundary_curvature(driver_face.Surface, edge, driver_tangents_dir)
            kappa = max(curvatures)

            if kappa * h * h > 0.25:
                needs_refinement = True

            if needs_refinement:
                layers = max(refinement.get("layers", 0), estimate_refinement_layers(h, kappa, tol=0.25))

                if layers == 0:
                    App.Console.PrintWarning("Refinement needed but no layers estimated or provided\n")
                    continue  # no refinement needed or possible

                if refinement_method == "weighted" and weighted_possible:
                    # compute weights based on driver tangents
                    refined |= insert_curvature_weighted_knots(bs, not along_v, boundary_index, 
                        edge, curvatures, layers, tol)
                else:
                    boundary_bias = refinement.get("boundary_bias", 1.2)                    
                    refined |= insert_boundary_biased_knots(bs, not along_v, boundary_index, 
                        layers, boundary_bias, tol)
                              
    # Debug: show modified surface after refinement
    # if refined:
    # show_surface(bs, "Refined_Surface")

    nu, nv = bs.NbUPoles, bs.NbVPoles
    is_rational = isRational(bs)

    # Save original poles and weights
    orig_poles = [[bs.getPole(i+1, j+1) for j in range(nv)] for i in range(nu)]
    orig_weights = None
    if is_rational:
        orig_weights = [[bs.getWeight(i+1, j+1) for j in range(nv)] for i in range(nu)]

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
    g1_intents = []   # list of dicts

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
            App.Console.PrintWarning("Failed to find shared edge with driver face\n")
            continue

        # Detect direction and boundary index 
        along_v, boundary_index = detect_edge_direction_and_boundary(bs, edge)

        row_len = nv if along_v else nu

        # Sample tangents from driver face along edge (reuse your sampling code)
        driver_tangents_dir = sample_driver_tangents(driver_face, edge, samples, row_len)

        # Compute new G1 row
        g1_row = compute_g1_row(orig_poles, driver_tangents_dir, along_v, boundary_index, beta)

        # Next to boundary row        
        adj_row = boundary_index + (1 if boundary_index == 0 else -1)

        # Capture desired transverse directions and intents to move poles
        for j in range(row_len):
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

            g1_intents.append({
                "pole": (i, jj),            # adjacent row pole
                "delta": delta,
                "boundary": boundary_index, # boundary row index (0-based)
                "along_v": along_v          # boundary direction
            })
            
            desired_dirs[key].append({
                "dir": dir,
                "mag": mag,
                "along_v": along_v,
                "boundary": boundary_index,
                "driver_idx": driver_idx
            })

    pole_constraints = resolve_constraints(desired_dirs)
    pole_basis = build_pole_basis(bs)
    
    def add_delta(i, j, delta,weight=None):
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

    # fixed poles set (boundary poles)
    fixed = set()
    for v in [0, nv-1]:
        for i in range(nu):
            fixed.add((i, v))
    for u in [0, nu-1]:
        for j in range(nv):
            fixed.add((u, j))

    # Accumulate deltas from all intents and store boundary poles
    for intent in g1_intents:
        i, j = intent["pole"]

        delta = App.Vector(intent["delta"])
        if is_rational:
            add_delta(i, j, delta, orig_weights[i][j])
        else:
            add_delta(i, j, delta)

    # Track moved (G1) poles for fairing, exclude boundaries
    # set of (i,j)
    moved = set()    
    for (i, j), acc in delta_accum.items():
        if acc["parallel"] or (acc["orthogonal"] and len(acc["orthogonal"]) > 0):
            moved.add((i, j))

    def add_fairing_delta(i, j, delta, weight=1.0, d = 0):
        if delta.Length < 1e-9:
            return

        # never move fixed / constrained poles
        if (i, j) in fixed or (i, j) in moved:
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

        # print(f"Applying propagated fairing delta to pole({i}, {j}) with magnitude {delta.Length}; decay weight {weight}; distance {d}")
               

        # fairing is always orthogonal to constraints
        delta_accum[key]["fairing"].append(App.Vector(delta).multiply(weight))

    # Apply accumulated deltas to adjacent row poles solving conflicts

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
    
    # show_surface(bs, "G1_Enforced_Surface_Raw")

    # Compute distance field from fixed poles (boundary and constrained)
    dist = compute_distance_to_fixed(nu, nv, moved)
    # print("Distance field:")
    # for j in range(nv):
    #     print([dist[i][j] for i in range(nu)])

    # gather fairing drivers
    drivers_fairing = []
    for drv in drivers:
        axis = drv.get('axis', 'U')
        side = drv.get('side', 1)
        spread = int(1.0 * max([d for rows in dist for d in rows]))
        strength = drv.get('fairing_strength', 0.8)

        driver = FairingDriver(
            axis=axis,
            side=side,
            max_dist=spread,
            strength=strength,
            decay=lambda d, md=spread: smooth_decay(d, md)
        )

        drivers_fairing.append(driver)

    propagate_g1_fairing(g1_intents, drivers_fairing, nu, nv, add_fairing_delta)
    # deltas = []
    for (i, j), acc in delta_accum.items():
        d = App.Vector(0, 0, 0)
        if acc["fairing"] and len(acc["fairing"]) > 0:
            d = d.add(v_avg(acc["fairing"]))

        if d.Length < 1e-9:
            continue
        # deltas.append(orig_poles[i][j].add(d))
        bs.setPole(i + 1, j + 1, orig_poles[i][j].add(d))
    
    # draw_points(deltas, (1.0, 0.0, 1.0), 0.2, "Faired_Poles")

    # Laplacian based fairing - save for now just in case

    # Run fairing multiple times for better effect
    # for _ in range(3):    
    #     poles = [[bs.getPole(i+1, j+1) for j in range(nv)] for i in range(nu)]

    #     fair_surface(poles, drivers_fairing, fixed, moved, dist, add_delta=add_fairing_delta)

    #     # apply fairing deltas, just an average
    #     for (i, j), acc in delta_accum.items():
    #         d = App.Vector(0, 0, 0)
    #         # --- fairing (soft, damped) ---
    #         if acc["fairing"] and len(acc["fairing"]) > 0:
    #             d = d.add(v_avg(acc["fairing"]))

    #         if d.Length < 1e-9:
    #             continue

    #         bs.setPole(i + 1, j + 1, orig_poles[i][j].add(d))
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
        'driver_face': f,   # tool face
        'samples': 30,      # (optional, default 30) number of samples along edge
        'beta': 1.0,        # blending factor (0.0 - 1.0). Higher = stronger driver influence.        
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