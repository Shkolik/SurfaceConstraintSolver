# SPDX-License-Identifier: LGPL-2.1-or-later

__title__ = "Surface tangency solver"
__author__ = "Andrew Shkolik"
__license__ = "LGPL 2.1+"

import FreeCAD as App
import Part

def draw_vector(point, vector, color=(1.0,0,0), scale=10.0, name_prefix="Vec"):
    """
    Method to draw a vector as a line in FreeCAD for visualization.
    
    Args:
    point: Starting point (App.Vector)
    vector: Vector to draw (App.Vector)
    color: Line color (r,g,b)
    scale: Scale factor for vector length if vector is unit length
    name_prefix: Name prefix for the object
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

    obj = App.ActiveDocument.addObject("Part::Feature", name_prefix)
    obj.Shape = edge
    obj.ViewObject.ShapeColor = color
    obj.ViewObject.LineColor = color
    obj.ViewObject.PointColor = color
    return obj

def draw_points(points, color=(1.0,0,0), size=2.0, name="Pt"):    
    """
    Method to draw points as spheres in FreeCAD for visualization.

    Args:
    points: List of points (App.Vector)
    color: Sphere color (r,g,b)
    size: Sphere radius
    name: Name for the object
    """
    obj = App.ActiveDocument.addObject("Part::Feature", name)
    obj.Shape = Part.makeCompound([Part.makeSphere(size, point) for point in points])
    obj.ViewObject.ShapeColor = color
    obj.ViewObject.LineColor = color
    obj.ViewObject.PointColor = color
    return obj


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

def insert_boundary_biased_knots(bs, along_v, boundary_index, target_poles=8, bias=3.0, tol=1e-9):
    """
    Insert knots clustered near a boundary to increase DOF
    without deforming the surface.

    Args:
        bs (Part.BSplineSurface)
        along_v (bool): True → boundary is U=const (refine V)
        boundary_index (int): 0 or max index of boundary row/column
        target_poles (int): minimum number of poles in that direction        
        bias (float): >1 clusters knots near boundary
        tol (float) : knot insertion tolerance
    Returns:
        bool: True if the surface was modified, False otherwise
    """
    modified = False

    # --- Direction setup ---
    if along_v:
        # Boundary is U=const → refine V
        get_nb_poles = lambda: bs.NbVPoles
        get_knots   = lambda:  bs.getVKnots()
        get_mults   = lambda: bs.getVMultiplicities()
        insert_knot = lambda u: bs.insertVKnot(u, 1, tol)
        param_min, param_max = bs.bounds()[2:4]
    else:
        # Boundary is V=const → refine U
        get_nb_poles = lambda: bs.NbUPoles
        get_knots   = lambda: bs.getUKnots()
        get_mults   = lambda: bs.getUMultiplicities()
        insert_knot = lambda u: bs.insertUKnot(u, 1, tol)
        param_min, param_max = bs.bounds()[0:2]

    cur_poles = get_nb_poles()
    if cur_poles >= target_poles:
        return  modified # nothing to do
    
    layers = target_poles - cur_poles

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

def ensure_min_dof(target_face, drivers, min_u=8, min_v=8, boundary_bias=1.0, tol=1e-9):
    """
    Ensure the BSpline surface has at least min_u and min_v poles.
    Args:
        target_face: Part.Face
        drivers: List of driver dictionaries
        min_u: Minimum number of U poles
        min_v: Minimum number of V poles
        boundary_bias: Bias toward boundary
        tol: Tolerance for knot insertion
    Returns: Part.BSplineSurface with ensured minimum degrees of freedom
    """
    bs = target_face.Surface.copy()
    modified = False

    if boundary_bias > 1.0:
        # Refine near boundaries based on drivers
        for drv in drivers:
            driver_face = drv['driver_face']
            edge = find_shared_edge(target_face, driver_face)
            if edge is None:
                continue
            
            along_v, boundary_index = detect_edge_direction_and_boundary(bs, edge)
            target_poles = min_u if along_v else min_v
            modified |= insert_boundary_biased_knots(bs, not along_v, boundary_index, bias=boundary_bias, target_poles=target_poles, tol=tol)
    else:
        nu, nv = bs.NbUPoles, bs.NbVPoles
        if nu < min_u:
            modified |= insert_uniform_knots(bs, along_v=False, target_poles=min_u, tol=tol)

        if nv < min_v:
            modified |= insert_uniform_knots(bs, along_v=True, target_poles=min_v, tol=tol)

    ''' Debug: visualize modified surface'''
    '''
    if modified:
        obj = App.ActiveDocument.addObject("Part::Feature", "ModifiedSurface")
        obj.Shape = Part.Face(bs.copy())
    '''
    return bs

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
        row_len: Number of poles in the target row

    Returns:
        tangent_at_pole_dir: list of App.Vector (unit) or None
    """
    u0, u1 = edge.ParameterRange
    surface = driver_face.Surface

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

def compute_g1_row(orig_poles, driver_tangents_dir, along_v, boundary_index, beta,
    desired_ref_dirs=None, desired_ref_weights=None):
    """
    Compute new pole row enforcing G1 continuity.

    Args:
        orig_poles: 2D list [nu][nv] of App.Vector
        driver_tangents_dir: list of App.Vector or None (per pole)
        along_v: True if boundary is U=const, False if V=const
        boundary_index: index of boundary row/column (0-based)
        beta: magnitude blending factor [0..1]
        desired_transverse_dirs: optional 2D list [nu][nv] of App.Vector (unit) or None
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

        # Store desired transverse direction at boundary pole
        if along_v:
            i, jj = adj_row, j   # G¹ affects the adjacent row
        else:
            i, jj = j, adj_row

        w = beta  # or driver weight

        if desired_ref_dirs[i][jj] is None:
            desired_ref_dirs[i][jj] = App.Vector(Tblend)
            desired_ref_weights[i][jj] = w
        else:
            desired_ref_dirs[i][jj] = (
                desired_ref_dirs[i][jj].add(App.Vector(Tblend).multiply(w))
            )
            desired_ref_weights[i][jj] += w

        # Target control net spacing for magnitude
        d = d_orig

        Tfinal = Tblend.multiply(d)
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
        refinement: dict with optional keys 'enabled', 'min_poles_u', 'min_poles_v', 'boundary_bias'

    Returns:
        True if all edges processed, False if some failed.
    """
    # Get BSpline surface from target face   
    if not isinstance(target_face.Surface, Part.BSplineSurface):
        raise RuntimeError("Target face is not a BSpline surface")

    bs = None
    if refinement is not None and refinement.get("enabled", False): # refinement requested, but not guaranteed
        min_u = refinement.get("min_poles_u", 8)
        min_v = refinement.get("min_poles_v", 8)
        boundary_bias = refinement.get("boundary_bias", 1.0)
        tol = refinement.get("tol", 1e-9)

        bs = ensure_min_dof(target_face, drivers, min_u=min_u, min_v=min_v, boundary_bias=boundary_bias, tol=tol)
    else:
        bs = target_face.Surface

    nu, nv = bs.NbUPoles, bs.NbVPoles
    is_rational = bs.isURational() or bs.isVRational()

    # Save original poles and weights
    orig_poles = [[bs.getPole(i+1, j+1) for j in range(nv)] for i in range(nu)]
    orig_weights = None
    if is_rational:
        orig_weights = [[bs.getWeight(i+1, j+1) for j in range(nv)] for i in range(nu)]

    
    desired_ref_dirs = [[None for _ in range(nv)] for _ in range(nu)]
    desired_ref_weights = [[0.0 for _ in range(nv)] for _ in range(nu)]

    def get_ref_dir(i, j):
        d = desired_ref_dirs[i][j]
        w = desired_ref_weights[i][j]

        if d is None or w < 1e-9:
            return None

        v = App.Vector(d)
        if v.Length < 1e-9:
            return None

        v.normalize()
        return v
    
    # Initialize accumulators for pole deltas and weights
    delta_parallel = [[None for _ in range(nv)] for _ in range(nu)] # Store parallel deltas, keep strongest
    delta_orthogonal = [[[] for _ in range(nv)] for _ in range(nu)] # Store orthogonal deltas
    delta_avg = [[[] for _ in range(nv)] for _ in range(nu)] # Store full deltas

    # For weights, if rational, accumulate similarly (optional)
    weights = None
    if is_rational:
        weights = [[[] for _ in range(nv)] for _ in range(nu)]

    # Helper: set pole delta
    def add_delta(i, j, delta,weight=None):
        if delta.Length < 1e-9:
            return
        
        delta_avg[i][j].append(delta)
        if is_rational and weight is not None:
            weights[i][j].append(weight)

        ref_dir = get_ref_dir(i, j)
        if ref_dir is None or ref_dir.Length < 1e-9:
            # Fallback: accumulation
            delta_orthogonal[i][j].append(delta)
            return
        ref = App.Vector(ref_dir)
        ref.normalize()
    
        # Decompose delta
        d_par = App.Vector(ref).multiply(delta.dot(ref))
        d_orth = delta.sub(d_par)

        # --- PARALLEL: keep strongest ---
        cur = delta_parallel[i][j]
        if cur is None or d_par.Length > cur.Length:
            delta_parallel[i][j] = d_par

        # --- ORTHOGONAL: store all ---
        delta_orthogonal[i][j].append(d_orth)
    
    # Process each driver
    for drv in drivers:
        driver_face = drv['driver_face']
        samples = drv.get('samples', 30)
        beta = drv.get('beta', 1.0)
        spread_rows = drv.get('spread_rows', 4)
        falloff = drv.get('falloff', 0.6)

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
        g1_row = compute_g1_row(orig_poles, driver_tangents_dir, along_v, boundary_index, beta, desired_ref_dirs, desired_ref_weights)

        ''' Debug: visualize G1 row '''
        # draw_points(g1_row, color=(1.0,0.0,0.0), size=0.5, name="G1_Row_Points")

        # Spread corrections over multiple rows
        rows_indices = get_spread_rows(boundary_index, spread_rows, nu, nv, along_v)

        # rows_indices = [adjacent_row, next_row, next_row, ...]        
        adj_row = rows_indices[0]

        # --- APPLY G1 ROW TO ADJACENT ROW ---
        for j in range(row_len):
            if along_v:
                i, jj = adj_row, j
                P_orig = orig_poles[i][jj]
            else:
                i, jj = j, adj_row
                P_orig = orig_poles[i][jj]

            delta = App.Vector(g1_row[j]).sub(P_orig)

            if is_rational:
                add_delta(i, jj, delta, orig_weights[i][jj])
            else:
                add_delta(i, jj, delta)

        for idx_row, r in enumerate(rows_indices):
            if r == adj_row:
                continue  # adjacent row already handled by G1 row

            k = idx_row  # distance from adjacent row
            w = falloff ** k

            for j in range(row_len):

                # Get original poles at target row and adjacent row
                if along_v:
                    P_orig = orig_poles[r][j]
                    P_adj  = orig_poles[adj_row][j]
                else:
                    P_orig = orig_poles[j][r]
                    P_adj  = orig_poles[j][adj_row]

                # Compute delta from adjacent row original to G1 adjusted row
                delta = App.Vector(g1_row[j]).sub(P_adj)

                # Scale by falloff weight for this row distance
                applied_delta = delta.multiply(w)

                # Map indices correctly for delta accumulation depending on along_v
                i, jj = (r, j) if along_v else (j, r)

                # Rational: keep original weights (no blending here)
                if is_rational:
                    add_delta(i, jj, applied_delta, orig_weights[i][jj])
                else:
                    add_delta(i, jj, applied_delta)

        
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
        
        # Compute final poles by averaging deltas
        if collision_mode == 'average':
            for i in range(nu):
                for j in range(nv):
                    if len(delta_avg[i][j]) > 0:
                        new_pole = orig_poles[i][j].add(v_avg(delta_avg[i][j]))
                        bs.setPole(i + 1, j + 1, new_pole)
                        if is_rational:
                            avg_weight = orig_weights[i][j]
                            if weights and len(weights[i][j]) > 0:
                                avg_weight = w_avg(weights[i][j])
                            bs.setWeight(i + 1, j + 1, avg_weight)
        # Compute final poles by derectionaly averaging deltas
        elif collision_mode == 'directional_clamp':
            for i in range(nu):
                for j in range(nv):
                    d = App.Vector(0,0,0)

                    if delta_parallel[i][j]:
                        d = d.add(delta_parallel[i][j])

                    if len(delta_orthogonal[i][j]) > 0:
                        d = d.add(v_avg(delta_orthogonal[i][j]))

                    if d.Length > 1e-9:
                        new_pole = orig_poles[i][j].add(d)
                        bs.setPole(i + 1, j + 1, new_pole)

                        if is_rational:
                            avg_weight = orig_weights[i][j]
                            if weights and len(weights[i][j]) > 0:
                                avg_weight = w_avg(weights[i][j])
                            bs.setWeight(i + 1, j + 1, avg_weight)
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
        'samples': 30,      # number of samples along edge
        'beta': 1.0,        # blending factor (0.0 - 1.0). Higher = stronger driver influence.
        'spread_rows': 5,   # number of rows to spread corrections over. Higher = smoother transition, but spreading will not go past half the surface.
        'falloff': 0.6      # falloff factor for spreading (0.0 - 1.0). Lower = faster falloff.
    })

# Optional: refinement before G1 enforcement
g1_refinement = {
    "enabled": True,
    "min_poles_u": 9,
    "min_poles_v": 12,
    "boundary_bias": 1.2,
    "tol": 1e-9
}

#bs = enforce_G1_multiple(target_face, drivers, collision_mode='average') # more gentle corners blending
bs = enforce_G1_multiple(target_face, drivers, collision_mode='directional_clamp', refinement=g1_refinement) # better G1 in tight corners
new_face = Part.Face(bs)

doc = App.ActiveDocument
obj = doc.addObject("Part::Feature", "G1_ConstrainedSurface")
obj.Shape = new_face
doc.recompute()