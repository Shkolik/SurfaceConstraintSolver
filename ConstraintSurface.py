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
from enum import Enum

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


#### Enums ####
class Direction(Enum):
    U = 1
    V = 2

class DecayType(Enum):
    SMOOTH = 1
    LINEAR = 2    
    GAUSSIAN = 3

class RefineType(Enum):
    UNIFORM = 1
    WEIGHTED = 2  

class CollisionMode(Enum):
    AVERAGE = 1
    DIR_CLAMP = 2  

class ConstraintMode(Enum):
    SINGLE = 1
    PLANE = 2  
    CLAMP = 3
#### End Enums ####


#### Utils ####
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

def linear_decay(d, max_d):
    return max(0.0, 1.0 - d / max_d)

def smooth_decay(d, max_d):
    t = min(1.0, d / max_d)
    return 1.0 - (3*t*t - 2*t*t*t)

def gaussian_decay(d, sigma):
    return math.exp(-(d*d) / (2*sigma*sigma))
#### End utils ####


#### Core functionality ####

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

def best_possible_transverse(surface, point, edge_tangent, driver_dir):
    """
    Returns the best achievable transverse direction for G1 continuity.
    Args:
        surface: Part.Surface
        point: App.Vector
        edge_tangent: App.Vector (unit)
        driver_dir: App.Vector (unit)
    Returns:
        transverse tangent App.Vector (unit) or None
    """    
    
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

def find_transverse_tangent(surface, point, edge_tangent):
    """    
    Compute surface tangent transverse to the edge using projection onto basis.

    Args:
        surface: (Part.Surface) surface to analyze
        point: App.Vector
        edge_tangent: App.Vector (unit)

    Returns:
        transverse tangent App.Vector (unit) or None
    """
    u, v = surface.parameter(point)
    N = surface.normal(u, v)

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

    d_dir = best_possible_transverse(surface, point, edge_tangent, driver_dir)

    return d_dir

def get_surface_boundary_tangents(surface, driver, pole_params = None):
    """
    Get transverse tangents at boundary poles.

    Args:
        surface: Part.Surface
        driver: (Driver) driver surface data
        pole_params: list of float (surface parameters u or v of boundary poles)
    Returns:
        list of App.Vector (unit) per boundary pole
    """
    # Sample driver tangents along edge
    sample_dirs = []

    for t in driver.sample_params:
        p = driver.edge.Curve.value(t)

        # Edge tangent
        Te = driver.edge.Curve.tangent(t)[0]
        if Te.Length < 1e-9:
            continue
        Te.normalize()

        # Transverse tangent on surface
        Td = find_transverse_tangent(surface, p, Te)
        if Td is None or Td.Length < 1e-9:
            continue

        Td.normalize()

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
        pe = max(min(driver.sample_params), min(max(driver.sample_params), p))

        # Find bracketing samples
        for k in range(len(driver.sample_params) - 1):
            t0 = driver.sample_params[k]
            t1 = driver.sample_params[k + 1]

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
                pole_locs.append(driver.edge.Curve.value(pe))
                pole_dirs.append(Td)
                break
        else:
            # Outside range (numerical edge case)
            if pe < driver.sample_params[0]:
                pole_dirs.append(App.Vector(sample_dirs[0]))
            else:
                pole_dirs.append(App.Vector(sample_dirs[-1]))
    return pole_dirs


class HotRegion:
    def __init__(self, u0, u1, stiffness_ratio, severity, curvature_energy):
        self.u0 = u0
        self.u1 = u1
        self.stiffness_ratio = stiffness_ratio
        self.severity = severity
        self.curvature_energy = curvature_energy

    def collect_spans(self, knots, periodic=False):
        spans = []
        def collect_linear(u0, u1):
            for i in range(len(knots) - 1):
                a, b = knots[i], knots[i + 1]
                span_len = b - a

                if span_len <= 1e-9:
                    continue

                if b <= u0 or a >= u1:
                    continue

                overlap = min(b, u1) - max(a, u0)
                if overlap <= 0:
                    continue

                coverage = overlap / span_len
                spans.append((i, coverage, span_len))

        if periodic and self.u1 < self.u0:
            # wrapped region: [u0, max] U [min, u1]
            umin, umax = knots[0], knots[-1]
            collect_linear(self.u0, umax)
            collect_linear(umin, self.u1)
        else:
            collect_linear(self.u0, self.u1)

        return spans

    def estimate_parallel_layers(self, max_layers=8):
        """
        Estimate number of parallel knots needed from region severity.
        """

        # severity ≈ how many "curvature quanta"
        layers = int(math.ceil(self.severity))

        return min(max(1, layers), max_layers)

class Refinement:
    def __init__(self, type, min_poles_u = 8, min_poles_v = 8, knot_tol=0.05, energy_tol=0.15, angle_tol=0.2):
        self.type = type
        self.min_poles_u = min_poles_u
        self.min_poles_v = min_poles_v
        self.knot_tol = knot_tol
        self.energy_tol = energy_tol
        self.angle_tol = angle_tol

class Driver:
    def __init__(self, face, samples = 30, beta=1.0, decay_type=DecayType.SMOOTH, fairing_distance=0.6, 
                 gaussian_sigma=2.0, energy_tol=0.15, angle_tol=0.04):
        self.face = face
        self.samples = samples
        self.edge = None
        self.boundary_direction = None
        self.boundary_side = None
        self.beta = beta
        self.decay_type = decay_type
        self.fairing_distance = fairing_distance
        self.driver_dirs = []                       # driver transverse tangents along boundary (sampled)
        self.target_dirs = []                       # target transverse tangents along boundary (sampled)
        self.sample_params = []                     # edge parameters per sample
        self.gaussian_sigma = gaussian_sigma        # for DecayType.GAUSSIAN
        self.energy_tol = energy_tol             
        self.angle_tol = angle_tol
        self.angles_mismatch = []                   # angle mismatches along boundary (radians) per sample
        self.mismatch_vectors = []                  # mismatch vectors per sample
        self.min_stiffness_ratio = 2.0              # minimum stiffness ratio to consider hot region
        self.spread_rows = []                       # rows to spread G1 influence into
        self.decay = None                           # decay function
        self.curvature_energy = 0.0                 # geometry-scaled curvature energy along boundary

    def _compute_curvature_energy(self):
        """
        Geometry-scaled curvature energy along boundary.
        """

        if self.samples < 2:
            return 0.0

        E = 0.0
        for i in range(len(self.sample_params) - 1):
            du = self.sample_params[i+1] - self.sample_params[i]
            a0 = self.angles_mismatch[i]
            a1 = self.angles_mismatch[i+1]

            # trapezoidal integration of squared angle
            E += 0.5 * (a0*a0 + a1*a1) * du

        # normalize by tolerance ~ 2-3 degrees
        E /= (self.angle_tol * self.angle_tol)
        return E
    
    def _compute_tangents(self, bs):
        """
        Compute driver and target transverse tangents along the shared edge.

        Args:
            bs: target BSpline surface
        Sets:
            self.driver_dirs: list of App.Vector (unit) per sample along boundary
            self.target_dirs: list of App.Vector (unit) per sample along boundary
            self.angles_mismatch: list of float (radians) per sample along boundary
            self.mismatch_vectors: list of App.Vector per sample along boundary
        """
        if not self.edge:
            raise RuntimeError("Shared edge not set for driver")

        u0, u1 = self.edge.ParameterRange
        eps = (u1 - u0) * 1e-9
    
        params = []
        for i in range(self.samples):
            if i == 0:
                t = u0 + eps
            elif i == self.samples - 1:
                t = u1 - eps
            else:
                t = u0 + (u1 - u0) * i / (self.samples - 1)
            params.append(t)
        
        self.sample_params = params
        # Sample driver and target tangents along edge
        self.driver_dirs = get_surface_boundary_tangents(self.face.Surface, self)
        self.target_dirs = [d.negative() for d in get_surface_boundary_tangents(bs, self)]

        draw_vectors([self.edge.valueAt(p) for p in self.sample_params], self.driver_dirs, color=(1.0, 0.0, 0.0), name="DriverDirs", show=True)
        draw_vectors([self.edge.valueAt(p) for p in self.sample_params], self.target_dirs, color=(0.0, 0.0, 1.0), name="TargetDirs", show=True)
        # Compute angles mismatch
        self.angles_mismatch = [
            math.acos(max(-1.0, min(1.0, abs(d.dot(t)))))
            for d, t in zip(self.driver_dirs, self.target_dirs)
        ]

        # Compute mismatch vectors        
        self.mismatch_vectors = [
            App.Vector(d).sub(t)
            for d, t in zip(self.driver_dirs, self.target_dirs)
        ]

        self.curvature_energy = self._compute_curvature_energy()

    def get_normalized_weights(self):
            weights_mismatch = [max(vec.Length, 1e-7) for vec in self.mismatch_vectors]
            s = sum(weights_mismatch)
            if s < 1e-12:
                return [1.0 / len(self.mismatch_vectors)] * len(self.mismatch_vectors)
            return [x / s for x in weights_mismatch]
    
    def set_shared_edge(self, bs, edge):
        """
        Detect if edge runs along U or V direction of the BSpline surface,
        and determine which boundary row/column it corresponds to.

        Args:
            bs: BSpline surface
            edge: Edge to analyze  

        Sets:
            self.boundary_direction: Direction of the boundary (Direction.U or Direction.V)
            self.boundary_side: Side of the boundary (1 for first row/column, -1 for last)
            self.edge: The edge being analyzed   
            self.driver_dirs: list of App.Vector (unit) per sample along boundary
            self.target_dirs: list of App.Vector (unit) per sample along boundary   
            self.sample_params: list of float edge parameters per sample
            self.angles_mismatch: list of float (radians) per sample along boundary
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

        if du < dv:
            # Edge runs along V, boundary is U = const
            direction = Direction.V

            if abs(us[0] - u_min) < abs(us[0] - u_max):
                side = 1      # first U row
            else:
                side = -1     # last U row
        else:
            # Edge runs along U, boundary is V = const
            direction = Direction.U

            if abs(vs[0] - v_min) < abs(vs[0] - v_max):
                side = 1      # first V row
            else:
                side = -1     # last V row

        self.boundary_direction = direction
        self.boundary_side = side
        self.edge = edge

        self._compute_tangents(bs)

    def estimate_transverse_refinement_need(self, bs, max_layers=8):
        """
        Estimate transverse refinement layers from predicted G1 energy
        Args:
            drv: Driver surface data
            bs: target BSpline surface
            max_layers: maximum allowed refinement layers
        """
        degree = bs.VDegree if self.boundary_direction == Direction.U else bs.UDegree
        
        if degree <= 0 or not self.driver_dirs:
            return 0

        max_mismatch = max([vec.Length for vec in self.mismatch_vectors])

        # predicted dimensionless energy
        E = self.beta * max_mismatch / degree

        print(f"Predicted transverse energy: {E:.6f}. Energy tolerance: {self.energy_tol:.6f}")

        if E < self.energy_tol:
            return 0
        
        if E > 6.0:
            App.Console.PrintWarning(f"High predicted G1 energy={E:.2f}. G1 may not be possible.\n")
            return max_layers

        layers = int(math.ceil(math.sqrt(E / self.energy_tol)))
        return max(0, min(layers, max_layers))
    
    # def localize_hot_region(self, bs):
    #     """
    #     Decide whether localized parallel refinement is needed.
    #     And calculate hot region data if so.
    #     Args:
    #         bs: target BSpline surface

    #     Returns:
    #         hot_region data or None
    #     """
    #     sample_params = []

    #     if self.boundary_direction == Direction.V:
    #         # U = const boundary, iterate along V
    #         v_min, v_max = bs.bounds()[2:4]

    #         for k in range(self.samples):
    #             s = k / (self.samples - 1)
    #             v = v_min + s * (v_max - v_min)
    #             sample_params.append(v)
    #     else:
    #         # V = const boundary, iterate along U
    #         u_min, u_max = bs.bounds()[0:2]

    #         for k in range(self.samples):
    #             s = k / (self.samples - 1)
    #             u = u_min + s * (u_max - u_min)
    #             sample_params.append(u)

    #     # Detect hot region
    #     hot_idx = [i for i,a in enumerate(self.angles_mismatch) if a >= self.angle_tol]
        
    #     region = None

    #     if hot_idx:
    #         i0, i1 = hot_idx[0], hot_idx[-1]
    #         u0, u1 = sample_params[i0], sample_params[i1]
    #         hot_width = max(u1 - u0, 1e-12)

    #         # Effective tangential stiffness proxy:
    #         # assume first interior row "sees" one knot span
    #         # we normalize by boundary param length = 1.0
    #         stiffness_ratio = 1.0 / hot_width

    #         if stiffness_ratio < self.min_stiffness_ratio:
    #             return None

    #         severity = max(self.angles_mismatch[i0:i1+1]) / self.angle_tol

    #         energy = 0.0

    #         for i in range(i0, i1):
    #             du = sample_params[i+1] - sample_params[i]
    #             a0 = self.angles_mismatch[i]
    #             a1 = self.angles_mismatch[i+1]

    #             # Trapezoidal integration of squared angle
    #             energy += 0.5 * (a0*a0 + a1*a1) * du

    #         # Normalize by tolerance and region width
    #         energy /= (self.angle_tol * self.angle_tol * hot_width)
            
    #         region = HotRegion(u0, u1, stiffness_ratio, severity, energy)

    #     return region
    
    
    def bending_capacity(self, bs):
        """
        Count number of spans along boundary direction.
        Args:
            bs: (Part.BSplineSurface) target surface
        Returns:
            int: number of spans along boundary
        """
        spans = 0
        degree = 0
        if self.boundary_direction == Direction.V:
            spans = bs.NbVKnots - (0 if bs.isVPeriodic() else 1)
            degree = bs.VDegree
        else:
            spans = bs.NbUKnots - (0 if bs.isUPeriodic() else 1)
            degree = bs.UDegree
        
        return spans * degree

    def compute_bending_energy(self, degree_parallel=1):
        """
        Measures total bending demand along boundary.
        Drives parallel refinement.
        """
        if self.samples < 2:
            return 0.0
        
        ds = self.edge.Length / (self.samples - 1)

        E = 0.0
        for v in self.mismatch_vectors:
            E += (v.Length ** 2) * ds

        return E / max(degree_parallel, 1)

    def estimate_parallel_refinement_need(self, bs):
        """
        Global bending demand estimator.
        """
        degree = bs.UDegree if self.boundary_direction == Direction.U else bs.VDegree
        V = 0.0
        for i in range(len(self.angles_mismatch)-1):
            da = abs(self.angles_mismatch[i+1] - self.angles_mismatch[i])
            ds = self.sample_params[i+1] - self.sample_params[i]
            if ds > 1e-9:
                V = max(V, da / ds)

        required_spans = math.ceil(self.edge.Length * V / (degree + 1))
        current_spans = bs.NbUKnots - (0 if bs.isUPeriodic() else 1) if self.boundary_direction == Direction.U else bs.NbVKnots - (0 if bs.isVPeriodic() else 1)
        return max(0, required_spans - current_spans)
    
    def estimate_surface_refinement(self, bs):
        """
        Estimates required transverse + parallel refinement layers and degree elevation
        to achieve G1 continuity along shared edge between driver and target surfaces.
        Args:
            bs: (Part.BSplineSurface) target surface
        Returns:
            estimation: (RefinementEstimation) estimation result
        """

        # Transverse refinement from tangent mismatch
        transverse_layers = self.estimate_transverse_refinement_need(bs)

        print(f"Pre-estimated {transverse_layers} transverse layers from tangent mismatch")

        # # Compute "hot" region where angle mismatch is significant
        # # That will give us region width, stiffness ratio, severity and curvature energy
        # hot_region = self.localize_hot_region(bs)

        parallel_layers = self.estimate_parallel_refinement_need(bs)
        print(f"Pre-estimated {parallel_layers} parallel layers from bending demand")

        transverse_layers = int(math.ceil(transverse_layers * math.sqrt(1.0 + self.curvature_energy)))
        print(f"Adjusted to {transverse_layers} transverse layers after curvature energy scaling")

        return RefinementEstimation(
            transverse_layers=transverse_layers,
            parallel_layers=parallel_layers,
            curvature_energy=self.curvature_energy
            )
    
    def boundary_pole_params(self, bs):
        """
        Return surface parameters corresponding to boundary poles.

        Args:
            bs (Part.BSplineSurface): target surface
        Returns:
            list[float]: parameter values along boundary direction
                        (U or V depending on along_v)
        """

        nu = bs.NbUPoles
        nv = bs.NbVPoles

        params = []

        if self.boundary_direction == Direction.V:
            # U = const boundary, iterate along V
            i = 1 if self.boundary_side > 0 else nu
            for j in range(1, nv + 1):
                P = bs.getPole(i, j)
                u, v = bs.parameter(P)
                params.append(v)
        else:
            # V = const boundary, iterate along U
            j = 1 if self.boundary_side > 0 else nv
            for i in range(1, nu + 1):
                P = bs.getPole(i, j)
                u, v = bs.parameter(P)
                params.append(u)

        return params

class MoveIntent:
    def __init__(self, direction, magnitude, driver_index, boundary_index, boundary_direction, priority=0):
        self.direction = direction
        self.magnitude = magnitude
        self.priority = priority
        self.driver_index = driver_index
        self.boundary_index = boundary_index
        self.boundary_direction = boundary_direction

class PoleConstraint:
    def __init__(self, mode, direction, max_mag, priority = 0, basis = None):
        """
        Args:
            mode: (ConstraintMode) mode
            basis: list of App.Vector (unit) spanning constraint directions
            direction: App.Vector (unit) representing the constraint direction
            max_mag: float, maximum allowed movement magnitude
            priority: int, higher priority constraints wins
        """
        self.basis = basis
        self.mode = mode
        self.direction = direction
        self.max_mag = max_mag
        self.priority = priority

class PoleDelta:
    def __init__(self):        
        self.raw = []
        self.parallel = None
        self.orthogonal = []        
        self.weight = []

class RefinementEstimation:
    def __init__(self, transverse_layers:int=0, parallel_layers:int=0, hot_region:HotRegion=None, 
                 bending_energy:float=0.0, curvature_energy:float=0.0):
        self.transverse_layers = transverse_layers
        self.parallel_layers = parallel_layers
        self.hot_region = hot_region
        self.bending_energy = bending_energy
        self.curvature_energy = curvature_energy

class RefinementPlan:
    def __init__(self):
        self.U = [] # list of U knots to insert
        self.V = [] # list of V knots to insert

    def add_parallel_knots(self, estimation:RefinementEstimation, driver:Driver):        
        if estimation.parallel_layers <= 0:
            return

        if estimation.hot_region:
            u0, u1 = estimation.hot_region.u0, estimation.hot_region.u1
        else:
            u0, u1 = 0.0, 1.0

        for i in range(estimation.parallel_layers):
            knot = u0 + (i + 1) / (estimation.parallel_layers + 1) * (u1 - u0)
            if driver.boundary_direction == Direction.V:
                self.V.append(knot)
            else:
                self.U.append(knot)

    def add_transverse_knots(self, bs, estimation:RefinementEstimation, driver:Driver):
        """
        Propose boundary-layer knots near a given boundary.
        """
        if estimation.transverse_layers <= 0:
            return
        
        # Direction setup
        if driver.boundary_direction == Direction.V:
            add = lambda t: self.U.append(t)
            surface_extent = bs.bounds()[1] - bs.bounds()[0]
            knots = bs.getUKnots()
            degree = bs.UDegree
            param_min, param_max = bs.bounds()[0:2]
        else:
            add = lambda t: self.V.append(t)
            surface_extent = bs.bounds()[3] - bs.bounds()[2]
            knots = bs.getVKnots()
            degree = bs.VDegree
            param_min, param_max = bs.bounds()[2:4]

        boundary_param = param_min if driver.boundary_side > 0 else param_max
        at_start = (driver.boundary_side > 0)

        # hot_width = estimation.hot_region.u1 - estimation.hot_region.u0 if estimation.hot_region else 1.0
        severity = estimation.hot_region.severity if estimation.hot_region else 1.0
        p = 1.0 + 0.5 * min(severity - 1.0, 2.0)
        
        def degree_adaptive_alpha(degree):
            """
            Boundary safety coefficient.
            Increases gently with degree.
            """
            # Clamp degree to a reasonable range
            d = max(1, min(degree, 9))
            return 0.25 + 0.02 * d

        def compute_min_depth(knots, degree, surface_extent, geometric_frac=0.02):
            """
            Compute safe minimum transverse knot depth.
            """

            # --- geometric fallback (always valid)
            geom_min = geometric_frac * surface_extent

            # --- find interior spans
            interior_spans = [
                knots[i+1] - knots[i]
                for i in range(len(knots) - 1)
                if knots[i] > 0.0 and knots[i+1] < 1.0
            ]

            if not interior_spans:
                # Only [0, 1] → fallback
                return geom_min

            span = min(interior_spans)
            alpha = degree_adaptive_alpha(degree)

            return max(geom_min, alpha * span)

        min_depth = compute_min_depth(knots, degree, surface_extent, 0.02)
        
        if estimation.hot_region:
            hot_width = estimation.hot_region.u1 - estimation.hot_region.u0
        else:
            hot_width = min(0.15 * surface_extent, 3.0 * min_depth)

        depths = []
        for l in range(estimation.transverse_layers):
            t = (l + 1) / (estimation.transverse_layers + 1)
            d = hot_width * (t ** p)
            d = max(d, min_depth)
            depths.append(d)

        eps = 1e-9 * surface_extent
        for d in depths:
            u = boundary_param + d if at_start else boundary_param - d
            u = max(param_min + eps, min(param_max - eps, u))
            add(u)

    
    def add_transverse_mismatch_weighted_knots(self, bs, driver:Driver, layers:int, tol=1e-9):
        """
        Insert transverse knots near boundary weighted by tangent mismatch.

        Args:
            bs (Part.BSplineSurface): target surface
            driver (Driver): driver surface data
            layers (int): number of layers to insert
            tol (float) : knot insertion tolerance
        """

        if not driver or not layers:
            return

        samples = min(len(driver.driver_dirs), len(driver.target_dirs))
        

        # Direction setup
        if driver.boundary_direction == Direction.V:
            add = lambda t: self.U.append(t)
            param_min, param_max = bs.bounds()[0:2]
        else:
            add = lambda t: self.V.append(t)
            param_min, param_max = bs.bounds()[2:4]

        u0, u1 = driver.edge.ParameterRange
        edge_params = [u0 + (u1 - u0) * i / (samples - 1)
                    for i in range(samples)]

        # Normalized mismatch weights
        weights = driver.get_normalized_weights()

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

        # Generate insertion parameters
        refine_params = weighted_refine_params(
            edge_params,
            weights,
            layers
        )

        # Map edge → surface parameter
        eps = (param_max - param_min) * 1e-6
        at_start = (driver.boundary_side > 0)
        boundary_param = param_min if at_start else param_max

        modified = False
        for t in refine_params:
            s = (t - u0) / (u1 - u0)
            d = s * (param_max - param_min) * 0.5

            u = boundary_param + d if at_start else boundary_param - d
            u = max(param_min + eps, min(param_max - eps, u))

            add(u)
            modified = True

        return modified

    # may be not needed for now
    @staticmethod
    def insert_or_nudge(u, existing, tol, direction):
        """
        direction: +1 pushes inward (toward 1), -1 pushes inward (toward 0)
        """
        existing = sorted(existing)

        for k in existing:
            if abs(u - k) < tol:
                # nudge inward
                u2 = k + direction * tol
                # clamp
                u2 = max(0.0, min(1.0, u2))

                # verify ordering
                if all(abs(u2 - kk) > tol for kk in existing):
                    return u2
                else:
                    return None

        return u

    def needs_degree_elevation(self, bs, direction):
        return False
        if direction == Direction.U:
            return len(self.U) > (bs.NbUPoles - bs.UDegree - 0 if bs.NbUKnots == 2 else 1)#bs.UDegree * 2
        if direction == Direction.V:
            return len(self.V) > (bs.NbVPoles - bs.VDegree - 0 if bs.NbVKnots == 2 else 1)#bs.VDegree * 2
        
    def resolve_U(self, existing, tol=0.05):
        new = []
        merged = sorted(existing)
        for u in sorted(self.U):
            if all(abs(u - k) > tol for k in merged):
                merged.append(u)
                new.append(u)
            else:
                print(f"U knot {u} skipped due to proximity")        
        return sorted(new)
    
    def resolve_V(self, existing, tol=0.05):
        new = []
        merged = sorted(existing)
        vMin, vMax = merged[0], merged[-1]
        mid = 0.5 * (vMin + vMax)
        for v in sorted(self.V):
            if v < mid:
                d = v-vMin
            else:
                d = vMax - v
            s = max(0.0, min(d / (vMax - vMin), 1.0))
            tol_scaled = 0.01 + (tol - 0.01) * s
            if all(abs(v - k) > tol_scaled for k in merged):
                merged.append(v)
                new.append(v)
            else:
                print(f"V knot {v} skipped due to proximity")        
        return sorted(new)
    
    def resolve(self, existing_u_knots, existing_v_knots, tol=0.05):
        u_new = self.resolve_U(existing_u_knots, tol=tol)
        v_new = self.resolve_V(existing_v_knots, tol=tol)
        return u_new, v_new


def insert_uniform_knots(bs, direction, target_poles=8, tol=1e-9):
    """
    Insert uniform knots into a BSplineSurface to reach at least target pole count.
    Shape is preserved exactly.

    Args:
        bs (Part.BSplineSurface): surface to refine (MODIFIED IN PLACE)
        direction (Direction): Direction to refine (U or V)
        target_poles (int): minimum number of poles in that direction
        tol (float): knot insertion tolerance
    Returns:
        bool: True if the surface was modified, False otherwise
    """
    modified = False

    if direction == Direction.V:
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


#### Old refinement functions ####
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
#### End old refinement functions ####



def compute_g1_row(orig_poles, driver_tangents_dir, driver:Driver):
    """
    Compute new pole row enforcing G1 continuity.

    Args:
        orig_poles: 2D list [nu][nv] of App.Vector
        driver_tangents_dir: list of App.Vector or None (per pole)
        driver: (Driver) driver surface data
    Returns:
        g1_row: list of App.Vector (new adjacent pole row)
    """
    g1_row = []

    row_len = len(driver_tangents_dir)
    step = driver.boundary_side
    direction = driver.boundary_direction
    if driver.boundary_side > 0:
        boundary_index = 0
    else:
        boundary_index = (len(orig_poles) - 1 if direction == Direction.V else len(orig_poles[0]) - 1)

    adj_row = boundary_index + step
    
    for j in range(row_len):
        if direction == Direction.V:
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
        Tblend = Norig.multiply(1.0 - driver.beta).add(T.multiply(driver.beta))
        if Tblend.Length < 1e-6:
            Tblend = T

        Tblend.normalize()

        Tfinal = Tblend.multiply(d_orig)
        # draw_vector(P1, Tfinal, (0.0,0.0,1.0), 10.0, "Tfinal_Blue")

        g1_row.append(P1.add(Tfinal))

    return g1_row

def resolve_constraints(desired_dirs):
    """
    Build per-pole motion constraints from desired directions.

    Args:
        desired_dirs: dict[(i,j)] -> list of (MoveIntent) intents 

    Returns:
        pole_constraints: dict[(i,j)] -> constraint dict
    """

    pole_constraints = {}

    for key, entries in desired_dirs.items():

        if not entries:
            continue

        # Priority filtering
        hard = [e for e in entries if e.priority == 0]

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
            d = App.Vector(e.direction)
            if d.Length < 1e-9:
                continue
            d.normalize()
            dirs.append(d)
            mags.append(abs(e.magnitude))
            priorities.append(e.priority)

        if not dirs:
            continue

        # Single driver
        if len(dirs) == 1:
            pole_constraints[key] = PoleConstraint(
                ConstraintMode.SINGLE, 
                dirs[0], 
                mags[0], 
                priorities[0]
            )
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

        # --- one dominant cluster ---
        if len(clusters) == 1:
            c = clusters[0]
            pole_constraints[key] = PoleConstraint(
                ConstraintMode.SINGLE,
                c["dir"], 
                max(c["mags"]), 
                min(c["priorities"])
            )
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

                pole_constraints[key] = PoleConstraint(
                    ConstraintMode.PLANE,
                    None,
                    min(max(c1["mags"]), max(c2["mags"])),
                    min(min(c1["priorities"]), min(c2["priorities"])),
                    basis=[u, v]
                )
                continue

        # --- fallback: clamp ---
        # Multiple conflicting directions, low DOF
        # Allow movement but cap magnitude and priority
        avg_dir = v_avg(dirs)

        if avg_dir.Length > 1e-9:
            avg_dir.normalize()

        pole_constraints[key] = PoleConstraint(
            ConstraintMode.CLAMP,
            avg_dir,
            min(mags),
            min(priorities)
        )

    return pole_constraints

def apply_constraint(delta, c:PoleConstraint):
    if delta.Length < 1e-9:
        return delta

    max_mag = c.max_mag if c.priority <= 0 else c.max_mag * (1.0 / (0.3 + c.priority))

    if c.mode == ConstraintMode.SINGLE:
        d = App.Vector(c.direction)
        mag = delta.dot(d)
        mag = max(-max_mag, min(max_mag, mag))
        return d.multiply(mag)

    if c.mode == ConstraintMode.PLANE:
        proj = App.Vector(0, 0, 0)
        for b in c.basis:
            proj = proj.add(App.Vector(b).multiply(delta.dot(b)))

        if proj.Length > max_mag:
            proj.normalize()
            proj = proj.multiply(max_mag)
        return proj

    if c.mode == ConstraintMode.CLAMP:
        if delta.Length > max_mag:
            delta = App.Vector(delta)
            delta.normalize()
            delta = delta.multiply(max_mag)
        return delta

    return delta

def build_pole_basis(bs):
    nu, nv = bs.NbUPoles, bs.NbVPoles
    basis = [[None for _ in range(nv)] for _ in range(nu)]

    poles = bs.getPoles()

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


def enforce_G1_multiple(target_face, drivers:list[Driver], collision_mode=CollisionMode.DIR_CLAMP, refinement:Refinement=None):
    """
    Enforce G1 continuity on target_face along multiple driver surfaces.

    Args:
        target_face: Part.Face to modify
        drivers: (list of Driver) driver faces
        collision_mode: (CollisionMode) collision handling mode. CollisionMode.DIR_CLAMP by default.
        refinement: (Refinement) refinement settings or None

    Returns:
        True if all edges processed, False if some failed.
    """
    # Get BSpline surface from target face   
    if not isinstance(target_face.Surface, Part.BSplineSurface):
        raise RuntimeError("Target face is not a BSpline surface")

    if not drivers:
        raise RuntimeError("No drivers provided for G1 enforcement")
    
    # DEBUG: Record the start time
    start_time = time()

    bs = target_face.Surface.copy()
    
    # track if surface modified during preprocessing - for debugging
    refined = False 

    # Ensure minimum degree 3 for G1 enforcement
    vDegree, uDegree = bs.VDegree, bs.UDegree
    if vDegree < 3 or uDegree < 3:
        App.Console.PrintMessage("Elevating target surface degree to minimum 3 for G1 enforcement\n")
        bs.increaseDegree(max(uDegree, 3), max(vDegree, 3))
        refined = True

    def is_single_span(bs, direction):
        kv = bs.getUKnots() if direction == Direction.U else bs.getVKnots()
        return len(kv) == 2

    def split_single_span(bs, direction, target_spans=2):
        # Uniform split
        for i in range(1, target_spans):
            u = i / target_spans
            if direction == Direction.U:
                bs.insertUKnot(u, 1, 1e-9)
            else:
                bs.insertVKnot(u, 1, 1e-9)

    # Preprocess drivers - detect directions and boundary indices
    for driver_idx, drv in enumerate(drivers):
        # Find shared edge
        edge = find_shared_edge(target_face, drv.face)
        if edge is None:
            App.Console.PrintWarning("Failed to find shared edge with driver face\n")
            continue

        # Detect direction and boundary index 
        drv.set_shared_edge(bs, edge)

        # transverse_dir = Direction.V if drv.boundary_direction == Direction.U else Direction.U
        # # For localized refinement, ensure target surface has >1 span in parallel and transverse direction
        # # fix fo single-span surfaces (bezier) that cannot be refined locally
        # if refinement and is_single_span(bs, drv.boundary_direction):
        #     App.Console.PrintMessage(f"Inserting uniform knots in {drv.boundary_direction.name} for single-span target surface before localized refinement\n")            
        #     split_single_span(bs, drv.boundary_direction, target_spans=2)
        #     refined = True
        # if refinement and is_single_span(bs, transverse_dir):
        #     App.Console.PrintMessage(f"Inserting uniform knots in {transverse_dir.name} for single-span target surface before localized refinement\n")
        #     split_single_span(bs, transverse_dir, target_spans=2)
        #     refined = True
    
    if refinement:        
        if refinement.type == RefineType.UNIFORM:   # for uniform refinement we do it independent of drivers
            min_u = refinement.min_poles_u
            min_v = refinement.min_poles_v
            tol = refinement.tol

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
                refined |= insert_uniform_knots(bs, Direction.U, target_poles=min_u, tol=tol)
            if nv < min_v:
                refined |= insert_uniform_knots(bs, Direction.V, target_poles=min_v, tol=tol)
        else:                                       # localized refinement based on drivers            
            refinement_plan = None
            passNum = 0
            while True:
                passNum += 1
                refinement_plan = RefinementPlan()
                for driver_idx, drv in enumerate(drivers):
                    # print(f"Pass {passNum}, Driver {driver_idx}: Estimating surface refinement\n")
                    estimate = drv.estimate_surface_refinement(bs)

                    # # DEBUG: Print relaxed estimation details
                    print(f"Estimated: {estimate.transverse_layers} transverse layers and {estimate.parallel_layers} parallel layers")

                    refinement_plan.add_parallel_knots(estimate, drv)
                    refinement_plan.add_transverse_knots(bs, estimate, drv)
                    # refinement_plan.add_transverse_mismatch_weighted_knots(bs, drv, estimate.transverse_layers)

                needs_Uelevation = refinement_plan.needs_degree_elevation(bs, Direction.U)
                needs_Velevation = refinement_plan.needs_degree_elevation(bs, Direction.V)

                if not needs_Uelevation and not needs_Velevation:
                    break
                if needs_Uelevation:
                    print(f"Elevating U degree from {bs.UDegree} to {bs.UDegree + 1}\n")
                    bs.increaseDegree(bs.UDegree + 1, bs.VDegree)
                    refined = True
                if needs_Velevation:
                    print(f"Elevating V degree from {bs.VDegree} to {bs.VDegree + 1}\n")
                    bs.increaseDegree(bs.UDegree, bs.VDegree + 1)
                    refined = True

            # Apply knot insertions
            u_new, v_new = refinement_plan.resolve(bs.getUKnots(), bs.getVKnots(), tol=refinement.knot_tol)

            for uK in u_new:
                bs.insertUKnot(uK, 1, 1e-9)
                refined = True

            for vK in v_new:
                bs.insertVKnot(vK, 1, 1e-9)
                refined = True

    # Debug: show modified surface after refinement
    if refined:
        show_surface(bs, "Refined_Surface")

    g1_poles = set()
    # Preprocess drivers - setup spreading and decay functions
    for driver_idx, drv in enumerate(drivers):
        direction = drv.boundary_direction
        side = drv.boundary_side

        # After refinement (if any) boundary index may have changed
        boundary_index = 0 if side == 1 else (bs.NbUPoles - 1 if direction == Direction.V else bs.NbVPoles - 1)
        
        g1_index = boundary_index + side

        if direction == Direction.V:
            for j in range(1, bs.NbVPoles-1):  # skip boundaries
                g1_poles.add((g1_index, j))
        else:
            for i in range(1, bs.NbUPoles-1):  # skip boundaries
                g1_poles.add((i, g1_index))

        spread = max(0, int(drv.fairing_distance * (bs.NbUPoles if direction == Direction.V else bs.NbVPoles)) - 2) # exclude boundaries and g1 row        
        spread_rows = [(boundary_index + side) + delta*side for delta in range(1, max(0, spread) + 1)] # rows to spread into

        if drv.decay_type == DecayType.LINEAR:
            decay = lambda d, md=spread: linear_decay(d, md)
        elif drv.decay_type == DecayType.GAUSSIAN:
            sigma = drv.gaussian_sigma
            decay = lambda d, s=sigma: gaussian_decay(d, s)
        else:
            decay = lambda d, md=spread: smooth_decay(d, md)

        drv.spread_rows = spread_rows                       # rows to spread into
        drv.decay = decay                                   # decay function
        
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
        boundary_index = 0 if drv.boundary_side > 0 else (bs.NbUPoles - 1 if drv.boundary_direction == Direction.V else bs.NbVPoles - 1)
                
        # Get exact boundary pole parameters
        pole_params = drv.boundary_pole_params(bs)

        # Sample driver tangents aligned to those poles
        driver_tangents_dir = get_surface_boundary_tangents(drv.face.Surface, drv, pole_params)

        row_len = bs.NbVPoles if drv.boundary_direction == Direction.V else bs.NbUPoles

        if row_len != len(driver_tangents_dir):
            App.Console.PrintWarning(f"Mismatch row length. tangents length: {len(driver_tangents_dir)}, row length: {row_len}\n")
            continue

        # Compute new G1 row
        g1_row = compute_g1_row(orig_poles, driver_tangents_dir, drv)

        # Next to boundary row        
        adj_row = boundary_index + drv.boundary_side

        # Capture desired transverse directions and intents to move poles
        for j in range(1, row_len-1):  # skip boundaries
            # Get adjacent row pole
            if drv.boundary_direction == Direction.V:
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

            desired_dirs[key].append(MoveIntent(dir, mag, driver_idx, boundary_index, drv.boundary_direction, 0))
            

            # Spread to interior rows
            for r in drv.spread_rows:
                if drv.decay is None:
                    continue
                
                d = abs(r - adj_row)
                d_mag = mag*drv.decay(d)

                if d_mag < 1e-6:
                    continue

                if drv.boundary_direction == Direction.V:
                    key = (r, j)
                else:
                    key = (j, r)
                
                if key in g1_poles:                    
                    continue  # do not override g1 poles

                desired_dirs[key].append(MoveIntent(dir, d_mag, driver_idx, boundary_index, drv.boundary_direction, d))


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
            delta_accum[key] = PoleDelta()

        acc = delta_accum[key]
        acc.raw.append(App.Vector(delta))

        if is_rational and weight is not None:
            acc.weight.append(weight)

        basis = pole_basis[i][j]

        # no basis: pure accumulation
        if not basis:
            acc.orthogonal.append(delta)
            return
        
        # Delta decomposition

        # project into basis
        d_proj = App.Vector(0, 0, 0)
        for b in basis:
            d_proj = d_proj.add(App.Vector(b).multiply(delta.dot(b)))

        d_orth = delta.sub(d_proj)

        # parallel: keep strongest projection
        cur = acc.parallel
        if cur is None or d_proj.Length > cur.Length:
            acc.parallel = d_proj

        # orthogonal: store residuals
        if d_orth.Length > 1e-9:
            acc.orthogonal.append(d_orth)
            
    # Accumulate deltas from all intents
    for key, items in desired_dirs.items():
        i, j = key        
        for item in items:
            delta = App.Vector(item.direction).multiply(item.magnitude)
            if is_rational:
                add_delta(i, j, delta, orig_weights[i][j])
            else:
                add_delta(i, j, delta)

    # Apply accumulated deltas to row poles solving conflicts

    # Compute final poles by averaging deltas
    if collision_mode == CollisionMode.AVERAGE:
        for (i, j), acc in delta_accum.items():            
            new_pole = orig_poles[i][j].add(v_avg(acc.raw))
            bs.setPole(i + 1, j + 1, new_pole)

            if is_rational:
                avg_weight = orig_weights[i][j]
                if acc.weight and len(acc.weight) > 0:
                    avg_weight = w_avg(acc.weight)
                bs.setWeight(i + 1, j + 1, avg_weight)
    # Compute final poles by derectionaly averaging deltas
    elif collision_mode == CollisionMode.DIR_CLAMP:
        for (i, j), acc in delta_accum.items():
            d = App.Vector(0, 0, 0)

            if acc.parallel:
                d = d.add(acc.parallel)

            if acc.orthogonal and len(acc.orthogonal) > 0:
                d = d.add(v_avg(acc.orthogonal))
            if d.Length < 1e-9:
                continue

            bs.setPole(i + 1, j + 1, orig_poles[i][j].add(d))

            if is_rational:
                avg_weight = orig_weights[i][j]
                if acc.weight and len(acc.weight) > 0:
                    avg_weight = w_avg(acc.weight)
                bs.setWeight(i + 1, j + 1, avg_weight)
    
    # draw_points(deltas, (1.0, 0.0, 1.0), 0.2, "Faired_Poles")

    # Record the end time
    end_time = time()

    # Calculate and print the elapsed time
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")

    return bs

# --------------------------------
# USER INPUT
# --------------------------------


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

selection_faces = getSelectedFaces()

if len(selection_faces) < 2:
    raise RuntimeError("Select at least two faces: target face first, driver face(s) second")
target_face = selection_faces[0]
driver_faces = selection_faces[1:]

# Optional: refinement before G1 enforcement
g1_refinement = Refinement(RefineType.WEIGHTED)

drivers = []
for f in driver_faces:
    drivers.append(Driver(
        face=f,
        samples=50,
        beta=1.0,
        decay_type=DecayType.SMOOTH,
        fairing_distance=0.6,
        gaussian_sigma=2.0,
        energy_tol=g1_refinement.energy_tol, # May expose for each driver
        angle_tol=g1_refinement.angle_tol    # May expose for each driver
    ))



#bs = enforce_G1_multiple(target_face, drivers, collision_mode=CollisionMode.AVERAGE) # more gentle corners blending
bs = enforce_G1_multiple(target_face, drivers, CollisionMode.DIR_CLAMP, g1_refinement) # better G1 in tight corners

new_face = Part.Face(bs)

doc = App.ActiveDocument
obj = doc.addObject("Part::Feature", "G1_ConstrainedSurface")
obj.Shape = new_face
doc.recompute()