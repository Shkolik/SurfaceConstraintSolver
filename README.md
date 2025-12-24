# SurfaceConstraintSolver
POC for surface-to-surface tangency solver

Enforces G1 (or near G1) continuity between target surface and 1 or more driver surfaces. Driver surfaces should share boundary edge with target surface.

In FreeCAD select target face, then CTRL+select all driver faces. Run macro.
As result new Feature will be added to the TreeView with a shape of modified surface.

For each driver surface following settings exposed (has to be changed in a code):
- samples: number of samples along edge - used for sampling tangents.
- beta: blending factor (0.0 - 1.0). Higher = stronger driver influence. 1.0 - exact or as close as possible to G1
- spread_rows: number of rows to spread corrections over. Help not to kink or distort surface 
- falloff: falloff factor for spreading (0.0 - 1.0). Lower = faster falloff.

If target surface is low on DOF G1 may be impossible. To help with that case there is additional setting "g1_refinement":
- enabled: should be True to allow relaxing surface DOF
- min_poles_u: depends on curvature, but good starting point is 8
- min_poles_v: depends on curvature, but good starting point is 8
- boundary_bias: (>=1.0) if more than 1, algorithm will insert additional knots closer to the boundary edge, else knots will be inserted uniformely. Good starting point is about 1.1-1.2 
- tol: by default 1e-9. Optional. Used in knot insertion.

