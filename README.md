# Turbulent-Bump-NASA-TMR-test-case

steady Turbulent bump in channel problem from NASA turbulence modelling resource with Spalart - Allmaras turbulence model 


Instructions 

1.- Download all files.

1.- Open TurbulentBump.m

2.- Select the mesh or generate a new one with meshbumpX_new_modified.m

3.- Run TurbulentBump.m and save the solution

4.- Ramp the solution towards a higher inlet velocity and let the solution converge.

5.- Interpolate to a finer mesh using InterpolationBetweenMeshes.m and save the interpolated solution.

6.- Disable all blocks in TurbulentBump.m till main solver , run testContinuity.m to pre-proccess the interpolated solution.

7.- Run TurbulentBump.m and repeat the process until desired Re (Maximum Re Tested is 3e6).
