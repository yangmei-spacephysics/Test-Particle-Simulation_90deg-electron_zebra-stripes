# Test-Particle-Simulation_90deg-electron_zebra-stripes
A 2-D test particle simulation script in MATLAB of tracing electrons in a dipole field in the presence of electrostatic fields.

Assumptions: 
  1. Dipole magnetic field
  2. An ad-hoc electrostatic field pointing dawn to dusk Drift
  3. Electrons are equatorially trapped
  4. The first adiabatic invariant is conserved

Updates 11/16/2025:
  A MATLAB function conducting test particle simulations and sampling simulated electron fluxes on a virtual satellite track is now available. By specifying detailed information of a virtual satellite during a inbound or outbound pass, one can use the function to simulate electron flux spectra being "measured" by the virtual satellite. This function can be used to reproduce Figure 8 in Mei et al.(2025).

Reference: Mei, Y., Li, X., O Brien, D., Xiang, Z., Zhao, H., Sarris, T., et al. (2025). Characteristics of “zebra stripes” of relativistic electrons unveiled by CIRBE/REPTile-2 measurements and test particle simulations. Journal of Geophysical Research: Space Physics, 130, e2024JA033187. https://doi.org/10.1029/2024JA033187

Note: A function "histcnd" is used in the script, which can be found in Mathew (2025) (https://www.mathworks.com/matlabcentral/fileexchange/29435-n-dimensional-histogram-count).
