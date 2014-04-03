Spectrum-GUI
============

A GUI to compute eigenvalue spectra for small finite difference models, motivated by schemes and boundaries used in audio and room acoustics.  The usefulness probably does not extend beyond reproducing results we're trying to publish, but it may be helpful for interested readers.  Please first read the disclaimer at the bottom, and feedback is always welcome.

Usage
-----

I hope the buttons and labels are clear enough to someone looking at the paper, but the basic usage involves setting simulation parameters at the top (e.g. scheme, boundary condition type, reflection coefficient assigned to the boundaries, geometry, problem size, Courant factor, and numerical precision), clicking the button that says "Set up matrix and compute eigs," and then clicking the buttons at the bottom to produce various figures.  After changing any of the simulation parameters, you are prompted to recompute the eigenvalues of the new matrix.

I don't have everything set up for all schemes, so sometimes options or parameters will be disabled.  For instance, with non-box geometries, the interpolated schemes are not available.  And for the interpolated schemes, I don't have velocity-centered boundaries implemented.  The GUI probably does a pretty good job of disabling things that will return incorrect results.  In addition to the parameter options and figure buttons, there is a geometry button which produces a figure showing the geometry of the model.  The "Options" menu is used mostly for plotting parameters.

Simulation Parameters
---------------------

Brief descriptions of simulation parameters which may be changed:

- **Problem Size:** There are three edit fields where you can change the size of a "box" model.  These should be integer dimensions of the 3D box in units of number-of-nodes.  The outer layer of nodes will be boundary nodes.  The number next to the edit boxes is the size of the matrix that will result.  Depending on the computer, it may be wise to keep this number below 1500-2000.  For this, I am only using direct eigenvalue solvers, so it is computing the entire eigen-decomposition.  This can take quite some time.
- **Scheme:** This is the scheme used to update the interior of the model.  I have used the naming conventions in [(Kowalczyk, 2011)](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5440917).  For instance, SLF is the standard rectilinear scheme, IWB is interpolated wideband, and so on.
- **Boundaries:** Boundaries can be either pressure-centered or velocity-centered for the SLF scheme.  For all the others, boundaries are pressure-centered and equivalent to those in [(Kowalczyk, 2011)](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=5440917).
- **Reflection Coefficient:** Boundary absorption is applied uniformly across the domain, and here we specify it by the associated normal incidence pressure reflection coefficient.  Rigid boundaries (homogeneous Neumann) correspond to R = 1, and pressure-release boundaries (homogeneous Dirichlet) correspond to R = -1.  Anything outside this range corresponds to active boundaries, which will be evident in the eigenvalue spectra.
- **Courant factor:** The Courant factor governs time stepping, and it is limited to a maximum value for each scheme.  In this field, the fraction of the maximum stable Courant factor is specified.  For example, 1.0 corresponds to 100% of the maximum value, 0.9 is 90% of the maximum value, and so on.  This is automatically adjusted for each scheme.
- **Numerical Precision:** Construct matrices using either single or double precision coefficients to mimic simulations done in single and double precision.
- **Hard Source Position:** A hard source can be placed at the x-y-z-coordinates in the edit boxes, again integers.  If the boxes are empty, there will be no source node--the case for a soft or transparent source.  The effect is to zero the corresponding row of the matrix.  Soft and transparent sources do not affect the operator.
- **Pre-defined Geometries:** To me, there is no obvious (and sufficiently easy way) to code options for general geometries, so I included a list of pre-defined geometries.  Choosing 'Box' enables many of the previous options.  There are two L-shaped rooms, one seemingly unstable and the other stable.  They differ very slightly in geometry.  The 'Corner' model is a cube with a re-entrant corner.  There are also 'Sphere' and 'Cylinder' models.


Figures
-------

A brief description of figures.  Scales should adjust automatically depending on numerical precision.

- **Eigenvalues:** This is the basic eigenvalue plot in the complex plane.  For the 'Box' models I also plot analytical eigenfrequencies for the ideal problem with either rigid or pressure-release boundaries (this changes depending on whether R is positive or negative).  There is a small threshold to allow for rounding errors, eigenvalues beyond which are colored red.  Once I used a plot like this with inverted colors, so there is an option for this in Options --> Eigenvalue Figures.  Be careful with trusting the exact position of the unit circle; it is just a polygon.

- **Eigenvectors:** This is another GUI-like tool that pops up in a figure window.  When you click near an eigenvalue, it will plot a representation of the real part of the associated eigenvector in a new figure window.  It will find the nearest eigenvalue to your mouse click.  You can zoom and pan, but to be able to click again, de-select the zoom or pan tool.  There are different plot types for the eigenvectors in Options --> Eigenvector Figures.  For the default format, there is another option in Options --> Z-scale, which lets you adjust the relative height of the slices.

- **Eigenvalue Condition Numbers:** This gives you a plot, arranged linear in frequency, of estimated perturbations of each eigenvalue due to rounding errors.  It is based on the eigenvalue condition numbers returned from Matlab's condeig function.  Note the logarithmic scale.
 
- **Run EigTool:** [EigTool](http://www.cs.ox.ac.uk/projects/pseudospectra/eigtool/) is a really neat tool for computing pseudospectra, and this button will load the matrix in EigTool if possible.  It looks for a directory called 'eigtoollib' in the path or current working directory.  Download it [here](http://www.cs.ox.ac.uk/pseudospectra/eigtool/download/).  It's pretty cool stuff!

- **Zoom to DC/Nyquist:** These are plots focused on the DC and Nyquist eigenvalues.  It should just plot DC and Nyquist eigenvalues on the scales in which they are typically perturbed.  I used these for the zoomed portions of the interpolated operator figures in the paper.  The figure will be zoomed to this region whether or not there are eigenvalues to show there.
 
- **Spectral Radii:** This prints to the console the spectral radius as well as the magnitudes of eigenvalues with largest and smallest real parts--in other words the radii of DC and Nyquist eigenvalues.  There is also a plot of eigenvalue magnitudes arranged on a linear frequency axis.


Disclaimer
----------

This GUI is hacked together from several different scripts used to produce figures for a paper or two.  The code is quite ugly and littered with kludges, but to the best of my knowledge it does what it's supposed to.  However, it is definitely breakable.  In the likely event that something doesn't work correctly, I appreciate it if you contact me.

