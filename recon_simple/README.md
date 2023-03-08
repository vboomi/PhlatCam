# About
A simple reconstruction code for lensless reconstruction. Performs Tikhonov regularized least squares or simply Weiner deconvolution.

# PSF
The psf file is present in 'recon_simple/psfs/' directory. The shape of psf is 912x1368x4x49.  The four channels in order are R,G,G,B. This is the bayer splitting of the captured psf. 
Can disregard one of the Gs and use it as color PSF.

The mask to sensor distance of PhlatCam is about 2 mm.

# Running the reconstruction
An example capture is placed in 'recon_simple/cap/' directory.
Run MATLAB script 'recon_simple/reconTikh.m' for weiner deconvolution reconstruction.
Currently reconstructs the capture file 'recon_simple/cap/hold_tigerr.mat'.

# Sample output
Output from above code as-is is present in 'recon_simple/cap/reconWeiner/'.
