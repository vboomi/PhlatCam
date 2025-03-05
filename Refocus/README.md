# About
Reconstruction refocus stack from single capture.

# PSF Stack
The refocus psf stack can be downloaded from here: https://www.dropbox.com/sh/rms5pamqomjbtgy/AAA_xgEQ9Phbv9FcJYRNnCTya?dl=0

Place the mat file in 'refocus/psfs' directory.

# Contents of PSF stack mat file
- psfs - 912x1368x4x3.  The four channels in order are R,G,G,B. This is the bayer splitting of the captured psf. Can disregard one of the Gs and use it as color PSF.
- drng - depth from mask in centimeters. 4.5cm, 7 cm, 11 cm.

The mask to sensor distance of PhlatCam is about 2 mm.

# Running the reconstruction
An example capture is placed in 'refocus/cap/' directory.
Run MATLAB script 'refocus/reconADMM_TV.m' for refocus reconstruction.

# Sample output
Output from above code as-is is present in 'refocus/cap/selectRecon/'
