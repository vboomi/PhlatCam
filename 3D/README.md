# About
Reconstruct 3D image stack from a single capture.

# 3D PSF
The mat file of 3D PSF can be downloaded here: https://www.dropbox.com/sh/rx4hwuygm3ndv7a/AAB6c6DG28txWrKNKfLAH9fTa?dl=0

Place the mat file in '3D/psfs' directory.

# Contents of 3D PSF mat file
- psfs - 912x1368x4x49.  The four channels in order are R,G,G,B. This is the bayer splitting of the captured psf. Can disregard one of the Gs and use it as color PSF.
- drng - depth from mask in centimeters. Starts from 1 cm to 20cm.
- refDepth - Reference depth to which relative magnifications and transforms are calculated
- mags - relative magnification of PSF wrt PSF at refDepth
- transForm - relative affine (I think) between PSFs and PSF at refDepth.

The mask to sensor distance of PhlatCam is about 2 mm.

# Running the reconstruction
An example capture is placed in '3D/cap/' directory.
Run MATLAB script '3D/reconADMM_3D_2DTV_Sparse.m' for 3D reconstruction.

# Sample output
Output from above code as-is can be downloaded here: https://www.dropbox.com/sh/ulff6868d9yzhrn/AADCOUxDgWnvmTqurxUDFEL4a?dl=0
