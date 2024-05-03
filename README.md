# PhlatCam: Designed Phase-Mask Based Thin Lensless Camera
Codebase for PhlatCam.

[V. Boominathan](https://vivekboominathan.com/), J. K. Adams, J. T. Robinson and A. Veeraraghavan, "PhlatCam: Designed Phase-Mask Based Thin Lensless Camera," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 42, no. 7, pp. 1618-1629, 1 July 2020. \[[paper](https://doi.org/10.1109/TPAMI.2020.2987489)\]

Free access: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7439257/

# Contents and Instructions

- ./design/ : Design of PhlatCam phase mask.
- ./recon_simple/ : A simple single-step reconstruction code.
- ./refocus/ : Reconstructing refocus stack from single lensless capture.
- ./3D/ : Reconstructing 3D image from signle lensless capture.

Instructions and details can be found inside each of the sub-directory.

# Data
- **3D psf** stack can be downloaded from here: https://www.dropbox.com/sh/rx4hwuygm3ndv7a/AAB6c6DG28txWrKNKfLAH9fTa?dl=0. Contents of the mat file are:
  - psfs - 912x1368x4x49. The four channels in order are R,G,G,B. This is the bayer splitting of the captured psf. Can disregard one of the Gs and use it as color PSF.
  - drng - depth from mask in centimeters. Starts from 1 cm to 20cm.
  - refDepth - Reference depth to which relative magnifications and transforms are calculated
  - mags - relative magnification of PSF wrt PSF at refDepth
  - transForm - relative affine (I think) between PSFs and PSF at refDepth.
  
- Phase mask design used in the paper can be downloaded here: https://www.dropbox.com/s/4mjzidwtfjhdy9a/PhlatCam_used_mask.zip?dl=0.
  - 'phHeight_perlin12_1.00umMaxH_2.00um_7.01mm_1869um_q200nm_lam532_as.png' is the height map stored as an image file. 255 corresponds to 1um.
  - 'phMaskDes_perlin12_2.00um_7.01mm_1869um_lam532_as.mat' has all the data regarding the same mask
  
Links to download sample lensless captures can be found in individual directories.

# Citation

If you use this code, please cite our work:
```
@ARTICLE{PhlatCamVB,
  author={Boominathan, Vivek and Adams, Jesse K. and Robinson, Jacob T. and Veeraraghavan, Ashok},
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence}, 
  title={PhlatCam: Designed Phase-Mask Based Thin Lensless Camera}, 
  year={2020},
  volume={42},
  number={7},
  pages={1618-1629},
  doi={10.1109/TPAMI.2020.2987489}}
```

# Additional publications with this design

### Compact microscope:

Wu, J., Chen, Y., Veeraraghavan, A. et al. "Mesoscopic calcium imaging in a head-unrestrained male non-human primate using a lensless microscope." Nature Communications (2024). \[[paper](https://doi.org/10.1038/s41467-024-45417-6)\]

J.K. Adams, D. Yan, J. Wu, V. Boominathan, et al. "In vivo lensless microscopy via a phase mask generating diffraction patterns with high-contrast contours." Nature Biomedical Engineering 6.5 (2022). \[[paper](https://doi.org/10.1038/s41551-022-00851-z)\]

### Deep learning reconstruction:

Jimin Wu, Vivek Boominathan, Ashok Veeraraghavan, and Jacob T. Robinson, "Real-time, deep-learning aided lensless microscope," Biomed. Opt. Express (2023). \[[paper](https://opg.optica.org/boe/fulltext.cfm?uri=boe-14-8-4037)\]

S. S. Khan, V. Sundar, V. Boominathan, A. Veeraraghavan and K. Mitra, "FlatNet: Towards Photorealistic Scene Reconstruction From Lensless Measurements," in IEEE Transactions on Pattern Analysis and Machine Intelligence (2022) \[[paper](https://arxiv.org/abs/2010.15440)\] \[[project](https://siddiquesalman.github.io/flatnet/)\]

# Review article

V. Boominathan, et al. "Recent advances in lensless imaging." Optica 9.1 (2022). \[[paper](https://doi.org/10.1364/OPTICA.431361)\]
