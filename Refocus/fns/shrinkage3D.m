function [dx_nxt, dy_nxt, dz_nxt] = shrinkage3D(dx, dy, dz, tau)

% Perform isotropic soft thresholding
mag = sqrt(dx.^2+dy.^2+dz.^2);

magt = shrinkageOp(mag,tau);

mmult = magt./mag;
mmult(mag==0) = 0;

dx_nxt = dx.*mmult;
dy_nxt = dy.*mmult;
dz_nxt = dz.*mmult;