function [dx_nxt, dy_nxt] = shrinkIso2DTV(dx, dy, tau)

% Perform isotropic soft thresholding
mag = sqrt(dx.^2+dy.^2);

magt = shrinkageOp(mag,tau);

mmult = magt./mag;
mmult(mag==0) = 0;

dx_nxt = dx.*mmult;
dy_nxt = dy.*mmult;