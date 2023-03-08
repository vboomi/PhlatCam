function Uout = prop2D(Uin,Lxy,lambda,z,method)
%% Author:
% Vivek Boominathan
% Rice University
% vivekb@rice.edu

%%
% propagation - transfer function approach
% 2D
% uniform sampling
% Inputs
% Uin - source plane field
% Lxy = [Lx,Ly] - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% method - 'as' (angular spectrum) or 'fp' (fresnel propagation)
% Output
% Uout - observation plane field

if nargin < 5
    method = 'as';
end
 
[Ny,Nx] = size(Uin);
Lx = Lxy(1); Ly = Lxy(2);

fs_x = Nx/Lx;
fs_y = Ny/Ly;

fx = fs_x * ( (-Nx+1)/2 : (Nx-1)/2 ) / Nx;
fy = fs_y * ( (-Ny+1)/2 : (Ny-1)/2 ) / Ny; fy = fy';

[Fx,Fy] = meshgrid(fx,fy);
bp = sqrt(Fx.^2+Fy.^2) < (1/lambda); % Medium's Bandpass

% Angular spectrum
if strcmpi(method,'as')
    H = bp .* exp(1j*2*pi*(z/lambda)* bp .* sqrt(1-lambda^2*(Fx.^2+Fy.^2)));
end

% Fresnel
if strcmpi(method,'fp')
    H = bp .* exp(1j*2*pi*(z/lambda)* bp .* (1-0.5*lambda^2*(Fx.^2+Fy.^2)));
end

A0 = fft2(Uin); A0 = fftshift(A0);
Az = A0 .* H;

Az = ifftshift(Az);
Uout = ifft2(Az);