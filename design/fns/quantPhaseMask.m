function [phMhq, phMq] = quantPhaseMask(phM,lambd,dn,dh)
%% Author:
% Vivek Boominathan
% Rice University
% vivekb@rice.edu

%%
% lambd - wavelength in um (usually 0.532 um)
% dn - change in refractive index (usually 0.5)
% dh - quantized height in um (usually 0.2 um)

%%
c = lambd/(2*pi*dn);
phMh = phM*c;

phMhq = dh*floor(phMh/dh);
phMq = phMhq/c;
figure, imagesc(phMhq), colormap jet, axis image

fprintf('Max height - \nContinuous: %0.2fum \nQuantized: %0.2fum\n',max(phMh(:)),max(phMhq(:)));