function [phM, Mm, MsA] = genPhaseMask(psf, lambd, pxSz, thickness, numIters, method)
%% Author:
% Vivek Boominathan
% Rice University
% vivekb@rice.edu

%% Input
% lambd - wavelength in um
% pxSz - pixel size in um
% thickness - thickness in um
% method - 'as' or 'fp'

%% Output
% phM - Phase map [0,2pi)
% Mm - Field at sensor plane
% MsA - Intensity or PSF at the sensor plane

%%
zMS = thickness;
Mamp = sqrt(psf);
Ms = Mamp;

netLenXY = size(psf)*pxSz;

fh = figure;

% figure(fh),
for ii=1:numIters
    Mm = prop2D(Ms,netLenXY,lambd,-zMS,method);
    Mm = Mm./abs(Mm);
    
    % handle nan
    Mm(isnan(Mm)) = 0;
    
    phM = angle(Mm);
    phM(phM<0) = 2*pi + phM(phM<0);
    if min(phM(:)) > pi
        phM = phM - pi;
    end
    
    figure(fh), subplot(1,2,1), imagesc(phM), axis image, colormap jet; colorbar;
    title(sprintf('Iter %d',ii))
    
    Ms = prop2D(Mm,netLenXY,lambd,zMS,'as');
    MsA = abs(Ms).^2;
    Ms = Mamp.*Ms./sqrt(MsA);
    
    subplot(1,2,2), imshow(MsA,[]);
    drawnow;
end
close(fh);
figure, imshow(MsA,[])

Mm = single(Mm);
MsA = single(MsA);

phM = angle(Mm);
phM(phM<0) = 2*pi + phM(phM<0);
% if min(phM(:)) > pi
%     phM = phM - pi;
% end