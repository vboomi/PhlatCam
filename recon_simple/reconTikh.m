%%
showFigs = 1;

%% Reconstruction code
% Set to working directory
% protoDir = 'G:\My Drive\Research\Lensless\PhaseMask\BlkFlyS_20MP\Perlin12_rect\fixed_v2_2018_05_27\2019_08_12\';
% protoPSFDir = 'F:\BoardLevel_PhlatScope\PSFs_bgSub\4D_hydra_depths\';
% protoPSFDir = 'G:\My Drive\Research\Lensless\PhaseMask\BlkFlyS_20MP\Perlin12_rect\fixed_v2_2018_05_27\2019_07_26\calib\';
psfExt = '.mat'; %  '.mat', '.tiff'
capExt = '.mat'; % '.mat', '.tiff', '.avi' % For avi, select only one file

%%
H = 3648; W = 5472; % BFS
% H = 3036; W = 4024; % Basler-ace

%% Params
gamm_l2 = 1000; %4000; %USAF - 500, Scene - 4000

%% Get capture file

[testFile,testCapDir] = deal('hold_tigerr.mat','cap\');
frameImg = importdata([testCapDir, filesep, testFile]);

Ichs = frameImg(:,:,[1,2,4]);

outDir = [testCapDir, filesep, 'reconWeiner/'];
mkdir(outDir)
[~,testName] = fileparts(testFile);

%% Get PSF and test image file names and directories
[psfFile,psfDir] = deal('psf_d20.mat','psfs\');

pfload = load([psfDir, filesep, psfFile]);
psf = pfload.psf;

psf = psf(:,:,[1,2,4]); % raw extract channels R, Gr, and B (skip Gb channel)

[Ny,Nx,Ch] = size(psf);

%% Filter operators
Fx = @(x) fft2(fftshift(x));
FiltX = @(H,x) real(ifftshift(ifft2(H.*Fx(x))));

%% process bulk PSF
py = 5; px = 5; 

% padding fns
pad2d = @(x) padarray(x,[py,px],0,'both');

% Smooth the boundary of psf
w_psfCut = 10;
kg = fspecial('gaussian',w_psfCut*[1,1],2);
crpSmth = zeros([Ny,Nx]);
crpSmth(w_psfCut+1:end-w_psfCut,w_psfCut+1:end-w_psfCut) = 1;
crpSmth = imfilter(crpSmth,kg,'same');
psf = bsxfun(@times, psf, crpSmth);

% Normalized psf with sum val
sms = sum(sum(psf,1),2);
psf = bsxfun(@rdivide, psf,sms);
psf = bsxfun(@times, psf, sms(:,:,:,1));

% norm and pad
% psf = psf*sqrt(Nz)/norm(psf(:));
psf = pad2d(psf);

%% Process input
b = pad2d(Ichs);
bR = b(:,:,1);
bG = b(:,:,2);
bB = b(:,:,3);

%% Setup
% RED
psfR = psf(:,:,1,pp);
% Filter in fourier domain
HsR = Fx(psfR);   %Compute 3D spectrum
HsR_conj = conj(HsR);
% Forward and adjoint of filter
HRfor = @(x) FiltX(HsR,x);
HRadj = @(x) FiltX(HsR_conj,x);
% H'H in fourier domain
HtHR = abs(HsR.*HsR_conj);

% GREEN
psfG = psf(:,:,2,pp);
% Filter in fourier domain
HsG = Fx(psfG);   %Compute 3D spectrum
HsG_conj = conj(HsG);
% Forward and adjoint of filter
HGfor = @(x) FiltX(HsG,x);
HGadj = @(x) FiltX(HsG_conj,x);
% H'H in fourier domain
HtHG = abs(HsG.*HsG_conj);

% BLUE
psfB = psf(:,:,3,pp);
% Filter in fourier domain
HsB = Fx(psfB);   %Compute 3D spectrum
HsB_conj = conj(HsB);
% Forward and adjoint of filter
HBfor = @(x) FiltX(HsB,x);
HBadj = @(x) FiltX(HsB_conj,x);
% H'H in fourier domain
HtHB = abs(HsB.*HsB_conj);

%% Reconstruct
% Recon R
xFilt_mult = 1./(HtHR + gamm_l2);
numerator = HRadj(bR);
Rt_nxt = FiltX(xFilt_mult,numerator);

% Recon G
xFilt_mult = 1./(HtHG + gamm_l2);
numerator = HGadj(bG);
Gt_nxt = FiltX(xFilt_mult,numerator);

% Recon B
xFilt_mult = 1./(HtHB + gamm_l2);
numerator = HBadj(bB);
Bt_nxt = FiltX(xFilt_mult,numerator);

%% Display
Xt_nxt = cat(3,Rt_nxt,Gt_nxt,Bt_nxt);
maxVal = max(Xt_nxt(:));
Xt_nxt = single(Xt_nxt/maxVal);
figure, imshow(Xt_nxt)

%% Save
imwrite(Xt_nxt,sprintf('%s/%s_%gl2_%02d.png',outDir,testName,gamm_l2,pp));