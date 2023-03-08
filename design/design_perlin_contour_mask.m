%% Author:
% Vivek Boominathan
% Rice University
% vivekb@rice.edu
% "PhlatCam: Designed Phase-Mask Based Thin Lensless Camera", TPAMI 2020

%% Paths
addpath('fns/');

masksOutDir = 'masks_out';
mkdir(masksOutDir)
%% %%%%%% Contour PSF generation
pxSz = 2; %um Pixel size of PSF and phase mask

%% Create Perlin contour pattern

pName = 'perlin12_20_example';
minFeature = 12; %um Width of contour

% Fix random seed
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% Perlin noise generator
pattern_size = [300, 300]; % pixels
In = perlin2D(pattern_size,[20,20]); % subsmpl - [20,20] was prototype params

% Contour PSF creation
In = In - min(In(:));
In = In/max(In(:));
Mp = edge(In,'canny');
M = double(Mp);

% Resize to get desired contour width
M = imresize(M,minFeature/pxSz,'nearest');
figure, imshow(M)

%% %%%%%% Phase mask design

%% Camera thickness (mask-sensor distance)
T = 1869; %um 
% After fabrication, the optimal mask-sensor distance can vary a little. In
% my case ~1900 um for prototype

%% Generate phase mask
zMS = T;
lambd = 0.532; %um
method = 'as'; % as - angular spectrum
numIters = 20; % 20 to 100
netLenXY = size(M)*pxSz; %um % Physical size of mask

[phMm, Mm, MsA] = genPhaseMask(M, lambd, pxSz, zMS, numIters, method); % phMm is the phase map

%% Save continuos phase profile
save(fullfile(masksOutDir,sprintf('phMaskDes_%s_%0.2fum_%0.2fmm_%dum_lam%d_%s.mat',pName,pxSz,netLenXY(1)/1000,T,lambd*1000,method)),'Mm','MsA','phMm','M',...
    'pxSz','netLenXY','zMS','lambd','-v7.3');

%% %%%%%% Quantizing phase mask height for fabrication

%% Quantization
dn = 0.5103; % refractive index difference of material wrt air
dh = 0.2; %um % quantized height step size

netLenXY = size(phMm)*pxSz;
[phMmHq,phMmq] = quantPhaseMask(phMm,lambd,dn,dh);
Ms = prop2D(exp(1i*phMmq),netLenXY,lambd,zMS,'as');
Ms = abs(Ms).^2;
figure, imshow(Ms,[]);

maxH = max(phMmHq(:));

%% Save quantized height map as image
imwrite(phMmHq/max(phMmHq(:)),fullfile(masksOutDir,sprintf('phHeight_%s_%0.2fumMaxH_%0.2fum_%0.2fmm_%dum_q%dnm_lam%d_%s.png',pName,maxH,pxSz,netLenXY(1)/1000,zMS,dh*1000,lambd*1000,method)));