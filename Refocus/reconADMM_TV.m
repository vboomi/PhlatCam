%%
addpath('fns/');
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
gamm_TV = 0.0005; %0.0001 -> works
mu_Z = 1; %1;
mu_W = 1; %1;

maxIters = 50;

resid_tol = 1.5;
tau_inc = 1.1;
tau_dec = 1.1;

params = struct('gamm_TV', gamm_TV, ...
    'mu_Z', mu_Z, 'mu_W', mu_W, ...
    'maxIters', maxIters, ...
    'resid_tol', resid_tol, ...
    'tau_inc', tau_inc, ...
    'tau_dec', tau_dec);

%% Get capture file

[testFile,testCapDir] = deal('avgCap20.mat','cap\');
frameImg = importdata([testCapDir, filesep, testFile]);

Ichs = frameImg(:,:,[1,2,4]);

outDir = [testCapDir, filesep, 'reconADMM_TV/'];
mkdir(outDir)
[~,testName] = fileparts(testFile);

%% Get PSF and test image file names and directories
[psfFile,psfDir] = deal('psfs_refocus.mat','psfs\');

% [psfFile,psfDir] = uigetfile([protoPSFDir, '\*.mat'], 'Select Test');
pfload = load([psfDir, filesep, psfFile]);

psfs = pfload.psfs;
drng = pfload.drng;

psfs = psfs(:,:,[1,2,4],:);

% idx = drng >= 3.5 & drng <= 10; %drng==5; % drng > 3.5 & drng < 8; %==20;
% idx = drng == 4.5 | drng == 7 | drng == 11;
% drng = drng(idx);
% psfs = psfs(:,:,:,idx);

[Ny,Nx,Ch,Nz] = size(psfs);

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
psfs = bsxfun(@times, psfs, crpSmth);

% Normalized psf with sum val
sms = mean(sum(sum(psfs,1),2),3);
psfs = bsxfun(@rdivide, psfs,sms);
psfs = bsxfun(@times, psfs, sms(1));

% norm and pad
nrms = sqrt(mean(sum(sum(psfs.^2,1),2),3));
psfs = sqrt(Nz)*bsxfun(@rdivide, psfs, nrms);
psfs = pad2d(psfs);

%% Gradient function-- 2D gradient
% x-grad filter
kx = zeros(Ny+2*py,Nx+2*px); % [Ny+2*py,Nx+2*px] is the 2D size of padded PSF and measurement
kx(1,1) = 1;
kx(1,2) = -1;
Kx = fft2(kx);

% y-grad filter
ky = zeros(Ny+2*py,Nx+2*px);
ky(1,1) = 1;
ky(2,1) = -1;
Ky = fft2(ky);

% Psi - Forward gradient (finite difference)
Psi = @(x) deal(FiltX(Kx,x), FiltX(Ky,x));

% Adjoint of gradient
PsiT = @(P1,P2) FiltX(conj(Kx),P1) + FiltX(conj(Ky),P2);

% Psi'Psi in fourier domain
lapl = zeros(Ny+2*py,Nx+2*px);    %Compute laplacian in closed form. This is the kernal to compute Psi'Psi
lapl(1,1,:) = 4;
lapl(1,2,:) = -1;
lapl(2,1,:) = -1;
lapl(1,end,:) = -1;
lapl(end,1,:) = -1;
PsiTPsi = real(fft2(lapl));   %Compute power spectrum of laplacian

%% Process input
b = pad2d(Ichs);
bR = b(:,:,1);
bG = b(:,:,2);
bB = b(:,:,3);

for pp = 1:Nz
    % RED
    psfR = psfs(:,:,1,pp);
    % Filter in fourier domain
    HsR = Fx(psfR);   %Compute 3D spectrum
    HsR_conj = conj(HsR);
    % Forward and adjoint of filter
    HRfor = @(x) FiltX(HsR,x);
    HRadj = @(x) FiltX(HsR_conj,x);
    % H'H in fourier domain
    HtHR = abs(HsR.*HsR_conj);

    % GREEN
    psfG = psfs(:,:,2,pp);
    % Filter in fourier domain
    HsG = Fx(psfG);   %Compute 3D spectrum
    HsG_conj = conj(HsG);
    % Forward and adjoint of filter
    HGfor = @(x) FiltX(HsG,x);
    HGadj = @(x) FiltX(HsG_conj,x);
    % H'H in fourier domain
    HtHG = abs(HsG.*HsG_conj);

    % BLUE
    psfB = psfs(:,:,3,pp);
    % Filter in fourier domain
    HsB = Fx(psfB);   %Compute 3D spectrum
    HsB_conj = conj(HsB);
    % Forward and adjoint of filter
    HBfor = @(x) FiltX(HsB,x);
    HBadj = @(x) FiltX(HsB_conj,x);
    % H'H in fourier domain
    HtHB = abs(HsB.*HsB_conj);

    %% Initialization
    Xt_prv = zeros(Ny+2*py,Nx+2*px,3);
    [Rt_prv, Gt_prv, Bt_prv] = deal(Xt_prv(:,:,1),Xt_prv(:,:,2),Xt_prv(:,:,3));

    [rhoZRdy_prv,rhoZRdx_prv] = deal(0);
    [rhoZGdy_prv,rhoZGdx_prv] = deal(0);
    [rhoZBdy_prv,rhoZBdx_prv] = deal(0);

    [rhoWR_prv,rhoWG_prv,rhoWB_prv] = deal(0);

    fdataFid = [];
    fregPen = [];
    fobj = [];
    fobj_alt = [];

    primal_res_Z = [];
    dual_res_Z = []; 
    primal_res_W = [];
    dual_res_W = [];

    [dyR_prv, dxR_prv] = Psi(Rt_prv);
    [dyG_prv, dxG_prv] = Psi(Gt_prv);
    [dyB_prv, dxB_prv] = Psi(Rt_prv);

    if showFigs
        fh1 = figure; title 'recon';
%         fh2 = figure; title 'primal residues';
%         fh3 = figure; title 'dual residues';
        fh4 = figure; title 'obj fn';
        drawnow
    end

    %% Iterate
    iter = 0;
    while iter < maxIters
        tic;

        iter = iter + 1;

        %% Primal update
        RFilt_mult = 1./(HtHR + mu_W*PsiTPsi + mu_W);
        GFilt_mult = 1./(HtHG + mu_W*PsiTPsi + mu_W);
        BFilt_mult = 1./(HtHB + mu_W*PsiTPsi + mu_W);

        % RED
%         ZRdy_nxt = shrinkageOp(dyR_prv + rhoZRdy_prv/mu_Z, gamm_TV/mu_Z);
%         ZRdx_nxt = shrinkageOp(dxR_prv + rhoZRdx_prv/mu_Z, gamm_TV/mu_Z);
        [ZRdx_nxt, ZRdy_nxt] = shrinkIso2DTV(dxR_prv + rhoZRdx_prv/mu_Z, dyR_prv + rhoZRdy_prv/mu_Z, gamm_TV/mu_Z);

        WR_nxt = max(Rt_prv + rhoWR_prv/mu_W,0);

        numr = HRadj(bR) + PsiT(mu_Z*ZRdy_nxt - rhoZRdy_prv, mu_Z*ZRdx_nxt - rhoZRdx_prv)...
            + mu_W*WR_nxt - rhoWR_prv;
        Rt_nxt = FiltX(RFilt_mult, numr);

        % GREEN
%         ZGdy_nxt = shrinkageOp(dyG_prv + rhoZGdy_prv/mu_Z, gamm_TV/mu_Z);
%         ZGdx_nxt = shrinkageOp(dxG_prv + rhoZGdx_prv/mu_Z, gamm_TV/mu_Z);
        [ZGdx_nxt, ZGdy_nxt] = shrinkIso2DTV(dxG_prv + rhoZGdx_prv/mu_Z, dyG_prv + rhoZGdy_prv/mu_Z, gamm_TV/mu_Z);

        WG_nxt = max(Gt_prv + rhoWG_prv/mu_W,0);

        numr = HGadj(bG) + PsiT(mu_Z*ZGdy_nxt - rhoZGdy_prv, mu_Z*ZGdx_nxt - rhoZGdx_prv)...
            + mu_W*WG_nxt - rhoWG_prv;
        Gt_nxt = FiltX(GFilt_mult, numr);

        % BLUE
%         ZBdy_nxt = shrinkageOp(dyB_prv + rhoZBdy_prv/mu_Z, gamm_TV/mu_Z);
%         ZBdx_nxt = shrinkageOp(dxB_prv + rhoZBdx_prv/mu_Z, gamm_TV/mu_Z);
        [ZBdx_nxt, ZBdy_nxt] = shrinkIso2DTV(dxB_prv + rhoZBdx_prv/mu_Z, dyB_prv + rhoZBdy_prv/mu_Z, gamm_TV/mu_Z);

        WB_nxt = max(Bt_prv + rhoWB_prv/mu_W,0);

        numr = HBadj(bB) + PsiT(mu_Z*ZBdy_nxt - rhoZBdy_prv, mu_Z*ZBdx_nxt - rhoZBdx_prv)...
            + mu_W*WB_nxt - rhoWB_prv;
        Bt_nxt = FiltX(BFilt_mult, numr);

        %% Next derivatives
        HR_nxt = HRfor(Rt_nxt);
        [dyR_nxt, dxR_nxt] = Psi(Rt_nxt);

        HG_nxt = HGfor(Gt_nxt);
        [dyG_nxt, dxG_nxt] = Psi(Gt_nxt);

        HB_nxt = HBfor(Bt_nxt);
        [dyB_nxt, dxB_nxt] = Psi(Bt_nxt);

        %% Dual update

        % RED
        rpZRdy = dyR_nxt - ZRdy_nxt;
        rpZRdx = dxR_nxt - ZRdx_nxt;
        rhoZRdy_nxt = rhoZRdy_prv + mu_Z*rpZRdy;
        rhoZRdx_nxt = rhoZRdx_prv + mu_Z*rpZRdx;

        rpWR = Rt_nxt - WR_nxt;
        rhoWR_nxt = rhoWR_prv + mu_W*rpWR;

        % GREEN
        rpZGdy = dyG_nxt - ZGdy_nxt;
        rpZGdx = dxG_nxt - ZGdx_nxt;
        rhoZGdy_nxt = rhoZGdy_prv + mu_Z*rpZGdy;
        rhoZGdx_nxt = rhoZGdx_prv + mu_Z*rpZGdx;

        rpWG = Gt_nxt - WG_nxt;
        rhoWG_nxt = rhoWG_prv + mu_W*rpWG;

        % BLUE
        rpZBdy = dyB_nxt - ZBdy_nxt;
        rpZBdx = dxB_nxt - ZBdx_nxt;
        rhoZBdy_nxt = rhoZBdy_prv + mu_Z*rpZBdy;
        rhoZBdx_nxt = rhoZBdx_prv + mu_Z*rpZBdx;

        rpWB = Bt_nxt - WB_nxt;
        rhoWB_nxt = rhoWB_prv + mu_W*rpWB;

        %% Objective fn value
        fdataFid(iter) = 0.5*(norm(bR-HR_nxt,'fro')^2 + norm(bG-HG_nxt,'fro')^2 + norm(bB-HB_nxt,'fro')^2);
        fregPen(iter) = gamm_TV * ( sum(abs(dyR_nxt(:))) + sum(abs(dxR_nxt(:))) ...
            + sum(abs(dyG_nxt(:))) + sum(abs(dxG_nxt(:))) + sum(abs(dyB_nxt(:))) + sum(abs(dxB_nxt(:))) );
        fobj(iter) = fdataFid(iter) + fregPen(iter);
        fobj_alt(iter) = fdataFid(iter) + gamm_TV * ( sum(abs(ZRdy_nxt(:))) + sum(abs(ZRdx_nxt(:))) ...
            + sum(abs(ZGdy_nxt(:))) + sum(abs(ZGdx_nxt(:))) + sum(abs(ZBdy_nxt(:))) + sum(abs(ZBdx_nxt(:))) );;

        %% Residuals
%         primal_res_Z(iter) = sqrt(norm(rpZRdy(:))^2 + norm(rpZRdx(:))^2 ...
%             + norm(rpZGdy(:))^2 + norm(rpZGdx(:))^2 + norm(rpZBdy(:))^2 + norm(rpZBdx(:))^2 );
%         dual_res_Z(iter) = mu_Z*sqrt(norm(dyR_nxt(:)-dyR_prv(:))^2 + norm(dxR_nxt(:)-dxR_prv(:))^2 ...
%             + norm(dyG_nxt(:)-dyG_prv(:))^2 + norm(dxG_nxt(:)-dxG_prv(:))^2 ...
%             + norm(dyB_nxt(:)-dyB_prv(:))^2 + norm(dxB_nxt(:)-dxB_prv(:))^2);
% 
%         primal_res_W(iter) = sqrt(norm(rpWR(:))^2 + norm(rpWG(:))^2 + norm(rpWB(:))^2);
%         dual_res_W(iter) = mu_W*sqrt(norm(Rt_nxt(:)-Rt_prv(:))^2 + norm(Gt_nxt(:)-Gt_prv(:))^2 + ...
%             norm(Bt_nxt(:)-Bt_prv(:))^2);
% 
%         %% Update mu (augmented penalties)
%         [mu_Z, muZ_update] = ...
%             penaltyUpdater(mu_Z,primal_res_Z(iter),dual_res_Z(iter),resid_tol,tau_inc,tau_dec);
%         [mu_W, muW_update] = ...
%             penaltyUpdater(mu_W,primal_res_W(iter),dual_res_W(iter),resid_tol,tau_inc,tau_dec);
% 
%         if muZ_update, disp('mu_Z was updated'); end
%         if muW_update, disp('mu_W was updated'); end
%         fprintf('mu_Z: %0.2f, mu_W: %0.2f.\n',mu_Z,mu_W);

        %% Update previous estimate as the current
        [Rt_prv, Gt_prv, Bt_prv] = deal(Rt_nxt, Gt_nxt, Bt_nxt);

        [rhoZRdy_prv,rhoZRdx_prv] = deal(rhoZRdy_nxt,rhoZRdx_nxt);
        [rhoZGdy_prv,rhoZGdx_prv] = deal(rhoZGdy_nxt,rhoZGdx_nxt);
        [rhoZBdy_prv,rhoZBdx_prv] = deal(rhoZBdy_nxt,rhoZBdx_nxt);

        [rhoWR_prv,rhoWG_prv,rhoWB_prv] = deal(rhoWR_nxt,rhoWG_nxt,rhoWB_nxt);

        [dyR_prv, dxR_prv] = deal(dyR_nxt, dxR_nxt);
        [dyG_prv, dxG_prv] = deal(dyG_nxt, dxG_nxt);
        [dyB_prv, dxB_prv] = deal(dyB_nxt, dxB_nxt);

        toc,

        %%
        Xt_nxt = cat(3, Rt_prv, Gt_prv, Bt_prv);
        Xt_nxt = Xt_nxt/max(Xt_nxt(:));

        %% Display
        if ~mod(iter,10) && showFigs
            figure(fh1), imshow(Xt_nxt); title 'recon';
%             figure(fh2), plot(1:iter,[primal_res_Z',primal_res_W']);  title 'primal residues';
%             figure(fh3), semilogy(1:iter,[dual_res_Z',dual_res_W']); title 'dual residues';
            figure(fh4), semilogy(1:iter,[fdataFid',fregPen',fobj',fobj_alt']); title 'obj fn';
            legend('DataFid', 'Penalty', 'Orig obj', 'Alt obj');
            drawnow;
        end
    end
    
    %% Concat
    Xt_nxt = cat(3, Rt_prv, Gt_prv, Bt_prv);
    maxVal = max(Xt_nxt(:));
    Xt_nxt = single(Xt_nxt/maxVal);
    figure, imshow(Xt_nxt)

    %%
    save(sprintf('%s/%s_iters%03d_TV%g_%02d_d%g.mat',outDir,testName,maxIters,gamm_TV,pp,drng(pp)), ...
        'Xt_nxt','gamm_TV','maxIters','params','maxVal');
    imwrite(Xt_nxt,sprintf('%s/%s_iters%03d_TV%g_%02d_d%g.png',outDir,testName,maxIters,gamm_TV,pp,drng(pp)));
end