%%
addpath('../fns/');
showFigs = 1;

%% Reconstruction code
psfExt = '.mat'; %  '.mat', '.tiff'
capExt = '.mat'; % '.mat', '.tiff', '.avi' % For avi, select only one file

%%
H = 3648; W = 5472; % BFS
% H = 3036; W = 4024; % Basler-ace

%% Params
gamm_TV2D = 0.001;
gamm_sp = 0.005;%  0.5
mu_V = 1;
mu_Z = 1;
mu_W = 1;
mu_T = 1; %1

strSparseSet = {'grpSparse','depthSparse','fullSparse'};
strSparse = 'fullSparse'; 
switch strSparse
    case strSparseSet{1}
        disp('Executing group sparsity on depth');
    case strSparseSet{2}
        disp('Executing depth sparsity');
    case strSparseSet{3}
        disp('Executing full sparsity');
end

color = 2; % 1- RED, 2- GREEN, 3- BLUE

maxIters = 50;

resid_tol = 1.5;
tau_inc = 1.5; %1.1;
tau_dec = 1.5; %1.1;

params = struct('gamm_TV2D',gamm_TV2D,'gamm_sp', gamm_sp,...
    'mu_V', mu_V, 'mu_Z', mu_Z, ...
    'mu_W', mu_W, 'mu_T', mu_T, ...
    'maxIters', maxIters, ...
    'resid_tol', resid_tol, ...
    'tau_inc', tau_inc, ...
    'tau_dec', tau_dec);

%% Get capture file
[testFile,testCapDir] = deal('avgCap4.mat','cap\');

% [testFile,testCapDir] = uigetfile([protoDir, '\*.mat'], 'Select Test');
frameImg = importdata([testCapDir, filesep, testFile]);

frameImg = frameImg(:,:,[1,2,4]);
Ichs = frameImg(:,:,color);

%% Get PSF and test image file names and directories
[psfFile,psfDir] = deal('psfs_save_magfs.mat','psfs\');

% [psfFile,psfDir] = uigetfile([protoPSFDir, '\*.mat'], 'Select Test');
pfload = load([psfDir, filesep, psfFile]);

psfs = pfload.psfs;
drng = pfload.drng;
magfs = pfload.mags;

psfs = psfs(:,:,[1,2,4],:);
%%
% idx = drng >= 4 & drng <= 20; %drng >= 2.5 & drng <= 14; %drng >= 3.5 & drng <= 8; %==20;
% idx = drng >= 4 & drng <= 10;
% idx = drng >= 2.5 & drng <= 14; % <-- ruler
% idx = drng >=4.5 & drng <=20; % <-- blocks + man
idx = drng <= 3; % <--- LIGHT
% idx = drng <= 5; % 
drng = drng(idx);
psfs = psfs(:,:,:,idx);
magfs = magfs(idx);

% [Ny,Nx,Ch,Nz] = size(psfs);

psfs = squeeze(psfs(:,:,color,:));
[Ny,Nx,Nz] = size(psfs);

invZ2 = 1./drng.^2;
invZ2 = Nz*invZ2/sum(invZ2);
wght = zeros(1,1,Nz);
wght(:) = invZ2;
% no depth preference
% wght(:) = 1;

% magfs(:) = 1; %<-- debug, comment out
dpth_scl = 1./magfs;
% dpth_scl = dpth_scl/sum(dpth_scl); %<-- debug, uncomment
dpth_scl = dpth_scl/max(dpth_scl);

dpth_scl_temp = zeros(1,1,Nz);
dpth_scl_temp(:) = dpth_scl(:);
dpth_scl = dpth_scl_temp;

%% Out dir
outDir = [testCapDir,'/reconADMM_magf_wt_2DTV_',strSparse,'-clr',num2str(color),filesep];
mkdir(outDir)
[~,testName] = fileparts(testFile);

%% Filter operators
fftshift2 = @(x) fftshift(fftshift(x,1),2);
ifftshift2 = @(x) ifftshift(ifftshift(x,1),2);

Fx2 = @(x) fft2(fftshift2(x));
FiltX2 = @(H,x) real(ifftshift2(ifft2(H.*Fx2(x))));

Fxn = @(x) fftn(fftshift(x));
FiltXn = @(H,x) real(ifftshift(ifftn(H.*Fxn(x))));

%% Process PSFs
py = 5; px = 5; 

% padding fns
pad2d = @(x) padarray(x,[py,px],0,'both');

% Smooth the boundary of psf
w_psfCut = 100; %10
kg = fspecial('gaussian',w_psfCut*[1,1],w_psfCut/10); %2);
crpSmth = zeros(Ny,Nx);
crpSmth(w_psfCut+1:end-w_psfCut,w_psfCut+1:end-w_psfCut) = 1;
crpSmth = imfilter(crpSmth,kg,'same');
psfs = bsxfun(@times, psfs, crpSmth);

% Normalized psf with sum val - color
% sms = mean(sum(sum(psfs,1),2),3);
% psfs = bsxfun(@rdivide, psfs,sms);
% psfs = bsxfun(@times, psfs, sms(:,:,:,1));

% % Normalized psf with sum val - mono
sms = sum(sum(psfs,1),2);
psfs = bsxfun(@rdivide, psfs,sms);
psfs = 500*psfs;
% psfs = bsxfun(@times, psfs, sms(:,:,1));
psfs = pad2d(psfs);

% % norm and pad - color
% nrms = sqrt(mean(sum(sum(psfs.^2,1),2),3));
% psfs = sqrt(Nz)*bsxfun(@rdivide, psfs, nrms);
% psfs = pad2d(psfs);

% % norm and pad - mono
% nrms = norm(psfs(:));
% psfs = sqrt(Nz)*bsxfun(@rdivide, psfs, nrms);
% psfs = pad2d(psfs);

% % Normalized psf with norm val and pad - mono
% nrms = sqrt(sum(sum(psfs.^2,1),2));
% psfs = bsxfun(@rdivide, psfs,nrms);
% psfs = pad2d(psfs);

% Filter in fourier domain
% Hd = Fx2(psfs);
Hd = bsxfun(@times, Fx2(psfs), dpth_scl);
Hd_conj = conj(Hd);
% Forward and adjoint of filter
Hdfor = @(x) FiltX2(Hd,x);
Hdadj = @(x) FiltX2(Hd_conj,x);
% H'H in fourier domain
HdtHd = abs(Hd.*Hd_conj);

%% Summation operator
Sfor = @(x) sum(x,3); % Summation forward
Sadj = @(x) repmat(x,[1,1,Nz]); % Summation adjoint

ss = zeros(size(psfs));
ss(1,1,:) = 1;
StS = real(fftn(ss));

%% Gradient function-- 2D gradient
% x-grad filter
kx = zeros(size(psfs));
kx(1,1,:) = 1;
kx(1,2,:) = -1;
Kx = fft2(kx);

% y-grad filter
ky = zeros(size(psfs));
ky(1,1,:) = 1;
ky(2,1,:) = -1;
Ky = fft2(ky);

% Psi - Forward gradient (finite difference)
Psi = @(x) deal(FiltX2(Kx,x), FiltX2(Ky,x));

% Adjoint of gradient
PsiT = @(P1,P2) FiltX2(conj(Kx),P1) + FiltX2(conj(Ky),P2);

% Psi'Psi in fourier domain
lapl = zeros(size(psfs));    %Compute laplacian in closed form. This is the kernal to compute Psi'Psi
lapl(1,1,:) = 4;
lapl(1,2,:) = -1;
lapl(2,1,:) = -1;
lapl(1,end,:) = -1;
lapl(end,1,:) = -1;
PsiTPsi = real(fft2(lapl));   %Compute power spectrum of laplacian

%% Process input
b = pad2d(Ichs);
% bR = b(:,:,1);
% bG = b(:,:,2);
% bB = b(:,:,3);
Stb = Sadj(b);

%% Initialization
Xt_prv = zeros(size(psfs));

rhoV_prv = 0;
[rhoZdx_prv,rhoZdy_prv] = deal(0);
rhoW_prv = 0;
rhoT_prv = 0;

fdataFid = [];
fregPen = [];
fobj = [];
fobj_alt = [];

primal_res_V = [];
dual_res_V = [];
primal_res_Z = [];
dual_res_Z = []; 
primal_res_W = [];
dual_res_W = [];
primal_res_T = [];
dual_res_T = [];

HX_prv = Hdfor(Xt_prv);
[dx_prv, dy_prv] = Psi(Xt_prv);

if showFigs
    fh1 = figure; title 'recon';
    fh2 = figure; title 'primal & dual residues';
    fh4 = figure; title 'obj fn';
    drawnow;
end

%% Iterate
iter = 0;
while iter < maxIters
    tic;

    iter = iter + 1;

    %% Primal update
    vFilt_mult = 1./(StS + mu_V);
    numr = Stb + mu_V*HX_prv + rhoV_prv;
    V_nxt = FiltXn(vFilt_mult,numr);
    
    Zdx_nxt = shrinkageOp(dx_prv + rhoZdx_prv/mu_Z, gamm_TV2D/mu_Z);
    Zdy_nxt = shrinkageOp(dy_prv + rhoZdy_prv/mu_Z, gamm_TV2D/mu_Z);

    W_nxt = max(Xt_prv + rhoW_prv/mu_W,0);

    switch strSparse
        case strSparseSet{1}
            T_nxt = proximal_21_grpSparse(Xt_prv + rhoT_prv/mu_T, gamm_sp/mu_T, wght);
        case strSparseSet{2}
            T_nxt = proximal_12(Xt_prv + rhoT_prv/mu_T, gamm_sp/mu_T, wght);
        case strSparseSet{3}
            T_nxt = proximal_L1(Xt_prv + rhoT_prv/mu_T, gamm_sp/mu_T, wght);
    end

    xFilt_mult = 1./(mu_V*HdtHd + mu_Z*PsiTPsi + mu_W + mu_T);
    numr = Hdadj(mu_V*V_nxt - rhoV_prv) ...
        + PsiT(mu_Z*Zdx_nxt - rhoZdx_prv, mu_Z*Zdy_nxt - rhoZdy_prv) ...
        + mu_W*W_nxt - rhoW_prv + mu_T*T_nxt - rhoT_prv;
    Xt_nxt = FiltX2(xFilt_mult, numr);

    %% Next derivatives
    HX_nxt = Hdfor(Xt_nxt);
    [dx_nxt, dy_nxt] = Psi(Xt_nxt);

    %% Dual update
    rpV = HX_nxt - V_nxt;
    rhoV_nxt = rhoV_prv + mu_V*rpV;

    rpZdx = dx_nxt - Zdx_nxt;
    rpZdy = dy_nxt - Zdy_nxt;
    rhoZdx_nxt = rhoZdx_prv + mu_Z*rpZdx;
    rhoZdy_nxt = rhoZdy_prv + mu_Z*rpZdy;

    rpW = Xt_nxt - W_nxt;
    rhoW_nxt = rhoW_prv + mu_W*rpW;
    
    rpT = Xt_nxt - T_nxt;
    rhoT_nxt = rhoT_prv + mu_T*rpT;

    %% Objective fn value
    switch strSparse
        case strSparseSet{1}
            penSparse = @(x) sum(squeeze(sqrt(sum(sum(x.^2,1),2)).*wght));
        case strSparseSet{2}
            penSparse = @(x) norm(sum(bsxfun(@times,abs(x),wght),3), 'fro')^2;
        case strSparseSet{3}
            penSparse = @(x) sum(squeeze(sum(sum(abs(x),1),2).*wght));
    end
    
    fdataFid(iter) = 0.5*(norm(b-Sfor(HX_nxt),'fro')^2);
    fregPen(iter) = gamm_sp * penSparse(Xt_nxt) ...
                     + gamm_TV2D * ( sum(abs(dx_nxt(:))) + sum(abs(dy_nxt(:))) );
    fobj(iter) = fdataFid(iter) + fregPen(iter);
    fobj_alt(iter) = 0.5*(norm(b-Sfor(V_nxt),'fro')^2) ...
        + gamm_sp *  penSparse(T_nxt) ...
        + gamm_TV2D * ( sum(abs(Zdx_nxt(:))) + sum(abs(Zdy_nxt(:))) );

    % Residuals
    primal_res_V(iter) = norm(rpV(:));
    dual_res_V(iter) = mu_V*norm(HX_nxt(:)-HX_prv(:));
    
    primal_res_Z(iter) = sqrt(norm(rpZdx(:))^2 + norm(rpZdy(:))^2);
    dual_res_Z(iter) = mu_Z*sqrt( norm(dx_nxt(:)-dx_prv(:))^2 ...
                    + norm(dy_nxt(:)-dy_prv(:))^2 );

    primal_res_W(iter) = norm(rpW(:));
    dual_res_W(iter) = mu_W*norm(Xt_nxt(:)-Xt_prv(:));
    
    primal_res_T(iter) = norm(rpT(:));
    dual_res_T(iter) = mu_T*norm(Xt_nxt(:)-Xt_prv(:));

    % Update mu (augmented penalties)
    [mu_V, muV_update] = ...
        penaltyUpdater(mu_V,primal_res_V(iter),dual_res_V(iter),resid_tol,tau_inc,tau_dec);
    [mu_Z, muZ_update] = ...
        penaltyUpdater(mu_Z,primal_res_Z(iter),dual_res_Z(iter),resid_tol,tau_inc,tau_dec);
    [mu_W, muW_update] = ...
        penaltyUpdater(mu_W,primal_res_W(iter),dual_res_W(iter),resid_tol,tau_inc,tau_dec);
    [mu_T, muT_update] = ...
        penaltyUpdater(mu_T,primal_res_T(iter),dual_res_T(iter),resid_tol,tau_inc,tau_dec);
    
%     if mu_V < 1
%         mu_V = 1;
%     end
    
%     if muV_update, disp('mu_V was updated'); end
%     if muU_update, disp('mu_U was updated'); end
%     if muZ_update, disp('mu_Z was updated'); end
%     if muW_update, disp('mu_W was updated'); end
%     if muT_update, disp('mu_T was updated'); end
    fprintf('mu_V: %0.2f, mu_Z: %0.2f, mu_W: %0.2f, mu_T: %0.2f.\n',...
            mu_V, mu_Z, mu_W, mu_T);

    %% Update previous estimate as the current
    Xt_prv = Xt_nxt;
    HX_prv = HX_nxt;  
    [dx_prv, dy_prv] = deal(dx_nxt, dy_nxt);

    rhoV_prv = rhoV_nxt;
    [rhoZdx_prv,rhoZdy_prv] = deal(rhoZdx_nxt,rhoZdy_nxt);
    rhoW_prv = rhoW_nxt;
    rhoT_prv = rhoT_nxt;

    toc,            
    %% Display
    if ~mod(iter,1)
        Xt_tmp = Xt_nxt;
%         figure(fh1), imshow(max(Xt_tmp/max(Xt_tmp(:)),[],3)); title 'recon';

        figure(fh2), title 'primal & dual residues';
        subplot(4,1,1),semilogy(1:iter,[primal_res_V',dual_res_V'])
        subplot(4,1,2),semilogy(1:iter,[primal_res_Z',dual_res_Z'])
        subplot(4,1,3),semilogy(1:iter,[primal_res_W',dual_res_W'])
        subplot(4,1,4),semilogy(1:iter,[primal_res_T',dual_res_T'])

        figure(fh4), semilogy(1:iter,[fdataFid',fregPen',fobj',fobj_alt']); title 'obj fn';
        legend('DataFid', 'Penalty', 'Orig obj', 'Alt obj');
        drawnow;
    end
    if ~mod(iter,10)
        implay(Xt_nxt/max(Xt_nxt(:)));
        maxVal = max(Xt_nxt(:));
        Xt_Stack = Xt_nxt/maxVal; %Scale through depth
        save(sprintf('%s/%s_iters%03d_%gTV2D_%g%s_%g-%g.mat',outDir,testName,iter,gamm_TV2D,gamm_sp,strSparse,drng(1),drng(end)), 'Xt_Stack','maxVal','drng','gamm_TV2D','gamm_sp','strSparse','iter','params','-v7.3');
    end
            
end

%% Normalize to image range
Xt_nxt = single(Xt_nxt);

maxVal = max(Xt_nxt(:));
Xt_Stack = Xt_nxt/maxVal; %Scale through depth
implay(Xt_Stack);

%% Saves etc.
% for yy = 1:size(psfs,3)
%     imwrite(Xt_Stack(:,:,yy),sprintf('%s/%s_iters%03d_%gTV2D_%gSparse_psf%gf.png',outDir,testName,maxIters,gamm_tikh,gamm_sp,drng(yy)));   
% end
save(sprintf('%s/%s_iters%03d_%gTV2D_%g%s_%g-%g.mat',outDir,testName,iter,gamm_TV2D,gamm_sp,strSparse,drng(1),drng(end)), 'Xt_Stack','maxVal','drng','gamm_TV2D','gamm_sp','strSparse','iter','params','-v7.3');