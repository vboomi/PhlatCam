function x = proximal_12(z, gamma, wt)
tic;
%   prox_L12(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * || z ||_12^2
%
%   where 
%
%       ' || x ||_12 =  sqrt ( sum_j ( sum_i w(i,j)*|x(i,j)|)^2  )

% l1 on 3rd (depth) dimension
% l2 on (x,y) dim together

% z - i x j x k
% wght - 1 x 1 x k

[H,W,B] = size(z);

z = reshape(z,[H*W, B]);
z = z'; % k x ij
wt = reshape(wt, [B,1]);
wt = repmat(wt,[1, H*W]);

r = bsxfun(@rdivide, abs(z), wt);

[rk, Irk] = sort(r,1,'descend');

wk = zeros(size(wt));
for ij = 1:H*W
    wk(:,ij) = wt(Irk(:,ij),ij);
end

wk2 = wk.^2;
wk2rk = wk2.*rk;

% Find Lk
wk2cs = cumsum(wk2,1);
wk2rkcs = cumsum(wk2rk,1);

rcmpt = gamma*wk2rkcs./(1+gamma*wk2cs);

rcmp = (rcmpt - rk) >= 0;

[~,Lk] = max(rcmp,[],1);
Lk = Lk-1;
Lk(Lk==0) = 1;

% find lllzklll_wk and Lwk
zk_wk = zeros(1,H*W);
Lwk = zeros(1,H*W);
for ij = 1:H*W
    zk_wk(ij) = sum(wk2rk(1:Lk(ij),ij));
    Lwk(ij) = sum(wk2(1:Lk(ij),ij));
end

% Thresholds
xi_k = bsxfun(@times, gamma*wt, zk_wk./(1+gamma*Lwk));
x = shrinkageOp(z, xi_k);

x = x'; % ij x k
x = reshape(x, [H,W,B]);
toc,
end



