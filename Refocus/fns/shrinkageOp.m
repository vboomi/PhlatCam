function s = shrinkageOp(z,lmbd)
%% Author:
% Vivek Boominathan
% Rice University
% vivekb@rice.edu

%%
s = max(abs(z)-lmbd,0).*sign(z);