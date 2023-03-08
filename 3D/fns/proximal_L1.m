function x = proximal_L1(z, gamma, wt)

x = max(bsxfun(@minus,abs(z),gamma*wt),0).*sign(z);