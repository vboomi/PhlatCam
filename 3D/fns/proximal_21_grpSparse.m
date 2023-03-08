function x = proximal_21_grpSparse(z, gamma, wt)

nrm = sqrt(sum(sum(z.^2,1),2));
shrink_nrm = shrinkageOp(nrm, wt*gamma);

zn = bsxfun(@rdivide, z, nrm);
zn(isnan(zn)) = 0;

x = bsxfun(@times, zn, shrink_nrm);