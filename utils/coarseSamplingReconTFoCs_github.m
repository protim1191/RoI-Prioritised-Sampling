function [xrC,y] = coarseSamplingReconTFoCs_github(M, solver,W, imSize, tfocs, xorg,ndev, meas_supplied,y)
% M : coarse sampling matrix
% y1: coarse measurements
% l1 or l2, 1= l1, 2 = l2
% W: sparsity basis, SPOT operator
% imSize: size of image to be reconstructed, 1x2 array
minArgs=8;
maxArgs=9;
narginchk(minArgs,maxArgs)
global Bcg
if ~size(xorg,3) == 1
    error('Image needs to be 2D')
end
xorg = xorg(:);
if meas_supplied
    if exist('y')
        fprintf('Measurements supplied. \n')
    else
        error('Meas_supplied = 1 but measurements not supplied')
    end
else
    fprintf('Generating measurements. \n')
    y0 = M*xorg(:);
    y = y0 + ndev*randn(size(M,1),1);
    tfocs.EPS = norm(y - y0,Inf);
end
switch solver
    case 1
        fprintf('Solving low resolution with l1 .... \n ')
        tfocs.opts.normA2 = linop_normest(linop_spot(M)).^2;
        tfocs.opts.normW2 = linop_normest(linop_spot(W)).^2;
        x0=[]; tfocs.opts.maxIts = 1000;
%         x0 = M'*y;
        xrC = reshape(solver_sBPDN_W( linop_spot(M), linop_spot(W), y, tfocs.EPS, tfocs.mu, x0,[],tfocs.opts),imSize(1),imSize(2));
    case 2
        fprintf('Solving low resolution for l2 .... \n ')
        Bcg = M;
        xrC = reshape(lsqr(Bcg,y), imSize(1), imSize(2));
end

