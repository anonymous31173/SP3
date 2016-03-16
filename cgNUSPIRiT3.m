function [sol, RESVEC] = cgNUSPIRiT3(u, aid, NFFT, SP3, nit, lmd, x0)
%
% [res, RESVEC] = cgNUSPIRiT3(u, dw, NFFT, SP3, nit, lmd, lev, mth, thr, x0);
%
% A conjugate gradient implementation of 3D image-domain Non-uniform SPIRiT 
%
% Inputs
%   u   - undersampled data
%   dw  - density weighting
%   NFFT- Non-uniform FFT operator
%   SP3 - 3D SPIRiT kernel
%   nit - number of iterations
%   lmd - splitting parameter for the above expression. Based on what is
%       used in cgNUSPIRiT, could use 1 to start with
%   x0  - inital image (3D)
%
% Outputs
%   sol - image solution
%   RESVEC - Residuals
%
% Yulin V. Chang, 20150714
%


dsz = size(u);
N   = numel(x0);
msz = size(x0);

naq = nnz(aid);

b = [u(aid); zeros(N,1)];
switch getMth(SP3)
    case 'image'
        [sol, ~, ~, ~, RESVEC, ~] = lsqr(@(x,tflag)afun(x,aid,NFFT,SP3,dsz,msz,naq,lmd,tflag), b, [], nit, speye(N,N), speye(N,N), x0(:));
        sol = reshape(sol, msz);
    case 'conv'
        [sol, ~, ~, ~, RESVEC, ~] = lsqr(@(x,tflag)bfun(x,aid,NFFT,SP3,dsz,msz,naq,lmd,tflag), b, [], nit, speye(N,N), speye(N,N), x0(:));
        sol = reshape(sol,msz);
        sol = ifft3vc(sol);
end
function [y, tflag] = afun(x, aid, NFFT, GOP, dsz, msz, naq, lmd, tflag)
if strcmp(tflag,'transp') 
    x1 = zeros(dsz);
    x1(aid) = x(1:naq);
    x2 = reshape(x(naq+1:end), msz);
    y  = NFFT'*ifftc(x1,3) + lmd*(GOP'*x2);
    y  = y(:);
else
    x  = reshape(x,msz);
    y1 = fftc(NFFT*x,3);
    y2 = GOP*x;
    y  = [y1(aid); y2(:)*lmd];
end

function [y, tflag] = bfun(x, aid, NFFT, GOP, dsz, msz, naq, lmd, tflag)
if strcmp(tflag,'transp')
    x1 = zeros(dsz);
    x1(aid) = x(1:naq);
    x2 = reshape(x(naq+1:end), msz);
    y  = fft2vc(NFFT'*x1) + lmd*(GOP'*x2);
    y  = y(:);
else
    x  = reshape(x,msz);
    y1 = NFFT*(ifft2vc(x));
    y2 = GOP*x;
    y  = [y1(aid); y2(:)*lmd];
end

