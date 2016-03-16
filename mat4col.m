function nm = mat4col(m,km)
%
% function ncol = mat4col(m,km)
%
% MATLAB function to convert a matrix to a large column (or a new matrix)
% based on the mask and the focal point within the mask
%
% Based off: mat3col
%
% Input
%   - m  -- original matrix, at least 4D
%   - km -- 4D kernel mask
%
% Output
%   - nm -- New matrix. The number of columns equals the number of non-zero
%           elements in km; the number of rows equals the number
%           of non-equivalent kernel-sized blocks m can fit.
%
%   Yulin V. Chang, 20150610
%   Modified to use nnz instead of sum(sum(sum(...)))              20150915

[N1,N2,N3,N4] = size(m);
[k1,k2,k3,k4] = size(km);

nnzel = nnz(km);

nm = zeros((N1-k1+1)*(N2-k2+1)*(N3-k3+1)*(N4-k4+1),nnzel);

for ll = 1:N4-k4+1
    for kk = 1:N3-k3+1
        for jj = 1:N2-k2+1
            for ii = 1:N1-k1+1
                movb = m(ii:ii+k1-1,jj:jj+k2-1,kk:kk+k3-1,ll:ll+k4-1);
                nm((ll-1)*(N3-k3+1)*(N2-k2+1)*(N1-k1+1)+(kk-1)*(N2-k2+1)*(N1-k1+1)+(jj-1)*(N1-k1+1)+ii,:) = movb(km==1).';
            end
        end
    end
end
