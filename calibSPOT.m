function kn = calibSPOT(cal, ksh, tik, mth)
%
% kn = calibSPOT(cal, ksh, tik, mth);
%
% Calibrate SPOT kernel. Works for both 2D and 3D
%
% Inputs
%   cal -- calibration data. Must be 3D for a 2D kernel, and 4D for a 3D
%           kernel
%   ksh -- kernel shape/geometry. 2D or 3D
%   tik -- Tikhonov regularization
%   mth -- method for solving kernel equations. Default in pinv
%           - 'cg'   = congjugate gradient (bicgstab)
%           - 'pinv' = pseudo inversion
%
% Output
%   kn  -- SPOT kernel
%
% (c) Yulin V. Chang, University of Pennsylvania, 20151004
%

ksz = size(ksh);
csz = size(cal);
dim = length(ksz);
assert((dim==2 && length(csz)==3)||(dim==3 && length(csz)==4), 'calibration data must be 3D or 4D, ksh must be 2D or 3D');

Q   = csz(1,end);
nzl = nnz(ksh);
ctr = (nzl+1)/2;

switch dim
    case 2
        soub = mat3col(cal, repmat(ksh,[1 1 Q]));
    case 3
        soub = mat4col(cal, repmat(ksh,[1 1 1 Q]));
end
fsm = soub'*soub;

kn = zeros([ksz Q Q]);

for ii = 1:Q
    switch dim
        case 2
            fkmsk = repmat(ksh,[1 1 Q]);
            fkmsk((ksz(1)+1)/2,(ksz(2)+1)/2,ii) = 0;
        case 3
            fkmsk = repmat(ksh,[1 1 1 Q]);
            fkmsk((ksz(1)+1)/2,(ksz(2)+1)/2,(ksz(3)+1)/2,ii) = 0;
    end
    tgt  = (ii-1)*nzl+ctr;
    smsk = ones(nzl*Q,1);
    smsk(tgt)=0;
    smx = fsm(smsk==1,:);
    smx = smx(:,smsk==1);
    tmx = fsm(:,tgt);
    tmx(tgt)=[];
    sgkn = zeros([ksz Q]);
    lambda = norm(smx,'fro')/size(smx,1)*tik;
    rmx = smx+eye(size(smx))*lambda;
    switch mth
        case 'cg'
            sgkn(fkmsk==1) = bicgstab(rmx, tmx, [], 100);
        case 'pinv'
            sgkn(fkmsk==1) = pinv(rmx)*tmx;
        case '\'
            sgkn(fkmsk==1) = rmx\tmx;
    end
    switch dim
        case 2
            kn(:,:,:,ii) = sgkn;
        case 3
            kn(:,:,:,:,ii) = sgkn;
    end
end

