function kn = compkn(d,k,w,mdh,ksh,tik,mth)
%
% Compute 3D SPIRiT kernel for multiple series
%
%
% YV Chang, Jun 2015
%
%

nint = size(d,2);
npar = size(d,3);
Q    = size(d,4);
nmea = size(d,5);

crg  = mdh.PlnRef + .025;
csz  = round(mdh.MatrixResolutionReadout * crg);
clmx = crg*.5;
clin = locate(abs(k(:,1)),clmx);
kcl  = k(1:clin,:)*max(abs(k(end,:)))/max(abs(k(clin,:)));
wcl  = w(1:clin,:);
for ii = round(clin*.7):clin-10
    if (wcl(ii+9,1)-wcl(ii,1))/wcl(ii,1) > .0075
        break;
    end
end
for jj = 1:nint
    wcl(ii+6:end,jj) = wcl(ii+5,jj);
end
M    = [csz,csz];
NFFTc= NUFFT(kcl,wcl,[0 0],M);
acc  = mdh.ParAcc;
ncbl = mdh.RefBlc;
if isfield(mdh,'SlicePF') && (npar*mdh.SlicePF < npar/2+ncbl/2*acc+1)
    zcl = npar/2-ncbl/2*acc+1:npar*mdh.SlicePF;
else
    zcl = npar/2-ncbl/2*acc+1:npar/2+ncbl/2*acc+1;
end
cal  = d(1:clin,:,zcl,:,:);
dcl  = fft2vc(NFFTc'*(cal.*repmat(sqrt(wcl),[1,1,numel(zcl),Q,nmea])));
kn   = zeros([size(ksh) Q Q nmea]);
for ii = 1:nmea
    disp(['measurement = ' num2str(ii)]);
    kn(:,:,:,:,:,ii) = calibSPOT(dcl(:,:,:,:,ii), ksh, tik, mth);
end

