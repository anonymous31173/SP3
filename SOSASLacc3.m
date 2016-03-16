function [im, RV] = SOSASLacc3(d,k,w,mdh,kn,nit,lmd,mt0)
%
% im = SOSASL3dl1(d,k,w,mdh,ksz,tik,nit,lmd,lev,mth,thr,im0);
%
% Input
%   d   -- k-space data. d should be entirely in k-space
%   k   -- k-space trajectory
%   w   -- weights
%   mdh -- data info
%   tik -- Tikhonov regularization
%   nit -- number of iterations
%   lmd -- ratio between data consistency and calibration consistency
%   mt0 -- Method of SPIRiT3 operator construction
%
% Output
%   im  -- image with    parallel imaging reconstruction
%   RV  -- RESVEC, the residual output from lsqr
%
% Yulin V Chang, 20150715
%

if nargin < 8
    mt0 = 'image';
end

N = [mdh.MatrixResolutionReadout,mdh.MatrixResolutionReadout];
npar = size(d,3);
Q    = size(d,4);
nmea = size(d,5);

NFFT = NUFFT(k,w,[0 0],N);

aid  = reshape(any(permute(d,[5 1 2 3 4])),[],size(d,2),npar,Q);
im   = zeros([N npar Q nmea]);
RV   = zeros(nit+1,nmea);
for ii = 1:nmea
    disp(['measurement = ' num2str(ii)]);
    switch mt0
        case 'image'
            m0 = NFFT'*(ifftc(d(:,:,:,:,ii),3).*repmat(sqrt(w),[1 1 npar Q]));
        case 'conv'
            m0 = fft2vc(NFFT'*(d(:,:,:,:,ii).*repmat(sqrt(w),[1 1 npar Q])));
    end
    disp('generating SPIRiT operator...');
    switch mt0
        case 'conv'
            GOP = SPIRiT3(kn(:,:,:,:,:,ii),'conv',[N npar]);
        case 'image'
            GOP = SPIRiT3(kn(:,:,:,:,:,ii), 'image', [N npar]);
    end
    disp('Reconstructing images...');
    [im(:,:,:,:,ii),RV(:,ii)] = cgNUSPIRiT3(d(:,:,:,:,ii).*repmat(sqrt(w),[1 1 npar Q]),aid,NFFT,GOP,nit,lmd,m0);
end

