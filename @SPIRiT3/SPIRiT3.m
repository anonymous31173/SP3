function sp3 = SPIRiT3(kn3, mth, msz)
%
% spt = SPIRiT3(kernel [,method, imSize])
%
% Construction of 3D SPIRiT operator
%
% Inputs:
%   kn3: [kx,ky,kz,Q,Q], 3D SPIRiT convolution kernel
%   mth: method
%   msz: image size
%
% Output:
%   spt: 3D SPIRiT kernel
%
% Yulin Chang, 20150614
%

if nargin < 3
    mth = 'conv';
    K3  = [];
end

if strcmp(mth,'conv')==1
    K3  = [];
end

if strcmp(mth, 'fft')==1 || strcmp(mth, 'image')==1
    Q = size(kn3,4);
    K3 = ifft3vc(zpad(kn3(end:-1:1,end:-1:1,end:-1:1,:,:)*sqrt(prod(msz)),msz(1), msz(2), msz(3), Q, Q));
end

sp3.kernel  = kn3;
sp3.adjoint = 0;
sp3.K3      = K3;
sp3.mth     = mth;
sp3.msz     = msz;
sp3         = class(sp3,'SPIRiT3');

