function [w, n, kk] = yc_DCFac(k)
%
% [w, n, kk] = yc_DCFac(k);
%
% Calculation of density compensation function for k-space trajectories of
% Archimedean spirals
%
% Input:
%   k  -- k-space trajectory
%
% Output:
%   w  -- density compensation function
%   n  -- image size
%   kk -- normalized k-space trajectory within 1
%
% Notes:
%
% Yulin V. Chang, University of Pennsylvania, 20140918
%



[np, ni] = size(k); 
[theta, rho] = cart2pol(real(k(:,1)),imag(k(:,1)));
if ni>1
    rhtmp = angle(k(100,:));
    dif = round((rhtmp(2)-rhtmp(1))/(2*pi/ni));
    if dif==1 || dif==-(ni-1)
        rot = 1;
    elseif dif==-1 || dif==ni-1
        rot = -1;
    else
        error('In what sense does the spiral rotate? Are these regular spirals? Please check');
    end
end
nt = 110;
th = theta(end-nt+1:end);
rh = rho(end-nt+1:end);
kmax = rho(end);
[mth,ith] = min(th);
if mth<-3
    th(ith:end) = th(ith:end)+2*pi;
end
P = polyfit(rh,th,1);
arc = P(1);
off = P(2);
delk = 2*pi/arc;
n = round(4096*kmax*ni/2047/delk);

% Now we extend the current k-space to make density compensation easier to compute
P = polyfit((1:nt)',rh,1);
nxtp = round((2*pi/ni+pi/6)/P(1)/arc);
kxtr = zeros(nxtp, ni);
xtp = (nt+1:nt+nxtp)';
rhxtp = P(1)*xtp + P(2);
thxtp = rhxtp*arc + off;
[xx,yy] = pol2cart(thxtp, rhxtp);
kxtr(:,1) = xx+1i*yy;
if ni>1
    for ii = 2:ni
        kxtr(:,ii) = kxtr(:,1)*exp(2*rot*1i*(ii-1)*pi/ni);
    end
end

nk = zeros(np+nxtp,ni);
nk(1:np,:) = k;
nk(np+1:end,:) = kxtr;

w = voronoidens(nk);
w = w(1:np,:);

kk = k*2047/kmax/4096;



