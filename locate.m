function ind = locate(arr, val)
%
% in = where(arr,val)
%
% Find the location of 'val' in an array 'arr', assuming arr is in
% ascending order
%
% (c) Yulin V Chang, 20150414
%

len = length(arr);
ind = 1;
for ii = 1:len
    if arr(ii)<val
        ind = ind+1;
    else
        break;
    end
end

