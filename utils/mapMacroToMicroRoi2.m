function [indices,indicesrc ]= mapMacroToMicroRoi2(pixel,downsampleFactor,m,n)
%m,n size to which you want extrapolate, eg: pixel belongs to 32*32(8*8
%macro pixel) and you want to get indices for 64*64 then m,n = 64.
%similarly, downsampleFactor is on the basis of the target macro pixel
%resolution. so from 32 to 64 it would be 1.

vec = @(x) x(:);
aa = reshape(1:m*n,m,n);
% downsampleFactor = 3;
windowSize = 2.^downsampleFactor;
downsampleSize = 2.^(nextpow2(m) - downsampleFactor); %size of the macro pixel image from which to extrapolate or to which the windowSize belongs
% indices = [];
if rem(pixel,downsampleSize) == 0
    cl = pixel/downsampleSize - 1 ;
else
    cl = floor(pixel/downsampleSize);
end
if cl == 0
    colmn = 1:windowSize;
else
    colmn = windowSize*(cl-1)+1:windowSize*cl;
end
if pixel > downsampleSize
	if mod(pixel,downsampleSize) == 0
	rw = windowSize*(downsampleSize-1) + 1:windowSize*downsampleSize;	
else
    	rw = windowSize*(mod(pixel,downsampleSize)-1) +1:windowSize*mod(pixel,downsampleSize);
end
else
    rw = windowSize*(pixel- 1)+1:windowSize*pixel;
end

indices =vec(aa(rw,colmn));
indicesrc = vec([rw,colmn]);

