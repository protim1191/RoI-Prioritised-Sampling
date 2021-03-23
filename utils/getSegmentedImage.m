function [seg, loc, segOrg] = getSegmentedImage(xr1, windowSize, m, n, maxRegions, tau, seopt)

vec = @(x) x(:);
pp = nextpow2(m) - nextpow2(windowSize);
img = zeros(2.^pp);
kk = 0; ll = 0;
for ii = 1:windowSize:m
    kk = kk + 1;
    for jj = 1:windowSize:n
        ll = ll +1;
        img(kk,ll) = mean(vec(xr1(ii:ii+(windowSize-1),jj:jj+(windowSize-1))));
    end
    ll = 0;
end

downsampleFactor = nextpow2(windowSize);
% [seg,loc] = segmentBands3(img, maxRegions, tau, seopt);% input data has to be 2D
[seg,loc] = segmentBands5(img, maxRegions, tau, seopt);% input data has to be 2D
 segOrg = convertSegmentedImageToOriginalScale(loc,downsampleFactor,m,n);
% segOrg = convertSegmentedImageToOriginalScale2(loc, windowSize,m,n);

