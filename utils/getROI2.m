function [opR,roim,roin,maskC,rstrt,rstp,cstrt,cstp] = getROI2(region, windowSize, m, n)
% generates the ROI from locations given in region. Returns the restriction
% operator for the ROI and the size of the ROI

vec = @(x) x(:);
downsampleFactor = nextpow2(windowSize);
for ii = 1:length(region)
    [~,indicesrc(:,ii)] = mapMacroToMicroRoi2(region(ii),downsampleFactor,m,n);
end
%%
%getting original row and column start and stop value. need to modify them
%to make the size of the ROI to be a power of 2.
rstrtOld = min(vec(indicesrc(1:size(indicesrc,1)/2,:))); rstpOld = max(vec(indicesrc(1:size(indicesrc,1)/2,:)));
cstrtOld = min(vec(indicesrc(size(indicesrc,1)/2 + 1:end,:))); cstpOld = max(vec(indicesrc(size(indicesrc,1)/2 + 1:end,:)));
% calculate original size of the ROI and how much it differs from the
% nearest power of 2.
roimOld = abs(rstpOld - rstrtOld) + 1;
roinOld = abs(cstpOld - cstrtOld) + 1;
%ppr is the next largest index of 2 for the row dimension
ppr  = nextpow2(roimOld); ppc = nextpow2(roinOld);
%ensuring ppr and ppc have the same value. We want a square ROI.
if ppr < ppc
    ppr = ppc;
else
    ppc = ppr;
end

ppr1 = 2.^ppr - roimOld ; ppc1 = 2.^ppc - roinOld ;
rowAd1 = floor(ppr1/2);rowAd2 = ppr1 - rowAd1;
colAd1 = floor(ppc1/2); colAd2 = ppc1 - colAd1;

flg = [0,0,0,0];

if rstrtOld - rowAd1 <= 0
    rstrt = rstrtOld;
    flg = flg + [1,0,0,0];
else
    rstrt = rstrtOld - rowAd1;
end
if rstpOld + rowAd2 >= m
    rstp  = rstpOld;
    flg = flg + [0,1,0,0];
else
    rstp = rstpOld + rowAd2;
end
if cstrtOld - colAd1 <= 0
    cstrt = cstrtOld;
    flg = flg + [0,0,1,0];
else
    cstrt = cstrtOld - colAd1;
end

if cstpOld +colAd2 >= n
    cstp = cstpOld;
    flg = flg + [0,0,0,1];
else
    cstp = cstpOld + colAd2;
end
% fprintf('%d \n',flg)
id = find(flg == 1);
if isempty(id)
else
    for ii = 1: length(id)
        ix = id(ii);
        switch(ix) % flag to calculate the boundary cases, i.e., if ctsp ==m etc
            case 1
                rstrt = 1;
                rstp = rstp + abs(rstrtOld - rowAd1) +1;
            case 2
                rstp = m;
                rstrt = rstrt - abs(rowAd2 - (m - rstpOld));
            case 3
                cstrt = 1;
                cstp = cstp + abs(cstrtOld - colAd1)+1;
            case 4
                cstp = n;
                cstrt = cstrt - abs(colAd2 - (n-cstpOld));
                
        end
    end
end

if rstrt <= 0 || rstp > m || cstrt <= 0 || cstp >n
    error('ROI size problem')
end


roim =rstp - rstrt + 1; roin = cstp - cstrt + 1;
mask = zeros(m,n); mask(rstrt:rstp,cstrt:cstp) = 1;mask = logical(mask);maskC = not(mask);
opR = opRestriction(m*n,find(mask==1));
