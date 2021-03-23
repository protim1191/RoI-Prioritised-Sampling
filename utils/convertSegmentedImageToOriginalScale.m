function segOrg = convertSegmentedImageToOriginalScale(loc,downsampleFactor,m,n)

count = 0;
segOrg = zeros(m,n);
for ii = 1:length(loc{1,:})
    region = loc{1,:}{ii};
    for jj = 1:length(region)
        count = count + 1;
%         indices(:,count)= mapMacroToMicroRoi3(region(jj),downsampleFactor,m,n);
%  mapMacroToMicroRoi3 not working properly
        indices(:,count)= mapMacroToMicroRoi2(region(jj),downsampleFactor,m,n);
    end
    count = 0;
    segOrg(indices(:))=ii;
end


% figure(); imshow(segOrg)