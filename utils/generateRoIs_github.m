function seg = generateRoIs_github(lowResRecon, macroPixelSize, m, n, maxRegions, tau, radius, original, ensemble)
% generates RoIs from low resolution image.
% Inputs:

minArgs=8;
maxArgs=9;
narginchk(minArgs,maxArgs)
[~, loc, seg.segOrg] = getSegmentedImage(lowResRecon, macroPixelSize, m, n, maxRegions, tau, radius);
fprintf('Number of RoIs after merging = %d\n',length(loc{1,1}))
% figure(); imshow(data.segOrg)

seg.ROI = cell(length(loc{1,1}),2);
seg.maskC = cell(1,length(loc{1,1}));
seg.opR = cell(1,length(loc{1,1}));
seg.roim = zeros(1,length(loc{1,1}));
seg.roin = zeros(1,length(loc{1,1}));
seg.rstrt = zeros(1,length(loc{1,1}));
seg.cstrt = zeros(1,length(loc{1,1}));
seg.segROI = zeros(m,n);
% Below are flags to check RoI overlap
seg.expandedRoiIndices = cell(1,length(loc{1,1}));
seg.maxRegions = maxRegions;
seg.tau = tau;
seg.radius = radius;
count = 1;
removeFlag = 0;
for kk = 1:length(loc{1,1})
   
    region = loc{1,1}{1,kk};
    [seg.opR{1,kk},seg.roim(1,kk),seg.roin(1,kk),seg.maskC{1,kk},...
        seg.rstrt(1,kk), rstp, seg.cstrt(1,kk), cstp] = ...
        getROI2(region, macroPixelSize, m, n);
    
    % The next few lines check whether regions overlap. If 99% of the
    % pixels overlap they are merged. This is possible because getROI2 will
    % expand to the next dyadic size.
    seg.expandedRoiIndices{1,kk} = find(~seg.maskC{1,kk});
    
    if kk >1
        for jj = kk-1:-1:1
            ss = sum(1.*(ismember(seg.expandedRoiIndices{1,kk},seg.expandedRoiIndices{jj})));
            if ss >= 0.99*length(seg.expandedRoiIndices{1,kk})
                seg.expandedRoiIndices{1,kk}=[];
                fprintf('RoI number %d is a subset of RoI %d, Removing RoI %d \n', kk, jj, kk)
                removeFlag = 1;
                break
                
            end
        end
    end
    if removeFlag
        removeFlag = 0;
    else
        seg.ROI{count,1} = reshape(seg.opR{1,kk}*original(:),...
            seg.roim(1,kk),seg.roin(1,kk)); %original ROI
        if strcmp(ensemble,'Walsh')
            seg.ROI{count,2} = 1;
        else
            seg.ROI{count,2} = macroPixelSize;
        end
        seg.segROI(seg.rstrt(1,kk):rstp,seg.cstrt(1,kk):cstp) = kk;

    count = count + 1;
   
    end
   
end
emptyCells = cellfun(@isempty,seg.ROI);
seg.ROI(emptyCells)  =[];
seg.numRoIs = count-1;
seg.ROI = reshape(seg.ROI,seg.numRoIs,2);