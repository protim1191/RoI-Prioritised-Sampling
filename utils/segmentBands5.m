function [segmentedImage,locations,regionsMerged] = segmentBands5(data,maxRegions,tau,mergeRadius)
%segment and merge
% created new function mergeCheck to ensure pixels at the edges are also
% merged if possible.
% 28-01-2020
[m,n,p] = size(data);
vec = @(x) x(:);
contrast = @(x) ((max(x) - min(x))./(max(x) + min(x)));
loc = cell(1,maxRegions);
count = 1;
regionsMerged = [];
for pp = 1:p
    img = data(:,:,pp);
    [~,idx] = sort(img(:),'descend');
    for kk = 1:maxRegions
        if isempty(idx)
            break
        end
        id = idx(1);
        [row,col] = ind2sub([m,n],id); % do calculations with linear indices
        if row == m || col == n || row == 1 || col == 1
            fprintf('Center near the edge\n')
            loc{1,kk} = idx(1);
            idx(1) = [];
            if mergeRadius> 0
                if kk>1
                [regionsMerged, count] = mergeCheck(m,n,kk,mergeRadius,id,count, regionsMerged,loc);
                end
            end
            continue
        end
        ii = 1;
        neighbourhoodOld = getNeighbourIndices2(m,n,ii,id);
        c = contrast(img(neighbourhoodOld));
        while (c<tau)
            ii = ii + 1;
            neighbourhood = getNeighbourIndices2(m,n,ii,id);
            if neighbourhood(ii) > m*n || neighbourhood(2*ii) < 0 || mod(neighbourhood(3*ii),m)==0|| mod(neighbourhood(4*ii),m)==1
                fprintf('ROI exceeds Image boundaries \n')
                ii = ii - 1;
                break
            else
                c = contrast(img(neighbourhood));
                neighbourhoodOld = neighbourhood;
            end
        end
        loc{1,kk} = neighbourhoodOld;
        %        boundaryPixels{1,kk} = getROIBoundary(m,ii,neighbourhoodOld);
        %% merge
        if mergeRadius> 0
            if kk >1
            [regionsMerged, count] = mergeCheck(m,n,kk,ii+mergeRadius,id,count, regionsMerged, loc);
            end
        end
        [~,idIdx, ~] = intersect(idx,loc{1,kk});
        idx(idIdx) = [];
        
    end
end

%%
% region merging
if ~isempty(regionsMerged)
    newRegions = mergeRegions(regionsMerged);
    ff = checkIntersection(newRegions);
    while ff
        newRegions = mergeRegions(newRegions);
        ff = checkIntersection(newRegions);
    end
    bb = setdiff(1:maxRegions, [newRegions{:}]);
    nn = length(newRegions);
    for ii = 1:length(bb)
        newRegions{1,nn + ii} = bb(ii);
    end
    
    %%
    % merge region indices
    
    newLoc = cell(1,size(newRegions,2));
    for ii = 1:length(newLoc)
        for jj = 1:length(newRegions{1,ii})
            newLoc{1,ii} = [newLoc{1,ii} loc{1,newRegions{1,ii}(jj)}'];
        end
    end
else
    newLoc = loc;
end
%%
seg = zeros(m,n);
for kk = 1:length(newLoc)
    %     seg(newLoc{1,kk}) = 1;
    seg(newLoc{1,kk}) = kk;
end
segmentedImage(:,:,pp) = seg;
locations{1,pp} = newLoc;
end

%%
function mergedRegions = mergeRegions(regionsMerged)
%provides indices for new regions after they have been merged
newRegionCount = 1;
mergedRegions{1,newRegionCount} = regionsMerged{1,:};
for ii = 2:size(regionsMerged,2)
    independent = 0;
    for jj = 1:(newRegionCount+1)-1
        if ~isempty(intersect(regionsMerged{1,ii},mergedRegions{1,jj}))
            mergedRegions{1,jj} = unique([regionsMerged{1,ii},mergedRegions{1,jj}]);
        else
            independent = independent + 1;
        end
    end
    if independent == newRegionCount
        newRegionCount = newRegionCount +1 ;
        mergedRegions{1,newRegionCount} = regionsMerged{1,ii};
    end
end
end
function fl = checkIntersection(region)
fl = 0;
for ii = 1:size(region,2)
    for jj = ii +1 : size(region,2)
        if ~isempty(intersect(region{1,ii},region{1,jj}))
            fl = 1;
            break
        end
    end
end
end
function [regionsToBeMerged, count] = mergeCheck(m,n,kk,radius,id,count, regionsToBeMerged,loc)
%check whether it is possible to merge current region with previous
%regions, if yes record them



    searchAreaIndices = getNeighbourIndices2(m,n,radius,id); % get the search area around the current region where to find other possible ROIs.
    
    for jj = 1:kk-1
        fl(:,jj) = isempty(intersect(loc{1,jj},searchAreaIndices));
    end
    
    for jj = 1:length(fl)
        if ~fl(jj)
            fprintf('merging regions %d and %d \n',jj,kk)
            %                     regionsMerged(count,:) =[jj,kk];
            regionsToBeMerged{1,count} = [jj,kk];
            count = count + 1;
        end
    end

end
