close all
clear
clc
% add utils folder
addpath(genpath('C:\Users\Eigenaar\Documents\MATLAB\codes\RPS\utils\'))

%%  user data and initializations
dataset = 'cuprite';
imSize = 256;
[data.original,init] = loadDataset_github(dataset,imSize);
init.noiseDev = 1e-7;
%%
rps.totalMeasurements = 5300;
rps.lowRes.macroPixelSize = 8; %macro pixel size for low resolution reconstruction(step 1).
rps.lowRes.meas = 1000; % num meas for low resolution recon (step 1).
rps.seg.maxRegions = 10; % num. of RoIs in the segmentation (step 2).
rps.seg.tau = 0.1756; % contrast threshold for cut-off.
rps.seg.radius = 1; % radius for merging RoIs.

rps.lowRes.sparsityBasis = 'dct'; % 'dct', 'daubechies','haar', none
rps.low.sparsityDecompLevels = [];
rps.lowRes.solver = 1; % 1= l1, 2 = l2, only for first level

% rps.ri.MeasurementEnsemble = 'random'; % random 0/1, rademacher, walsh
% rps.ri.numRefinedMeas = 0.1;
% given in fraction between 0.01 and 1. 10% of num pixels in RoI = 0.1.
%for walsh measurements, it is a redundant paramter. It is fixed to 20%.
rps.ri.numDecompLevels_Walsh = 3; % number of levels for walsh decomposition.
rps.ri.sparsityBasis = 'daubechies';
rps.ri.sparsityDecompLevels = 3;
gg = 1;
cmap = hsv(rps.seg.maxRegions);
%% tfocs options.
tfocs.opts.maxIts = 20000; tfocs.opts.printEvery = 0;
tfocs.opts.continuation = 1; tfocs.opts.alg = 'N07';
tfocs.EPS = 0.1*init.noiseDev;
tfocs.alpha = 1; tfocs.beta = 0.6;
% for coarse reconstructions only 1000 iterations were used to generate
% results in the paper.
savestr = strcat(dataset,'_',num2str(imSize),'_Walsh_',datestr(now,'dd-mm-yyyy'),'_',datestr(now,'HH_MM'),'.mat');
clearvars dataset imSize

%% step 1: Low resolution reconstruction
load('seedToBeUsedForRandomMatrixGeneration.mat')
rps.lowRes.M = opSamplingMatrix_github(rps.lowRes.macroPixelSize,rps.lowRes.meas,init.m,init.n,seed,'random');
fprintf('Generated low resolution sampling matrix with macro pixel size = %d \n', rps.lowRes.macroPixelSize)
meas.lowRes = rps.lowRes.M*data.original(:) + init.noiseDev*randn(size(rps.lowRes.M,1),1);
fprintf('Generated observations \n')

if rps.lowRes.solver == 1
    rps.lowRes.W = generateSparsityOperator_github(init.m,init.n,rps.lowRes.sparsityBasis,rps.low.sparsityDecompLevels);
    tfocs.mu = 0.5*norm(rps.lowRes.W*data.original(:),2)/norm(data.original(:),2);
    
else
    rps.lowRes.W=[];
end

[results.lowResRecon] = coarseSamplingReconTFoCs_github(rps.lowRes.M,...
    rps.lowRes.solver, rps.lowRes.W,[init.m, init.n], tfocs, data.original,init.noiseDev,0);

results.addedImage = results.lowResRecon;
figure(); imshow(results.lowResRecon); title('Step 1: Low resolution sampling and reconstruction')
str = strcat('Macro pixel size = ', num2str(rps.lowRes.macroPixelSize),' solver = l', num2str(rps.lowRes.solver));
xlabel(str)
meas.remainingMeasurements = rps.totalMeasurements - rps.lowRes.meas;
fprintf('Remaining measurements after low resolution reconstruction = %d \n',meas.remainingMeasurements)

% step 2: Detection and Segmentation
rps.seg = generateRoIs_github(results.lowResRecon, rps.lowRes.macroPixelSize, init.m, init.n, rps.seg.maxRegions, rps.seg.tau, rps.seg.radius,data.original,'Walsh');
figure(); imshow(rps.seg.segROI,[])
title('Step 2: Detection and segmentation')
xlabel('Different colors are masks for different RoIs')
%% step 3: Prioritisation of RoIs with RI
rois_currentRes= cell(1,size(rps.seg.ROI,1));
count = 1;
for kk = 1:size(rps.seg.ROI,1)
    rois_currentRes{1,kk} = reshape(rps.seg.opR{1,kk}*vec(results.lowResRecon),rps.seg.roim(1,kk), rps.seg.roim(1,kk));
    results.NMSE{kk,count} = nmse(rois_currentRes{1,kk},rps.seg.ROI{kk,1});
    results.SSIM{kk,count} = ssim(rois_currentRes{1,kk},rps.seg.ROI{kk,1});
    count = count + 1;
end
V = zeros(1,size(rps.seg.ROI,1)); B = cell(1,size(rps.seg.ROI,1));
for kk = 1:size(rps.seg.ROI,1)
    [V(kk),B,yB,ms] = findRefinementIndicatorsWalsh_github(rois_currentRes{1,kk},rps.seg.ROI{kk,1},...
        rps.ri.numDecompLevels_Walsh, rps.seg.ROI{kk,2}, meas.remainingMeasurements,init.noiseDev);
    meas.ri.y{1,kk} = yB;
    rps.ri.A{1,kk} = B;
    meas.remainingMeasurements = meas.remainingMeasurements - ms;
end
fprintf('remaining measurements after initial RI calculations = %d \n', meas.remainingMeasurements)
results.RI{1,gg} = V;
gg = gg + 1;
clearvars B ms yB kk

%%

while meas.remainingMeasurements >=0
    [~,kk] = max(V);
    if sum(V) == 0
        fprintf('All acquired RoIs are at the highest possible resolution \n')
        break
    end
    % calculate the refined solution of the selected RoI.
    rois_currentRes{1,kk} = reconTfocsWaveletTV_github(rps.seg.roim(1,kk),rps.seg.roin(1,kk),...
        rps.ri.A{1,kk},rps.ri.sparsityBasis,rps.ri.sparsityDecompLevels,init.noiseDev,rps.seg.ROI{kk,1},tfocs,1,meas.ri.y{1,kk});
    results.ri.xr{kk,count} = rois_currentRes{1,kk};
    rps.seg.ROI{kk,2} = rps.seg.ROI{kk,2} + 1;
    %     calculate the refinement indicator for the next fine scale for the
    %     selected RoI, all other RoIs carry forward their values
    %     measReq = floor(rps.ri.numRefinedMeas*numel(rps.seg.ROI{kk,1}));
    %     if  meas.remainingMeasurements > measReq && ~(rps.seg.ROI{kk,2} == 1)
    %         [V(kk), B{1,kk}, yB] = findRefinementIndicators_github(rois_currentRes{1,kk}, rps.seg.ROI{kk,1},...
    %             rps.seg.ROI{kk,2}/2, measReq, rps.ri.MeasurementEnsemble, init.noiseDev);
    [V(kk),B,yB,ms] = findRefinementIndicatorsWalsh_github(rois_currentRes{1,kk},rps.seg.ROI{kk,1},...
        rps.ri.numDecompLevels_Walsh, rps.seg.ROI{kk,2}, meas.remainingMeasurements,init.noiseDev);
    meas.ri.y{1,kk} = [meas.ri.y{1,kk}; yB];
    rps.ri.A{1,kk} = [rps.ri.A{1,kk}; B];
    meas.remainingMeasurements = meas.remainingMeasurements - ms;
    %     else
    %         if meas.remainingMeasurements < measReq
    %             fprintf('Not enough measurements available for ROI %d. Current refined resolution level %d\n',kk,rps.seg.ROI{kk,2})
    %         elseif ROI{kk,2} == 1
    %             fprintf('Current refined resolution level for ROI %d is already %d\n',kk,rps.seg.ROI{kk,2})
    %         end
    %         V(kk) = 0;
    %     end
    results.RI{1,gg} = V;
    gg = gg + 1;
    results.addedImage = rps.seg.maskC{1,kk}.*results.addedImage + reshape(rps.seg.opR{1,kk}'*rois_currentRes{1,kk}(:),init.m,init.n);
    imshow(results.addedImage);hold on
    rectangle('position',[rps.seg.cstrt(1,kk),rps.seg.rstrt(1,kk),...
        rps.seg.roim(1,kk),rps.seg.roin(1,kk)],'EdgeColor', cmap(kk,:), 'LineWidth', 3)
    drawnow
    results.NMSE{kk,count} = nmse(rois_currentRes{1,kk},rps.seg.ROI{kk,1});
    results.SSIM{kk,count} = ssim(rois_currentRes{1,kk},rps.seg.ROI{kk,1});
    count = count + 1;
    fprintf('remaining measurements = %d \n',meas.remainingMeasurements)
end
% stpTime = toc;
% fprintf('Time required to complete the RPS procedure = %d s', strpTime - strtTime);
indicator = cell2mat(results.RI')';
figure(); imshow(results.addedImage); title('Final multi-resolution image');hold on
for kk = 1:length(rps.seg.roim)
    if isempty(rps.seg.expandedRoiIndices{1,kk})
        continue
    end
    rectangle('position',[rps.seg.cstrt(1,kk),rps.seg.rstrt(1,kk),rps.seg.roim(1,kk),rps.seg.roin(1,kk)],...
        'EdgeColor', cmap(kk,:), 'LineWidth', 3)
end
save(savestr,'-v7.3')

