function [Walsh_M, reqdMeas] = generateWalshMeasurementMatrix_github(m,n, numLevel, currNumLevel, remMeas)
if currNumLevel > numLevel
    error('Requested decomposition level exceeds desinged levels for Walsh transform not available.\n')
end
if m ~=n
    error('RoI is not a square')
end
switch(m)
    case 32
        om = load('chosenSamplingScheme32percent20_100mc.mat');
        omega = om.chosenSamplingScheme;
    case 64
        om = load('chosenSamplingScheme64percent20_20mc.mat');
        omega = om.chosenSamplingScheme;
    case 128
        om = load('chosenSamplingScheme128percent20_20mc.mat');
        omega = om.chosenSamplingScheme;
    otherwise
        if roiSize < 32
            fprintf('RoI size too small, Pattern not generated for %d x %d\n', roiSize, roiSize)
            omega = [];
        else
            fprintf('sampling pattern not generated for size %d x %d \n', roiSize, roiSize)
            omega = [];
        end
end

omegaCurrLevel = samplesInCurrentLevel(m,currNumLevel,omega);
%%
%mask = drawSamples(omega,m,0); %requires by separateWaveletBands for divided frequency space in dyadic levels
%[~,~,subBandIds] = separateWaveletBands(mask,numLevel);

% t=[];
% find if omega is present in current G,i.e. subBandIds{1,ws}
% for ii = 1:length(omega)
%     tt = ismember(omega(ii),subBandIds{1,currNumLevel} );
%     if tt
%         t = [t; omega(ii)];
%     end
% end
% omega = t;
%%
if remMeas >= length(omegaCurrLevel)
    R = opRestriction(m*m,omegaCurrLevel);
    walsh = generateWalshOperator_unNormalized(m*n);
    nrmFac = 1/sqrt(length(omegaCurrLevel));
    nrmFac = opDiag(opOnes(m*n,1)*nrmFac);
    walsh_nrm = opFoG(nrmFac,walsh);
    Walsh_M = opFoG(R,walsh_nrm);
    reqdMeas = length(omegaCurrLevel);
else
    error('Enough Measurements not available for for Walsh measurements at Level = %d\n',currNumLevel);
end

function omegaCurrLevel = samplesInCurrentLevel(m,numLevel,omega)

switch(m)
    case 128
        G128 = generateGrid(128);
        switch(numLevel)
            case 1
                subBandIds = vec(G128(1:45,1:45));
            case 2
                subBandIds = [vec(G128(46:64,1:64)); vec(G128(1:45,46:64))];
            case 3
                subBandIds = [vec(G128(65:128,1:128)); vec(G128(1:64,65:128))];
        end
    case 64
        G64 = generateGrid(64);
        switch(numLevel)
            case 1
                subBandIds = vec(G64(1:26,1:26));
            case 2
                subBandIds = [vec(G64(27:32,1:32)); vec(G64(1:26,27:32))];
            case 3
                subBandIds = [vec(G64(33:64,1:64)); vec(G64(1:32,33:64))];
        end
    case 32
        G32 = generateGrid(32);
        switch(numLevel)
            case 1
                subBandIds = vec(G32(1:11,1:11));
            case 2
                subBandIds = [vec(G32(12:16,1:16)); vec(G32(1:11,12:16))];
            case 3
                subBandIds = [vec(G32(17:32,1:32)); vec(G32(1:16,17:32))];
        end
        
        
end
omegaCurrLevel = [];
for ii = 1:length(omega)
    tt = ismember(omega(ii),subBandIds);
    if tt
        omegaCurrLevel = [omegaCurrLevel; omega(ii)];
    end
end
