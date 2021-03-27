function [V,B,yB,ms] = findRefinementIndicatorsWalsh_github(uC, org, numLevels, currNumLevels,remMeas,ndev)
% find refinement indicator of  region
% Inputs:
% uC : Coarse solution of ROI at current resolution level, has to be 2D
% org : Original ROI has to be 2D
% numLevels : Total number of decompositions level, a value of 3 is used througuout. 
% currNumLevels: The level at which refined solution is already calculated.
% remMeas: Remaining Measurements
% ndev: Noise std
% 
% Outputs:
% V : refinement indicator value for uC
% B : the measurement matrix for refinement
% yB: new measurements
% ms: Remaining measurements after calculation of RI
% uses opSamplingMatrix from utils
% Protim 2 April 2020
if currNumLevels > numLevels
    V = 0;
    yB = [];
    B = [];
    ms = 0;
    fprintf('RoI in the native resolution of sensor. \n')
    return
end
[m,n] = size(uC);
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
        if m < 32
            fprintf('RoI size too small, Pattern not generated for %d x %d\n', m, m)
            omega = [];
        else
            fprintf('sampling pattern not generated for size %d x %d \n', m,m)
            omega = [];
        end
end
omegaCurrLevel = samplesInCurrentLevel(m,currNumLevels,omega);

if remMeas >= length(omegaCurrLevel)
    ms = length(omegaCurrLevel);
    R = opRestriction(m*n,omegaCurrLevel);
    walsh = generateWalshOperator_unNormalized(m*n);
    nrmFac = 1/sqrt(length(omegaCurrLevel));
    nrmFac = opDiag(opOnes(m*n,1)*nrmFac);
    walsh_nrm = opFoG(nrmFac,walsh);
    B = opFoG(R,walsh_nrm);
    yB = B*org(:) + ndev*randn(size(B,1),1);
    V = sum(abs(B*uC(:) - yB).^2);
 else
    V = 0;
    yB = [];
    B =[];
    ms = 0;
    fprintf('Not enough measurements for current ROI. \n')
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
