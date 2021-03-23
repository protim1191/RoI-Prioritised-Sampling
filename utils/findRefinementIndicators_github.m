function [V,B,yB] = findRefinementIndicators_github(uC, org, ws, extraMeas, ensemble,ndev)
% find refinement indicator of  region
% Inputs:
% uC : Coarse solution of ROI at current resolution level, has to be 2D
% org : Original ROI has to be 2D
% ws : window size of the next resolution level
% extraMeas : percent of extra measurements based on size of roi
% Outputs:
% V : refinement indicator value for uC
% B : the measurement matrix for refinement
% fB: new measurements
% uses opSamplingMatrix from utils


% nk = floor(extraMeas*numel(org));

B = opSamplingMatrix_github(ws, extraMeas, size(org,1),size(org,2),rng, ensemble);
yB = B*org(:)+ ndev*randn(size(B,1),1);
V = sum(abs(B*uC(:) - yB).^2);
end

    