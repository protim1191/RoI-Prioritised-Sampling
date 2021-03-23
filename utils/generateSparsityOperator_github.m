function W = generateSparsityOperator_github(m,n,sparsityBasis,numLevels)
% generates operator for the chosen sparsity basis
% only possible options are DCT2, haar and Daubechies-8 wavelet. Can be expanded!
% Toolbox required: SPOT operator toolbox.
sparsityBasis = lower(sparsityBasis);
% numLevels = 3; % Number of levels in the wavelet decomposition, value of
% 3 used in the paper results.
switch (sparsityBasis)
    case 'dct'
        W = opDCT2(m,n);
    case 'haar'
        W = opWavelet2(m, n, 'haar',[],numLevels);
    case 'daubechies'
        W = opWavelet2(m, n, 'daubechies',8,numLevels);
    otherwise
        error('Unknown sparsity operator. Possibilities DCT, haar or daubechies')
end
