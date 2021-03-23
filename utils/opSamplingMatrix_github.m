function M = opSamplingMatrix_github(windowSize,meas,m,n,seed, ensemble)
% generate a binary measurement matrix with given macro pixel size. Returns
% SPOT operator. both random 0/1 and rademacher systems can be generated.

rng(seed)
sensor = zeros(m,n);
M = zeros(meas,m*n);
for kk = 1:size(M,1)
    
    for ii = 1:windowSize:m - (windowSize -1)
        for jj = 1:windowSize:n-(windowSize -1)
            switch(ensemble)
                case 'random'
                    sensor(ii:ii+(windowSize-1),jj:jj+(windowSize-1)) = randi([0,1],1);
                case 'rademacher'
                    t = randn(1)<0;
                    if t
                        sensor(ii:ii+(windowSize-1),jj:jj+(windowSize-1)) = -1;
                    else
                        sensor(ii:ii+(windowSize-1),jj:jj+(windowSize-1)) = 1;
                    end
            end
        end
    end
    M(kk,:) = reshape(sensor,1,m*n);
end
for ii = 1:size(M,2)
    M(:,ii) = M(:,ii)./norm(M(:,ii));
end
M = opMatrix(M);