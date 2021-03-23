function walsh = generateWalshOperator_unNormalized(N)

% needs SPOT toolbox

h =  opHadamard(N);
HadIdx = 0:N-1;                             % Hadamard index
M = log2(N) + 1;                              % Number of bits to represent the index
binHadIdx = fliplr(dec2bin(HadIdx,M))-'0';  % Bit reversing of the binary index
binSeqIdx = zeros(N,M-1);                   % Pre-allocate memory
for k = M:-1:2
    binSeqIdx(:,k) = xor(binHadIdx(:,k),binHadIdx(:,k-1)); % Binary sequency index
end
SeqIdx = binSeqIdx*pow2((M-1:-1:0)');
walsh = h(SeqIdx+1,:);



end