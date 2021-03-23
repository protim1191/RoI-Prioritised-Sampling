function mask = drawSamples(samples,m,chk)
% samples need to a vector containing the indices of the sampled
% points
% m is the image size,scalar, 64 for a 64 x 64 image
coeffGrid = generateGrid(m);
mask = zeros(m);
for ii = 1:length(samples)
    id = find(samples(ii)==coeffGrid(:));
    mask(id) = 1;
end
if chk == 1
figure(); imshow(mask)
figure(); imagesc(mask)
end