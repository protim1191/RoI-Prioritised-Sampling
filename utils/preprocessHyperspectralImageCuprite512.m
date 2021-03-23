function [yy,y] = preprocessHyperspectralImageCuprite512

load('Cuprite_f970619t01p02_r02_sc03.a.rfl.mat', 'X')

x = squeeze(X(:,1:512,:));
xx = reshape(x,prod([size(x,1),size(x,2)]),size(x,3));
xx = double(xx);

%% removing bands, information obtained from ENVI tuitorial for this particular dataset(cuprite)

bandsToRemove = [1,2,3,98:128,148:170];
xx(:,bandsToRemove) = [];

for ii = 1:size(xx,2) % normalize along the spectra
y(:,ii) = (xx(:,ii) - min(xx(:,ii)))./(max(xx(:,ii)) - min(xx(:,ii)));
end


yy = reshape(y,size(x,2),size(x,2),167);