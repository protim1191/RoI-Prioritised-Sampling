function [yy,y] = preprocessHyperspectralImageCuprite1
load('Cuprite_f970619t01p02_r02_sc03.a.rfl.mat', 'X')
x = squeeze(X(80:335,50:305,:));
% x = squeeze(X(:,1:512,:));
xx = reshape(x,prod([size(x,1),size(x,2)]),size(x,3));
%%
xx = double(xx);

%% removing bands, information obtained from ENVI tuitorial for this particular dataset(cuprite)
bandsToRemove = [1,2,3,98:128,148:170];
xx(:,bandsToRemove) = [];

% for ii = 1:size(xx,1) % normalize along the spectra
% y(ii,:) = (xx(ii,:) - min(xx(ii,:)))./(max(xx(ii,:)) - min(xx(ii,:)));
% end
for ii = 1:size(xx,2) % normalize along the spectra
    y(:,ii) = (xx(:,ii) - min(xx(:,ii)))./(max(xx(:,ii)) - min(xx(:,ii)));
end


yy = [];
% yy = reshape(y,size(xx,1),size(xx,2),168);
% for pseudocolor visulaization use bands 1,10 and 30 after band removal.