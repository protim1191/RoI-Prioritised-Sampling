function [yy,y] = preprocessHyperspectralImageMexicoGulf1
load('mexicoGulfHyperSpectralData.mat')
x = squeeze(x1(180:435,50:305,:));
xx = reshape(x,prod([size(x,1),size(x,2)]),size(x,3));
%%
xx = double(xx);
%% removing bands, 
% removing the same bands cuprite using the ENVI file data for gulf data.
% No gurantee whether it is correct or not.

bandsToRemove = [169:218,248:279];
xx(:,bandsToRemove) = [];

% for ii = 1:size(xx,1) % normalize along the spectra
% y(ii,:) = (xx(ii,:) - min(xx(ii,:)))./(max(xx(ii,:)) - min(xx(ii,:)));
% end
for ii = 1:size(xx,2) % normalize along the spectra
y(:,ii) = (xx(:,ii) - min(xx(:,ii)))./(max(xx(:,ii)) - min(xx(:,ii)));
end
yy = reshape(y,size(x,2),size(x,2),278);
% for pseudocolor visulaization use bands 1,10 and 30 after band removal.