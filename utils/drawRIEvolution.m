function drawRIEvolution(results,rps,init,rr,cc,cmap)
%%Draws the evolution of the RoIs as the RPS progresses.
%%results should contain the recosntructed RoIs at various resolutions.
%%It should be like results.ri.xr geerated from the finalCode file.IT
%%shouls also contain the low resolution recon in results.lowResRecon.
%%rps is a structure and should contain ROI info in rps.seg.ROI. It is
%%same as the rps.seg in the finalCode file. 
%%init should contain init.m and init.n, that define the image size.
%%rr are the number of rows in the subplot and cc are the number of
%%columns.
%%The first two subplots contains the low resolution reconstrucion of the
%%entire scene and the coarse recosntruction of each region.
%%cmap is required for color the borders of the RoIs. Same as cmap in the
%%finalCode file.
%%Sometimes if you are loading from a saved .mat file, you may need to
%%generate the rps.seg.opR variables using the following command
%%rps.seg.opR{1,kk} = opRestriction(init.m*init.n,find(~rps.seg.maskC{1,kk}));

coarseRecon = results.lowResRecon;

for kk = 1:rps.seg.numRoIs
coarseRecon = rps.seg.maskC{1,kk}.*coarseRecon + reshape(rps.seg.opR{1,kk}'*results.ri.xr{kk,kk}(:),init.m,init.n);
end
fullImage = coarseRecon;
for kk = 1:rps.seg.numRoIs    
coarseRecon = insertShape(coarseRecon,'rectangle',[rps.seg.cstrt(1,kk),rps.seg.rstrt(1,kk),...
    rps.seg.roim(1,kk),rps.seg.roin(1,kk)],'Color', cmap(kk,:), 'LineWidth', 3);
end


figure(); subplot(rr,cc,1); imshow(results.lowResRecon); 
subplot(rr,cc,2); imshow(coarseRecon)
plot_count = 3;
for kk = rps.seg.numRoIs+1:size(results.ri.xr,2)
    t = find(~cellfun(@isempty,results.ri.xr(:,kk)));
    fullImage = rps.seg.maskC{1,t}.*fullImage + reshape(rps.seg.opR{1,t}'*results.ri.xr{t,kk}(:),init.m,init.n);
    im = insertShape(fullImage,'rectangle',[rps.seg.cstrt(1,t),rps.seg.rstrt(1,t),...
        rps.seg.roim(1,t),rps.seg.roin(1,t)],'Color',cmap(t,:),'LineWidth', 3);
    subplot(rr,cc,plot_count)
    imshow(im)
    plot_count = plot_count + 1;
    
end
