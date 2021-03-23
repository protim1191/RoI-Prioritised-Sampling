function [x1,init] = loadDataset_github(dataset,imSize)
%load data
% dataset  = cuprite; ImSize = 256(Vis), 3(IR), 512 (512x512 Vis)
%          = gulf; imSize = 256;
%          = pavia; imSize = 256;
%
init.dataset = strcat(dataset,'_',num2str(imSize));
switch(dataset)
    case 'cuprite'
        switch(imSize)
            case 512
                [~,x] = preprocessHyperspectralImageCuprite512;
                init.m = 512; init.n = 512;
                x1 = sum(x(:,1:41),2)/41;
            case 256
                [~,x] = preprocessHyperspectralImageCuprite1;
                init.m = 256; init.n = 256;
                x1 = sum(x(:,1:41),2)/41;
            case 3
                [~,x] = preprocessHyperspectralImageCuprite1;
                init.m = 256; init.n = 256;
                x1 = sum(x(:,61:160),2)/100;
            otherwise
                str = strcat('Size ', num2str(imSize),' for ', dataset,' unavailable');
                error(str)
        end
    case 'gulf'
        switch (imSize)
            case 256
                [~,x] = preprocessHyperspectralImageMexicoGulf1;
                init.m  = 256; init.n = 256;
                x1 = sum(x(:,1:54),2)/54;
            otherwise
                str = strcat('Size ', num2str(imSize),' for ', dataset,' unavailable');
                error(str)
        end
    case 'pavia'
        switch(imSize)
            case 256
                x1 = load('pavia256');
                x1 = x1.pavia256(:);
                init.m = 256; init.n = 256;
            otherwise
                str = strcat('Size ', num2str(imSize),' for ', dataset,' unavailable');
                error(str)
        end
    otherwise
        str = strcat(dataset,' unavailable');
        error(str)
end