%% Save images for Keras (save each individual image as a .mat file separately)

% first save original images
folder = '~/proj/fbsear/data/orig/';
categories = {'person','car','personcar','null'};


for ci = 1:length(categories)
    
    if ~isdir(fullfile(folder,categories{ci}))
        mkdir(fullfile(folder,categories{ci}));
    end
    
    cimgs = info.imgs.(categories{ci});
    
    for ii = 1:length(cimgs)
        cimg = cimgs{ii};
        
        save(fullfile(folder,categories{ci},sprintf('img%02.0f.mat',ii)),'cimg','-v7');
    end
end

% now save attend_person

% now save attend_car