%% Load file
load(fullfile('~/proj/fbsear/data/info_m.mat'));


categories = {'person','car','personcar','null'};

%% Save original images
for ci = 1:length(categories)
    imgs = info.imgs.(categories{ci});
    idir = fullfile('~/proj/fbsear/images/imgs/',categories{ci});
    if ~isdir(idir), mkdir(idir); end
    for ii = 1:length(imgs)
        cimg = imgs{ii};
        imwrite(cimg,fullfile(idir,sprintf('image%i.png',ii)),'PNG');
    end
end

clear imgs cimg idir

%% Save mask images
for ci = 1:length(categories)
    apimgs = info.mimgs.attend_person.(categories{ci});
    acimgs = info.mimgs.attend_car.(categories{ci});
    
    pdir = fullfile('~/proj/fbsear/images/mimgs/attend_person',categories{ci});
    if ~isdir(pdir), mkdir(pdir); end
    cdir = fullfile('~/proj/fbsear/images/mimgs/attend_car',categories{ci});
    if ~isdir(cdir), mkdir(cdir); end
    for ii = 1:length(apimgs)
        pimg = apimgs{ii};
        imwrite(pimg,fullfile(pdir,sprintf('image%i.png',ii)),'PNG');
        cimg = acimgs{ii};
        imwrite(cimg,fullfile(cdir,sprintf('image%i.png',ii)),'PNG');
    end
end
