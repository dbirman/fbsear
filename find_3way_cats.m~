%% Find 3-category images
% We want to find for combination of three categories (there are 80
% low-level categories: so 82,160 combinations) where the segmentations are
% at least 5% and less than 40% of the available space.
addpath(genpath('~/proj/fbsear'));

%% initialize COCO api (please specify dataType/annType below)
annTypes = { 'instances', 'captions', 'person_keypoints' };
dataType='train2017'; annType=annTypes{1}; % specify dataType/annType
annFile=fullfile('~/proj/fbsear/',sprintf('annotations/%s_%s.json',annType,dataType));
coco=CocoApi(annFile);

%% display COCO categories and supercategories
if( ~strcmp(annType,'captions') )
  catList = coco.loadCats(coco.getCatIds());
  nms={catList.name}; fprintf('COCO categories: ');
  fprintf('%s, ',nms{:}); fprintf('\n');
  nms=unique({catList.supercategory}); fprintf('COCO supercategories: ');
  fprintf('%s, ',nms{:}); fprintf('\n');
end

%% Split cats
high_cat = unique({catList.supercategory});
low_cat = {catList.name};

% %% get all images containing given categories, select one at random
% imgIds = [];
% while isempty(imgIds)
%     catIds = coco.getCatIds('catNms',low_cat(randi(length(low_cat),1,3)));
%     imgIds = coco.getImgIds('catIds',catIds);
% end
% 
% imgId = imgIds(randi(length(imgIds)));
% 
% %% load and display image
% img = coco.loadImgs(imgId);
% I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
% figure(1); imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[])
% 
% %% load and display annotations
% annIds = coco.getAnnIds('imgIds',imgId,'catIds',catIds,'iscrowd',[]);
% anns = coco.loadAnns(annIds); coco.showAnns(anns);

%% Get all category combinations
aCatIds = [catList.id];
cats = {}; imgs = {}; imgN = [];

allCombs = combinator(length(low_cat),3,'c');
N = size(allCombs,1);
disppercent(-1/N);
for ai = 1:N
    catIds = allCombs(ai,:);
    catIds = aCatIds(catIds);
    
    imgIds = coco.getImgIds('catIds',catIds);

    if length(imgIds)>100
        cats{end+1} = catIds;
        imgs{end+1} = imgIds;
        imgN(end+1) = length(imgIds);
    end
        
    disppercent(ai/N);
end
disppercent(inf);

%% Save data
data.cats = cats;
data.imgs = imgs;
data.imgN = imgN;

save('~/proj/fbsear/data.mat','data');

%% Load
load('~/proj/fbsear/data.mat');

%% Get category names
catNames = {};
for ai = 1:N
    catNames{end+1} = {catList(arrayfun(@(x) find(x==aCatIds,1),cats{ai})).name};
end

%% Go through 10 example images from each category
h = figure(1);
cN = length(cats);
reorder = randperm(cN);
for aii = 1:cN
    ai = reorder(aii);
    disp(strcat(catNames{ai}));
    % 4x4 array of random images
    n = length(imgs{ai});
    for x = 1:4
        for y = 1:4
            idx = (x-1)*4+y;
            subplot(4,4,idx);
            img = coco.loadImgs(imgs{ai}(randi(n)));
            I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
            imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[])

        end
    end
    pause(10);
end

%% Run through all images removing any that have more than one of the object and a segmentation of <3% or >33%

h = figure(1);
cN = length(cats);
fdata = {};

disppercent(-1/cN);
for ai = 1:cN
%     disp(sprintf('Category: %s, %s, %s',catNames{ai}{1},catNames{ai}{2},catNames{ai}{3}));
    % load all images
    cImgs = coco.loadImgs(imgs{ai});
    % get annotations separately for each category
    for ii = 1:length(cImgs)
        for ci = 1:3
            cAnns{ii,ci} = coco.loadAnns(coco.getAnnIds('imgIds',imgs{ai}(ii),'catIds',cats{ai}(ci)));
        end
    end
    % remove all images that had multiple annotations
    kIdxs = [];
    for ii = 1:length(cImgs)
        if all(cellfun(@(x) length(x)==1,cAnns(ii,:)))
            % keep
            kIdxs(end+1) = ii;
        end
    end
%     disp(sprintf('Keeping %2.0f%% for a total of %i',length(kIdxs)/length(cImgs)*100,length(kIdxs)));
    ncImgs = cImgs(kIdxs);
    ncAnns = cAnns(kIdxs,:);
    % run through and check segmentation sizes
    kIdxs = [];
    for ii = 1:length(ncImgs)
        totalArea = ncImgs(ii).width*ncImgs(ii).height;
        keep = true;
        for ci = 1:3
            % check each annotation
            percArea = ncAnns{ii,ci}.area/totalArea;
            if percArea<.03 || percArea>.33
                keep = false; break
            end
        end
        if keep
            kIdxs(end+1) = ii;
        end
    end
%     disp(sprintf('Keeping %2.0f%% for a total of %i',length(kIdxs)/length(ncImgs)*100,length(kIdxs)));
    fcImgs = ncImgs(kIdxs);
    fcAnns = ncAnns(kIdxs,:);
    % Display the first 16 remaining images
    % 4x4 array of random images
% %     n = length(fcImgs);
% %     clf
% %     for x = 1:4
% %         for y = 1:4
% %             idx = (x-1)*4+y;
% %             if idx<n
% %                 subplot(4,4,idx);
% %                 img = coco.loadImgs(fcImgs(idx).id);
% %                 I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
% %                 imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[])
% %             end
% %         end
% %     end
% %     pause(10);
    if length(fcImgs)>10
        dat.cats = catNames{ai};
        dat.imgs = fcImgs;
        dat.anns = fcAnns;
        fdata{end+1} = dat;
    end
    disppercent(ai/cN);
end
disppercent(inf);

%% Save data

save('~/proj/fbsear/fdata.mat','fdata');

%% Load
load('~/proj/fbsear/data.mat');

%% Get image lengths
lens = cellfun(@(x) length(x.imgs),fdata);

%% Find the category overlap between sets
catCount = zeros(size(catList));
for fi = 1:length(fdata)
    % this is one set of 3 objects in at least ten images, add each
    % category to its position
    for ci = 1:3
        cCat = fdata{fi}.cats{ci};
        n = length(fdata{fi}.imgs);
        for idx = 1:80
            if strcmp(cCat,catList(idx).name)
                catCount(idx) = catCount(idx) + n;
            end
        end
    end
end

% Names of the most common categories
[s,i] = sort(catCount,'descend');
i = i(s>0);
s = s(s>0);
fCatNames = {};
for si = 1:length(s)
    fCatNames{end+1} = catList(i(si)).name;
end