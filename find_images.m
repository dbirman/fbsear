%% Load all the (test?) images
addpath(genpath('~/proj/fbsear'));

annTypes = { 'instances', 'captions', 'person_keypoints' };
dataType='val2017'; annType=annTypes{1}; % specify dataType/annType
annFile=fullfile('~/proj/fbsear/',sprintf('annotations/%s_%s.json',annType,dataType));
coco=CocoApi(annFile);

if( ~strcmp(annType,'captions') )
  catList = coco.loadCats(coco.getCatIds());
  nms={catList.name}; fprintf('COCO categories: ');
  fprintf('%s, ',nms{:}); fprintf('\n');
  nms=unique({catList.supercategory}); fprintf('COCO supercategories: ');
  fprintf('%s, ',nms{:}); fprintf('\n');
end

high_cat = unique({catList.supercategory});
low_cat = {catList.name};
%% Reshape into long form ish and remove annotations > 20%
% create long form of original images and areas

% ID | AREA
imdata = zeros(length(coco.data.images),2);
% {[ANN ID | % AREA ]}
imanns = cell(1,length(coco.data.images));
% keep indexes
keep = zeros(1,length(coco.data.images));
for i = 1:length(coco.data.images)
    id = coco.data.images(i).id;
    area = coco.data.images(i).height*coco.data.images(i).width;
    imdata(i,:) = [id area];
    anns = coco.loadAnns(coco.getAnnIds('imgIds',id));
    anndata = zeros(length(anns),2);
    for ai = 1:length(anns)
        anndata(ai,:) = [anns(ai).id anns(ai).area/area];
    end
    imanns{i} = anndata;
    keep(i) = ~any(anndata(:,2)>0.33); % keep indexes that have no annotations greater than 20% of image
end

oidx = find(keep); % save the original indexes so we can quickly pull the full info from the coco dataset
imdata = imdata(logical(keep),:);
imanns = imanns(logical(keep));

%% Find images in the eight combination categories

base = {{'car','dog','person'}
        {'car','person'}
        {'dog','person'}
        {'car','dog'}
        {'person'}
        {'car'}
        {'dog'}
        {}};
allCats = {''};

imdata_ = [imdata zeros(size(imdata,1),8)];
% for each category combination load all the relevant images then add to
% the imdata array to tell which of the eight groups that image falls under
for bi = 1:length(base)
    catIds = coco.inds.catIds;
    baseidx{bi} = cellfun(@(x) find(cellfun(@(y) strcmp(y,x),low_cat),1),base{bi},'UniformOutput',false);
    imgIds = coco.getImgIds('catIds',catIds([baseidx{bi}{:}]));
    disp(sprintf('For base %s I found %i images',[base{bi}{:}],length(imgIds)));
    for i=1:size(imdata_,1)
        imdata_(i,2+bi) = any(imgIds==imdata_(i,1));
    end
    if bi>1
        imdata_(:,2+bi) = imdata_(:,2+bi) .* ~(sum(imdata_(:,3:(1+bi)),2)>0);
    end
    disp(sprintf('Removing images from previous categories left me with %i',sum(imdata_(:,2+bi))));
end

%% Display four random images from a category
category = 3; % category choice

imgIds = imdata_(imdata_(:,2+category)==1,1);

disp(sprintf('For base %s I found %i images',[base{category}{:}],length(imgIds)));

h = figure(1);
% 4x4 array of random images
n = length(imgIds);
for x = 1:2
    for y = 1:2
        idx = (x-1)*2+y;
        subplot(2,2,idx);
        img = coco.loadImgs(imgIds(randi(n)));
        I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
        imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[])
    end
end