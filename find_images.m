%% Load all the (test?) images
addpath(genpath('~/proj/fbsear'));

annTypes = { 'instances', 'captions', 'person_keypoints' };
dataType='train2017'; annType=annTypes{1}; % specify dataType/annType
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
N = length(coco.data.images);
disppercent(-1/N);
for i = 1:N
    id = coco.data.images(i).id;
    area = coco.data.images(i).height*coco.data.images(i).width;
    imdata(i,:) = [id area];
    anns = coco.loadAnns(coco.getAnnIds('imgIds',id));
    anndata = zeros(length(anns),2);
    for ai = 1:length(anns)
        anndata(ai,:) = [anns(ai).id anns(ai).area/area];
    end
    imanns{i} = anndata;
    keep(i) = isempty(anns) || ~any(anndata(:,2)>0.33); % keep indexes that have no annotations greater than 20% of image
    disppercent(i/N);
end
disppercent(inf);

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

%% add images that have NO categories in them
base{end+1} = {'nococo'};
nococo_imgs;

search_idxs = imdata_(:,end)==1;
s = size(imdata_,2);
s1 = size(imdata_,2)+1;
for ni = 1:length(nococo)
    cid = nococo(ni);
    imdata_idx = find(imdata_(:,1)==cid,1);
    imdata_(imdata_idx,s1) = 1;
    imdata_(imdata_idx,s) = 0;
end

%% Save imdata

save('~/proj/fbsear/imdata.mat','imdata_');









%% OLD CODE OLD CODE OLD CODE OLD CODE

% NOTE: You can't take out images, most of the images are random sizes but
% have one dimension either 640 or 500. 

%% Check how many images are left if we remove images that are vertical (keep only 640x480)

awidth = [coco.data.images.width];
aheight = [coco.data.images.height];

% figure; plot(awidth,aheight,'*');
% figure; hist([awidth aheight]);

cocoIds = [coco.data.images.id];

keep = zeros(1,size(imdata_,1));
N = length(cocoIds);
disppercent(-1/N);
for ii = 1:size(imdata_,1)
    cid = imdata_(ii,1);
    idx = find(cocoIds==cid,1);
    keep(ii) = awidth(idx)==640 || awidth(idx)==500;
    disppercent(ii/N);
end
disppercent(inf);

imdata_sq = imdata_(logical(keep),:);

disp(sum(imdata_sq(:,3:end)));

%% Display four random images from a category
category = 9; % category choice

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

%% Squareify?

category = 8; % category choice

imgIds = imdata_(imdata_(:,2+category)==1,1);

h = figure(1);
% 4x4 array of random images
n = length(imgIds);
for x = 1:2
    for y = 1:2
        idx = (x-1)*2+y;
        subplot(2,2,idx);
        img = coco.loadImgs(imgIds(randi(n)));
        I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
        
%         Isq = zeros(640,640,3);
%         pixels = I(:);
%         if size(I,1)<640 && size(I,2)<640
%             Isq = I;
%         elseif size(I,1)<640
%             dif = 640-size(I,1);
%             ld = floor(dif/2); rd = ceil(dif/2);
%             Isq = cat(1,reshape(pixels(randi(length(pixels),1,ld*640*3)),ld,640,3),I,zeros(rd,640,3));
%         else size(I,2)<640
%             Isq = I;
% %             dif = 640-size(I,2);
% %             ld = floor(dif/2); rd = ceil(dif/2);
% %             Isq = cat(2,reshape(pixels(randi(length(pixels),1,ld*640*3)),ld,640,3),I,zeros(rd,640,3));
%         end
        clear Isq
        if size(I,1)<640 || size(I,2)<640
            dif = 640-size(I,1);
            ld = floor(dif/2); rd = ceil(dif/2);
            dif = 640-size(I,2);
            ud = floor(dif/2); dd = ceil(dif/2);
            Isq = padarray(I,[ld ud],'replicate','pre');
            Isq = padarray(Isq,[rd dd],'replicate','post');
        else
            Isq = I;
        end
        imagesc(Isq); axis('image'); set(gca,'XTick',[],'YTick',[])
    end
end

%% For each image 