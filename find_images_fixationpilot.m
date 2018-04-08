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
idmap = [catList.id];
low_cat = {catList.name};
%% Reshape into long form ish and remove annotations <5 %
% we don't want images that have annotations that are just too small
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
    keep(i) = isempty(anns) || ~any(anndata(:,2)<0.01); % keep indexes that have no annotations greater than 20% of image
    disppercent(i/N);
end
disppercent(inf);

oidx = find(keep); % save the original indexes so we can quickly pull the full info from the coco dataset
imdata = imdata(logical(keep),:);
imanns = imanns(logical(keep));

%% For every image, load annotations and create an img * category matrix.

imcats = zeros(size(imdata,1),80);

for ii = 1:size(imdata,1)
    anns = coco.loadAnns(coco.getAnnIds('imgIds',imdata(ii,1)));
    for ai = 1:length(anns)
        imcats(ii,idmap==anns(ai).category_id)=1;
    end
end


%% Go through images and find images so that you have an even # of every category

pickIds = zeros(1080,1); % 
pickAnns = zeros(1080,80);


imcount = 1;
availidxs = ones(size(imdata,1),1);

for ii = 1:1080
    % pick which category we want to add
    cimgcat = sum(pickAnns(1:(ii-1),:)); % compute # in all categories
    cimgcat = cimgcat + [inf zeros(1,79)];
    [~,minidx] = min(cimgcat); % get the idx of the minimum category
    % get all images which include that category
    catidxs = imcats(:,minidx)==1;
    % pick one at random then remove
    avail = find(availidxs.*catidxs); % get all available indexes for this category
    avail = avail(randperm(length(avail))); % randomze
    pick = avail(1); % pick one
    availidxs(pick)=0;
    % save pick

    pickIds(ii) = imdata(pick,1);
    pickAnns(ii,:) = imcats(pick,:);
end

% pickIds = pickIds(1:(count-1));
% pickAnns = pickAnns(1:(count-1),:);

sum(pickAnns)

%% Save data picked for semantic processing
semdata = struct;
semdata.ids = pickIds;
semdata.anns = pickAnns;

save(fullfile('~/proj/fbsear/data/','semantic.mat'),'semdata');

%% Covariance matrix

imagesc(cov(ntdata))

%% Go through chosen images and re-pick images that suck
figure(1);
ii = randi(size(pickIds,1));
img = coco.loadImgs(pickIds(ii));
I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[]);
    

%%
newterms = {'outdoor','indoor','water','text','snow'};
ntdata = zeros(1080,length(newterms));
done = zeros(1080,length(newterms));

%% Other terms to add

figure(1);
for ni = 1:length(newterms)
    disp('*******************************************************');
    disp('*************** NEW TERM STARTING *********************');
    disp('*******************************************************');
    disppercent(-1/1080);
    for ii = 1:1080
        if ~done(ii,ni)
            img = coco.loadImgs(pickIds(ii));
            I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
            imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[]);
            pause(.25);
            key = waitForKeypress;

            if strcmp(key,'y')
                disp(sprintf('%s',newterms{ni}));
                ntdata(ii,ni) = 1;
                done(ii,ni) = 1;
            elseif strcmp(key,'n')
                disp(sprintf('not %s',newterms{ni}));
                ntdata(ii,ni) = 0;
                done(ii,ni) = 1;
            else
                disp('ignored');
            end
        end
        disppercent(ii/1080);
    end
    disppercent(inf);
end