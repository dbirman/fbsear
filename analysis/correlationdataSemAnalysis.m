%% Semantic analysis for fbsear_pilot
% We will build a semantic category matrix and use that to estimate the GLM
% beta weights 

%% Get the category label groups
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

% We'll use these later to label the categories when we add them to the
% stimulus files
idmap = [catList.id];
cats = {catList.name};

%% Build the semantic category matrix for all categories

% Load the image info
load(fullfile('~/proj/fbsear/data/info_m.mat'));

categories = {'person','car','personcar','null'};

cat_info = info;
cat_info.cats = struct;

for ci = 1:length(categories)
    % for each category, pull the images
    adat = info.anns.(categories{ci});
    curcats = zeros(length(adat),80);
    
    for ii = 1:length(adat)
        % for each image grab its annotations
        canns = adat{ii};
        
        % set the category id positions
        for ai = 1:length(canns)
            curcats(ii,idmap==canns(ai).category_id) = curcats(ii,idmap==canns(ai).category_id) + 1;
        end
    end
    
    cat_info.cats.(categories{ci}) = curcats;
end

%% Now for the actual images used in the experiment
% load each stimfile and pull the trial/cat order and then construct a
% design matrix based on these. We'll use the design matrices from all the
% runs combined to estimate which categories we should keep and which we
% should dump. Then we'll add the relevant categories back in as calculated
% variables using "addCalculatedVar". To avoid contaminating files we'll
% backup the entire /Etc folder before we start modifying anything.
cfolder = 's030020180124';
mrQuit
cd(fullfile('~/data/fbsear/',cfolder));
folders = dir(pwd);

%% Set files
stimfolder = fullfile(pwd,'Etc');
backupfolder = fullfile(stimfolder,'Backup');

%% backup
if ~isdir(backupfolder), mkdir(backupfolder); end
system(sprintf('cp %s/* %s/',stimfolder,backupfolder));

%% Load stimfiles and get trial/cat order

stimfiles = dir(fullfile(stimfolder,'*stim*.mat'));

allAnnShown = zeros(0,60,80);

keep = [];
for si = 1:length(stimfiles)
    if isempty(strfind(stimfiles(si).name,'original')) % ignore original files
        keep(end+1) = si;
        % load this file
        load(fullfile(stimfolder,stimfiles(si).name));
        run = stimulus.curRun;
        
        % get the index in the info file for each image that was shown
        idxShown = zeros(1,60);
        catShown = zeros(1,60);
        for ii = 1:60
            trial = run.trialOrder(ii);
            cat = run.catOrder(ii);
            idxShown(ii) = stimulus.runs.imageIndexes(cat,run.repeat,trial);
            catShown(ii) = cat;
        end
        % get the actual image numbers and their connected category
        % information
        imgShown = zeros(1,60);
        annShown = zeros(60,80);
        for ii = 1:60
            imgs = info.(categories{catShown(ii)});
            imgShown(ii) = imgs(idxShown(ii));
            anns = cat_info.cats.(categories{catShown(ii)});
            annShown(ii,:) = anns(idxShown(ii),:);
        end
        allAnnShown(end+1,:,:) = annShown;
    end
end
stimfiles = stimfiles(keep);

%% Collapse to only unique

uniqueAnnShown = allAnnShown(1:5:end,:,:);
uniqueAnnShown = reshape(uniqueAnnShown,120,80);
uniqueAnnShown = uniqueAnnShown>0;
% count instances of each condition
catInstances = sum(uniqueAnnShown);
% display categories that we have at least 8 showings
disp(cats(catInstances>7));
keep_idxs = catInstances>7;
keep_nidxs = find(keep_idxs);
keep_cat = cats(keep_idxs);
keep_ann = allAnnShown(:,:,keep_idxs);

% replace category names by single words
for ci = 1:length(keep_cat)
    if ~isempty(strfind(keep_cat{ci},' '))
        keep_cat{ci} = keep_cat{ci}(1:(strfind(keep_cat{ci},' ')-1));
    end
end
%% Build design matrix and addCalculatedVars
% we'll ignore the fact that multiple things often show up in the same
% scene for simplicity 

cd(stimfolder);

parfor si = 1:length(stimfiles)
    design = squeeze(keep_ann(si,:,:))>0;
    design2 = NaN(size(design));
    design2(design==1) = 1;

    for di = 1:size(design,2)
        addCalculatedVar(keep_cat{di},design2(:,di),stimfiles(si).name,'backup=0','force=1');
    end
    
end
disp('done');

%% generate string for analysis (put spaces between all the words

disp(strcat('{',sprintf('{''%s=1''},',keep_cat{:}),'}'));






%% RUN GLM ANALYSIS IN MLR




%% Pull GLM analysis

cd(fullfile('~/data/fbsear/',cfolder));
folders = dir(pwd);
skip = 1;
for fi = 1:length(folders)
    if ~isempty(strfind(folders(fi).name,'Concatenation')), skip = 0; end
end
if skip
    disp(sprintf('Data folder %s has not been prepared for analysis',cfolder));
    return
end

view = newView();
view = viewSet(view,'curGroup','Concatenation');


%% Load analysis
view = loadAnalysis(view,'glmAnalStats/GLM'); % check analysis name!

%% Pull ehdr for each run and r2
clear ehdr
for group = 1:5
    ehdr(group,:,:,:,:) = view.analyses{1}.d{group}.ehdr;
    r2(group,:,:,:) = view.analyses{1}.overlays(1).data{group};
end

%% For each voxel compute correlations between groups -- generate overlay for the differences
% Overlay 1: Person effect, corr(look p, con p) - corr(look p, con c)
% Overlay 2: Car effect, corr(look c, con c) - corr(look c, con p)
effect_p = zeros(size(ehdr,2),size(ehdr,3),size(ehdr,4));
effect_c = zeros(size(effect_p));
disppercent(-1/size(ehdr,2));
for xi = 1:size(ehdr,2)
    for yi = 1:size(ehdr,3)
        for zi = 1:size(ehdr,4)
            % compute matrix of correlations
            corrs = corrcoef(squeeze(ehdr(:,xi,yi,zi,:))');
            % compute effects
            effect_p(xi,yi,zi) = corrs(2,4)-corrs(2,5);
            effect_c(xi,yi,zi) = corrs(3,5)-corrs(3,4);
        end
    end
    disppercent(xi/size(ehdr,2));
end
disppercent(inf);


%% Load MLR and display these overlays

mrQuit;
mrLoadRet;

% viewSet(getMLRView,'curGroup','Concatenatation');
%%
mrDispOverlay(effect_p,2,view.curGroup,getMLRView,'overlayName=effectPeople');
mrDispOverlay(effect_c,2,view.curGroup,getMLRView,'overlayName=effectCar');

%% For each group, pull the beta weights for the >95%

ar2 = r2(:);
q95 = quantile(ar2,.95);
idxs = reshape(ar2>q95,size(r2));
ghdr = ehdr;
ghdr(repmat(idxs,1,1,1,1,size(ehdr,5))==1) = NaN;
% nanmean

for gi = 1:5
    data = squeeze(ghdr(gi,:,:,:,:));
    data = reshape(data,size(data,1)*size(data,2)*size(data,3),size(data,4));
    out(gi,:) = nanmean(data);
    outs(gi,:) = nanstd(data);
end
plot(1:11,out','o');
errbar(1:11,out',outs');
set(gca,'XTick',1:11,'XTickLabel',keep_cat);
legend(groups);

%% Compare betas
groups = {'fixate','look people','look cars','contrast people','contrast cars'};
plot(squeeze(ehdr(:,50,50,50,:))');
set(gca,'XTick',1:11,'XTickLabel',keep_cat);
legend(groups);