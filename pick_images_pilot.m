%% Load images that have outdoor information:
% Goal is to pick 30 images that have people, cars, and people/cars
% Then pick 90 images that have neither people or cars (but some of the
% other categories like dogs, or other random shit).
%
% To do this we'll pick random images and ask whether these are okay or
% not. Once images have all been filled out into the categories then the
% program will start asking to replace images (i.e. do you prefer image X
% or image Y, where X is the current image and Y is a random image taken
% from the outdoor set).
%
% The images will then be passed to the function parseImages, which will
% save the information about each category as well

% This script saves a struct *info* which contains information about the
% categories (person, car, personcar, and null) which each have 30/30/30/90
% images, respectively. 

%% Setup

categories = {'person','car','personcar','null'};
cat_base = [5,6,2,8];
cat_num = [60,60,60,60];

%% Load all images

load(fullfile('~/proj/fbsear/imdata_outdoor.mat'));
% imdata_ is now loaded

%% Load the info file
finfo = fullfile('~/proj/fbsear/data/info.mat');
if isfile(finfo)
    load(finfo);
else
    info = struct;
    for ci = 1:length(categories)
        info.(categories{ci}) = [];
    end
end

%% Get category numbers

for ci = 1:length(categories)
    c_cat(ci) = length(info.(categories{ci}));
end 

%% Go through each category and insert or replace images as necessary
figure;
for ci = 1:length(categories)
    disp(sprintf('Curating category: %s',categories{ci}));
    if length(info.(categories{ci}))<cat_num(ci)
        disp('Adding images');
        subplot(111);
        % we are inserting images
        while length(info.(categories{ci}))<cat_num(ci)
            % load a random image from this category
            b_cat = cat_base(ci);
            b_cat_idxs = find(imdata_(:,2+b_cat)==1);
            % avoid repeats
            cdat = imdata_(b_cat_idxs(randi(length(b_cat_idxs),1)),:);
            canns = coco.loadAnns(coco.getAnnIds('imgIds',cdat(1)));
            while any(cdat(1)==info.(categories{ci}))% || any([canns.iscrowd])
                cdat = imdata_(b_cat_idxs(randi(length(b_cat_idxs),1)),:);
            end
            % load image
            img = coco.loadImgs(cdat(1));
            I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
            imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[]);
            pause(.25);
            key = waitForKeypress;

            if strcmp(key,'y')
                disp('accepted');
                info.(categories{ci}) = [info.(categories{ci}) cdat(1)];
            else
                disp('ignored');
            end
        end
    else
        disp('Replacing images');
        % load all category info
        b_cat = cat_base(ci);
        b_cat_idxs = find(imdata_(:,2+b_cat)==1);
        
        % we are replacing images
        imgs = info.(categories{ci});
        for ii = 1:length(imgs)
            cont = true;
            brk = false;
            while cont
                img = imgs(ii);
                % get data
                cdat = imdata_(imdata_(:,1)==img,:);
                %
                % avoid repeats
                ndat = imdata_(b_cat_idxs(randi(length(b_cat_idxs),1)),:);
                while any(ndat(1)==info.(categories{ci}))
                    ndat = imdata_(b_cat_idxs(randi(length(b_cat_idxs),1)),:);
                end
                

                % load images
                cimg = coco.loadImgs(cdat(1));
                canns = coco.loadAnns(coco.getAnnIds('imgIds',cdat(1)));
                if ~any([canns.iscrowd])
                    nimg = coco.loadImgs(ndat(1));
                    cI = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,cimg.file_name)));
                    nI = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,nimg.file_name)));
                    subplot(121);
                    imagesc(cI); axis('image'); set(gca,'XTick',[],'YTick',[]);
                    title('current');
                    % check that no is crowd -- otherwise auto replace
                    if any([canns.iscrowd])
                        title('IS CROWD!!!');
                    end
                    subplot(122);
                    imagesc(nI); axis('image'); set(gca,'XTick',[],'YTick',[]);
                    title('replacement');
                    pause(.25);
                    key = waitForKeypress;

                    if strcmp(key,'y')
                        disp('replaced');
                        dat = info.(categories{ci});
                        dat(ii) = ndat(1);
                        info.(categories{ci}) = dat;
                    elseif strcmp(key,'q')
                        cont = false;
                    elseif strcmp(key,'b')
                        brk = true; break
                    else
                        disp('continuing');
                    end
                    imgs = info.(categories{ci});
                end
            end
            if brk
                break;
            end
        end
    end
end

%% Save info
save(finfo,'info');

%% Parse

parseImages;

%% View images
figure;
for ci = 1:length(categories)
    imgs = info.(categories{ci});
    for ii = 1:length(imgs)
        img = coco.loadImgs(imgs(ii));
        I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
        imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[]);
        pause(.25);
        key = waitForKeypress;
    end
end