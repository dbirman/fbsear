% Take images in info.cat and generate two sets of images saved in
% ~/proj/fbsear/data/cat and ~/proj/fbsear/data/con_cat. The second set are
% contrast enhanced to make the "search" category more visible

load(fullfile('~/proj/fbsear/data/info.mat'));

%% Use COCO to load all images
% figure(1);

categories = {'person','car','personcar','null'};
for ci = 1:2
    cidx_cell = cellfun(@(x) strcmp(x,categories{ci}),low_cat,'UniformOutput',false);
    cat_orig(ci) = find([cidx_cell{:}]);
end
cat_orig(3) = -1;
cat_orig(4) = -1;

cat_base = [5,6,2,8];
cat_num = [30,30,30,90];
info.imgs = struct;
info.anns = struct;
for ci = 1:length(categories)
    
    imgs = info.(categories{ci});
    info.imgs.(categories{ci}) = {};
    info.anns.(categories{ci}) = {};
    
    for ii = 1:length(imgs)
        img = coco.loadImgs(imgs(ii));
        imgname = sprintf('img_%i',imgs(ii));
        I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
        info.imgs.(categories{ci}){end+1} = I;
        
        anns = coco.loadAnns(coco.getAnnIds('imgIds',imgs(ii)));
        info.anns.(categories{ci}){end+1} = anns;
        
        % test code: display image and annotations
%         figure(1); imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[])
%         coco.showAnns(anns);
%         
% %         pause(1);
    end
end

%% Go through images and boost the contrast for people/cars: BOUNDING BOX MODE

figure(1);

% percent mask
mask_perc = 0.6;
mode = 'seg';

% seg blur kernel (5x5)
sigma = 15;
r = sigma * 3; % range
sx = -r:r;
sy = -r:r;
[SX,SY] = meshgrid(sx,sy);
skernel = normpdf(hypot(SX,SY),0,5);

info.mimgs = struct;
% This is complex, we have to get the annotation, build a blurry polygon
% mask and then use that to increase or decrease the contrast throughout
% the image.

% we only have to do this for the first two categories (since they are the
% only ones with cars and people)


for ci = 1:4
    % category of imgs to load
    imgs = info.imgs.(categories{ci});
    anns = info.anns.(categories{ci});
    
    
    for aci = 1:2
        % boost category, which category to pay attention to the
        % annotations for
        accat = cat_orig(aci);
    
        mimgs = cell(size(imgs));

        % load each image, then load its annotations
        for ii = 1:length(imgs)
            cimg = imgs{ii};
            canns = anns{ii};

            % build the meshgrid for this image
            x = size(cimg,2);
            y = size(cimg,1);
            [X,Y] = meshgrid(1:x,1:y);

            ann_cats = [canns.category_id];
            relev_idxs = find(ann_cats==accat);
            
            if ~isempty(relev_idxs)

                % get the number of gaussians that will be added
                mask = zeros(y,x,length(relev_idxs));
                crowd = zeros(1,length(relev_idxs));
                % for each annotation that matches the current category id (in
                % cat_orig) stack a gaussian at that location, then average
                for mi = 1:length(relev_idxs)
                    cidx = relev_idxs(mi);
                    cann = canns(cidx);

                    if ~cann.iscrowd
                        if strcmp(mode,'bbox')
                            % get bbox center
                            cx = cann.bbox(1) + cann.bbox(3)/2;
                            cy = cann.bbox(2) + cann.bbox(4)/2;
                            bb_diag = hypot(cann.bbox(3),cann.bbox(4))/2;

                            % compute mask
                            cmask = normpdf(hypot(X-cx,Y-cy),0,bb_diag);
                            % normalize to 0->1
                            cmask = cmask ./ max(cmask(:));
                        elseif strcmp(mode,'seg')
                            % get bw mask
                            s = cann.segmentation{1};
                            xx = s(1:2:end);
                            yy = s(2:2:end);
                            cmask = double(poly2mask(xx,yy,y,x));
                            cmask = imfilter(cmask,skernel);
                        end

                        % save
                        mask(:,:,mi) = cmask;
                    else
                        crowd(mi) = 1;
                    end
                end
                mask = mask(:,:,~crowd);
            else
                mask = zeros(y,x,1);
            end
            % sum and max at 2
            mask = min(2,nansum(mask,3));

            % combine images
            dI = double(cimg);
            dI = dI - 255/2;
            dI_mask = dI .* repmat(mask,1,1,3);
            % percentage combination
            dI = (1-mask_perc)*dI + mask_perc*dI_mask;
            dI = dI + 255;
            dI = dI / 2;
            mimg = uint8(dI);

            % test code:
    %         figure(1);
    %         subplot(311); imagesc(cimg);
    %         subplot(312); imagesc(mask);
    %         subplot(313); imagesc(mimg);
    %         pause(1);

            mimgs{ii} = mimg;
        end
        
        acat = sprintf('attend_%s',categories{aci});
        info.mimgs.(acat).(categories{ci}) = mimgs;
        disp(sprintf('Finished attend %s for images with %s',categories{aci},categories{ci}));
        
    end
    
end

%% Save out all images

save(fullfile('~/proj/fbsear/data/info_m.mat'),'info');