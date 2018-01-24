%% Load imdata
if isfile(fullfile('~/proj/fbsear/imdata_outdoor.mat'))
    load(fullfile('~/proj/fbsear/imdata_outdoor.mat'));
else
    load(fullfile('~/proj/fbsear/imdata.mat'));
end
%% Info
% num categories
disp(sum(imdata_(:,3:11),1));

%% If no 12th column add it

if size(imdata_,2)<12, imdata_(:,12) = nan; end

%% Start indoor/outdoor code
figure;

% find all nan values
nan_idxs = isnan(imdata_(:,12));

% Count outdoor per category and remove the top two categories that have the
% most
for i = 3:11
    cat_indexes{i-2} = find(nan_idxs .* imdata_(:,i)==1);
    dat = imdata_(cat_indexes{i-2},:);
    count(i-2) = nansum(dat(:,12)==1);
end

max_c = sort(count);
aidxs = [];
for i = 1:9
%     if ~((count(i)==max_c(end)) || (count(i)==max_c(end-1)) || (count(i)==max_c(end-2)) || (count(i)==max_c(end-3)) || (count(i)==max_c(end-4)) || (count(i)==max_c(end-5)) || (count(i)==max_c(end-6)))
        % remove
    cidxs = cat_indexes{i};
    cidxs = cidxs(randperm(length(cidxs)));
    disp(length(cidxs));
    aidxs = [aidxs cidxs(1:min(length(cidxs),10))'];
%     end
end

aidxs = aidxs(randperm(length(aidxs)));

for i = 1:length(aidxs)
    if isnan(imdata_(aidxs(i),12))
        img = coco.loadImgs(imdata_(aidxs(i),1));
        I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
        imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[]);
        pause(.25);
        key = waitForKeypress;

        if strcmp(key,'o')
            disp('outdoor');
            imdata_(aidxs(i),12)=1;
        elseif strcmp(key,'i')
            disp('indoor');
            imdata_(aidxs(i),12)=0;
        elseif strcmp(key,'c')
            disp('correction');
            imdata_(aidxs(i-1),12)=nan;
        else
            disp('ignore');
        end
    end
%     for i = 3:11
%         dat = imdata_(imdata_(:,i)==1,:);
%         disp(nansum(dat(:,12)==1));
%     end

end

%% Count outdoor per category
for i = 3:11
    dat = imdata_(imdata_(:,i)==1,:);
    disp(nansum(dat(:,12)==1));
end

%% Save
save(fullfile('~/proj/fbsear/imdata_outdoor.mat'),'imdata_');

%% Display images from all categories

imdata_select = imdata_(imdata_(:,12)==1,:);

% disp(sprintf('For base %s I found %i images',[base{category}{:}],length(imgIds)));
figure(1);
% 4x4 array of random images
n = length(imgIds);
for i = 1:9
    imdata_cat = imdata_select(imdata_select(:,2+i)==1,:);
    subplot(3,3,i);
    img = coco.loadImgs(imdata_cat(randi(size(imdata_cat,1)),1));
    I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
    imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[])
    title(sprintf('%s',[base{i}{:}]));
end