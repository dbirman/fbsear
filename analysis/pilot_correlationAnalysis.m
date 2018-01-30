function pilot_correlationAnalysis(cfolder,version)

%% Build the correlation overlay for the pilot experiment
% Five concat scans exist, one for each group, that are the correctly
% ordered individual scans. We have to pull these for all voxels and then
% in parallel build the correlation matrices. Then we build overlays to
% display the main effect comparisons (difference of
% corr(lookPeople,conPeople) and corr(lookPeople,conCars) and vice versa).
% The expectation is that look for people contrast people will be more
% correlated than the other contrast images (or fixation images). 

% This code performs the analysis step: pull the data for the whole brain
% (or individual ROIs), perform the cross correlation 

%% Move to Data WD
mrQuit
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

%% Setup a view + Load Concatenation
view = newView();
view = viewSet(view,'curGroup','Concatenation');

%%
allROIs = {'V1'};

clear rois desc

for curScan = 1:5
    desc{curScan} = viewGet(view,'description',curScan); % get description
%     rois{curScan} = loadROITSeries(view,allROIs,curScan,view.curGroup,'keepNAN=true');
    rois{curScan} = loadTSeries(view,curScan,[],[],[],[],[],view.curGroup);
    % reshape
    s = size(rois{curScan});
    rois{curScan} = reshape(rois{curScan},s(1)*s(2)*s(3),s(4));
end

%% Re-organize data

% roiData = zeros(size(rois{1}.tSeries,1),size(rois{1}.tSeries,2),5);
% 
% for curScan = 1:5
%     roiData(:,:,curScan) = rois{curScan}.tSeries;
% end

data = zeros(size(rois{1},1),size(rois{1},2),5);

for curScan = 1:5
    data(:,:,curScan) = rois{curScan};
end

%% Compute correlation matrix

corrMat = zeros(size(roiData,1),5,5);
parfor n = 1:size(data,1)
%     corrMat(n,:,:) = corrcoef(squeeze(roiData(n,:,:)));
    corrMat(n,:,:) = corrcoef(squeeze(data(n,:,:)));
end

%% Check fixation values
% We will remove all voxels that are not in the top 50% of average R^2
% values with the fixation task (i.e. that don't show a reasonable
% correlation with the fixation data, which would probably suggest that
% they were bogus voxels and not task responsive)
fixate = corrMat(:,1,2:5);
fixate = squeeze(mean(fixate,3));
% take only the top 5%
cutoff = quantile(fixate,.95);

fixidxs = fixate>cutoff;

%% Compute the average correlation matrix
avgCorrMat = squeeze(mean(corrMat(fixidxs,:,:)));
% set the diagonal to zero
avgCorrMat = avgCorrMat.*(~diag(ones(1,5)));
imagesc(avgCorrMat);
colormap('gray');
colorbar

%% Compute the difference scores
% We will compute three difference scores:
%   (1) Whether contrast runs look more like fixating or looking (ground
%   truth test)
%   (2) Whether look people is more correlated to contrast people, or
%   contrast cars
%   (3) Whether look cars is more correlated to contrast cars, or contrast
%   people

% Compute Look people
% row 2 (look people)
% row 4 - row 5 (contrast people - contrast cars)
peopleDiff = corrMat(:,2,4)-corrMat(:,2,5);

% Compute Look people
% row 3 (look people)
% row 4 - row 5 (contrast cars - contrast people)
carDiff = corrMat(:,3,5)-corrMat(:,3,4);

%% Load MLR and display these overalsy

mrQuit;
mrLoadRet;
%%
view = getMLRView;
% viewSet(view,'curGroup','Concatenatation');

peopleDiff = reshape(peopleDiff,s(1),s(2),s(3));
carDiff = reshape(carDiff,s(1),s(2),s(3));

mrDispOverlay(peopleDiff,2,view.curGroup,getMLRView,'overlayName=peopleDiff');
mrDispOverlay(carDiff,2,view.curGroup,getMLRView,'overlayName=carDiff');

% diff2
diff2 = peopleDiff - carDiff;
mrDispOverlay(diff2,2,view.curGroup,getMLRView,'overlayName=People_minus_Cars');