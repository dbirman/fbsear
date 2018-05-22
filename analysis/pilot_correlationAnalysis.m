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
% data = roiData;

data = zeros(size(rois{1},1),size(rois{1},2),5);

for curScan = 1:5
    data(:,:,curScan) = rois{curScan};
end

%% Compute correlation matrix

corrMat = zeros(size(data,1),5,5);
parfor n = 1:size(data,1)
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

%% Compute the average correlation matrix

h = figure;
qs = [0.5 0.95 0.99];

for qi = 1:length(qs)
    subplot(length(qs),1,qi);
    
    cutoff = quantile(fixate,qs(qi));

    fixidxs = fixate>cutoff;
    avgCorrMat = squeeze(mean(corrMat(fixidxs,:,:)));
    % set the diagonal to zero
    avgCorrMat = triu(avgCorrMat.*(~diag(ones(1,5))));
    imagesc(avgCorrMat);
    colormap('gray');
    colorbar
    set(gca,'XTick',1:5,'XTickLabel',desc,'XAxisLocation','top');
    set(gca,'YTick',1:5,'YTickLabel',desc);
end
savepdf(h,fullfile('~/proj/fbsear/analysis/figs','correlation_mat.pdf'));

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
fixate = reshape(fixate,s(1),s(2),s(3));

mrDispOverlay(peopleDiff,2,view.curGroup,getMLRView,'overlayName=peopleDiff');
mrDispOverlay(carDiff,2,view.curGroup,getMLRView,'overlayName=carDiff');
mrDispOverlay(fixate,2,view.curGroup,getMLRView,'overlayName=fix_cutoff');

%% Generate plots
% First generate plots of the look people corr people/car con histograms
% after using the fixation cutoffs (>.95 seems to work well)

h = figure;

corrMatFix = corrMat(fixate>quantile(fixate(:),.95),:,:);

lp_x_cp = corrMatFix(:,2,4);
lp_x_cc = corrMatFix(:,2,5);
lc_x_cp = corrMatFix(:,3,4);
lc_x_cc = corrMatFix(:,3,5);
vl = [min([lp_x_cc(:) ;lc_x_cc(:)]) max([lp_x_cp(:) ;lc_x_cp(:)])];
bin_centers = vl(1):(diff(vl)/10):vl(2);

% plot marginals
lp_x_cp_hist = hist(lp_x_cp,bin_centers);
subplot(3,6,[1 2]);
bar(bin_centers,lp_x_cp_hist,'FaceColor','k');
a = axis;
axis([vl a(3) a(4)]);
a = axis;
drawPublishAxis;

lp_x_cc_hist = hist(lp_x_cc,bin_centers);
subplot(3,6,[9 15]);
barh(bin_centers,lp_x_cc_hist,'FaceColor','k');
axis([a(3) a(4) a(1) a(2)]);
drawPublishAxis;


lc_x_cp_hist = hist(lc_x_cp,bin_centers);
subplot(3,6,[4 5]);
bar(bin_centers,lc_x_cp_hist,'FaceColor','k');
axis(a);
drawPublishAxis;

lc_x_cc_hist = hist(lc_x_cc,bin_centers);
subplot(3,6,[12 18]);
barh(bin_centers,lc_x_cc_hist,'FaceColor','k');
axis([a(3) a(4) a(1) a(2)]);
drawPublishAxis;

% plot correlation plot
subplot(3,6,[7 8 13 14]); hold on
plot(lp_x_cc,lp_x_cp,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',1);
plot(vl,vl,'--r');
title('Look for people (selected for visual-responsive voxels)');
xlabel('Corr with contrast cars');
ylabel('Corr with contrast people');
axis([vl vl]);
drawPublishAxis;

subplot(3,6,[10 11 16 17]); hold on
plot(lc_x_cc,lc_x_cp,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',1);
vl = [min(lc_x_cc(:)) max(lc_x_cp(:))];
plot(vl,vl,'--r');
title('Look for cars (selected for visual-responsive voxels)');
xlabel('Corr with contrast cars');
ylabel('Corr with contrast people');
axis([vl vl]);
drawPublishAxis;

savepdf(h,fullfile('~/proj/fbsear/analysis/figs','correlation_comp.pdf'));