function pilot_semanticAnalysis(cfolder,version)

%% Build the semantic overlay for the pilot experiment
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
cd(fullfile('~/data/fbsear_fixationpilot/',cfolder));
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

clear rois desc

desc = viewGet(view,'description',1); % get description
concatInfo = viewGet(view,'concatInfo');
stimfiles = viewGet(view,'stimfile');


% allROIs = {'V1'};
% scan = loadROITSeries(view,allROIs,1,view.curGroup,'keepNAN=true');

scan = loadTSeries(view,1,[],[],[],[],[],view.curGroup);
    
% reshape
% s = size(rois{curScan});
% rois{curScan} = reshape(rois{curScan},s(1)*s(2)*s(3),s(4));

%% Go through the stimfiles and check that group info is available (first pilot didn't include this)
if ~isfield(stimfiles{1}.stimulus.curRun,'group')
    for si = 1:length(stimfiles)
        for run = 1:4
            if ~isempty(strfind(stimfiles{si}.stimulus.curRun.text,num2str(run)))
                stimfiles{si}.stimulus.curRun.group = run; break
            end
        end
    end
end

%% Compute the r^2 on repeated runs for all voxels

scan = reshape(scan,size(scan,1)*size(scan,2)*size(scan,3),size(scan,4));

%   group    |    repeat    |   voxels    | timeseries
repdata = zeros(4,2,size(scan,1),size(scan,2)/length(stimfiles));

for group = 1:4
    for repeat = 1:2
        for si = 1:length(stimfiles)
            if (stimfiles{si}.stimulus.curRun.group==group) && (stimfiles{si}.stimulus.curRun.repeat==repeat)
                rt = concatInfo.runTransition(si,:);
                
                repdata(group,repeat,:,:) = scan(:,rt(1):rt(2));
                
            end
        end
    end
end

maxr2 = zeros(4,size(repdata,3));

for group = 1:4
    % compute the correlation matrix 
    data = squeeze(repdata(group,:,:,:));
    parfor vi = 1:size(data,2)
        c = corr(squeeze(data(:,vi,:))');
        maxr2(group,vi) = c(1,2);
    end
end

maxr2 = maxr2.^2;

%% Average across folds
maxr2_ = sort(nanmean(maxr2));
maxr2_ = maxr2_(~isnan(maxr2_));
%% Plot r2 quantiles
figure;
qs = 0.9:.001:1;
plot(qs*length(maxr2_),quantile(maxr2_,qs));


%% Missing: plot r^2 average across groups on surfaces

%% Semantic analysis

%% Get annotation data
clear aAnns
for si = 1:length(stimfiles)
    cAnns = stimfiles{si}.stimulus.curRun.anns;
    % multiply the spacing by three to account for TRs
    cAnns_ = zeros(size(cAnns,1)*3,size(cAnns,2));
    for ci = 1:size(cAnns,2)
        temp = cAnns(:,ci);
        temp = repmat(temp,1,3)'; 
        temp = temp(:);
        cAnns_(:,ci) = temp;
    end
    % flag if the scan length is different
    runTimes = concatInfo.runTransition(si,:);
    runLength = runTimes(2)-runTimes(1)+1;
    if runLength<size(cAnns_,1)
        warning('Runs were shorter than the stimulus -- this is a serious issue');
    end
    if runLength~=size(cAnns_,1)
        warning(sprintf('Run %i was length %i and stimulus length %i -- this could be important',si,runLength,size(cAnns_,1)));
    end
    % pad the end if the runLength is longer
    cAnns_(end:(end+runLength-size(cAnns_,1)),:) = 0;
    % copy data in
    aAnns(si,:,:) = cAnns_;
end

%% Build the design matrix for the 80 categories
[params,hrf] = hrfDoubleGamma([],0.5);
hrf = hrf ./ max(hrf) / 3; % div three so that outputs are max 1

design = zeros(size(aAnns,1)*size(aAnns,2),size(aAnns,3));

l = size(aAnns,2);

% convolve each column
for si = 1:size(aAnns,1)
    for ci = 1:size(aAnns,3)
        dat = aAnns(si,:,ci);
        dat_ = conv(dat,hrf);
        
        % remove trailing
        dat_ = dat_(1:length(dat));
        
        design(((si-1)*l+1):((si)*l),ci) = dat_;
        
        % test code
%         figure; hold on
%         plot(dat,'-r');
%         plot(dat_);
    end
end

%% For each voxel compute a lasso 

tSeries = scan.tSeries;

for vi = 1:size(tSeries,1)
    y = tSeries(vi,:)';
    
    % lasso version
    [b,stats] = lasso(design,y-1,'CV',10);
    
    % partial least squares verion
%     [xl,yl,xs,ys,beta,pctvar,mse] = plsregress(design,y-1,10,'CV',10);
%     b = B(:,fitInfo.IndexMinMSE);
%     r(vi) = corr(design*b,y-1);
    %% test code
    figure; hold on
    plot(y(1:500)-1);
    pred = design*b;
    plot(pred(1:500),'-r');
end

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