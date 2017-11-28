%% Test
img = coco.loadImgs(581061);
I = imread(fullfile('~/proj/fbsear/',sprintf('images/%s/%s',dataType,img.file_name)));
imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[]);

%% Take image and randomly gaussian sample from it
[xx,yy] = meshgrid(1:size(I,2),1:size(I,1));

n = 100;
gfilts = zeros(size(I,1),size(I,2),size(I,3),n);
crop = cell(1,n);
for g = 1:n
    x = randi(size(I,2));
    y = randi(size(I,1));
    sd = 50+randi(25);
    dist = hypot(xx-x,yy-y);
    gauss = normpdf(dist,0,sd);
    gauss = gauss ./ max(gauss(:));
    fwhm = ceil(sd*2);
    cropx = [max(1,x-fwhm) min(size(I,2),x+fwhm)];
    cropy = [max(1,y-fwhm) min(size(I,1),y+fwhm)];
    gfilts(:,:,:,g) = repmat(gauss,1,1,3) .* double(I);
    crop{g} = gfilts(cropy(1):cropy(2),cropx(1):cropx(2),:,g);
end

%% Look

for i = 1:n
%     imagesc(uint8(gfilts(:,:,:,i))); axis('image'); set(gca,'XTick',[],'YTick',[]);
    imagesc(uint8(crop{i})); axis('image'); set(gca,'XTick',[],'YTick',[]);
    pause(0.1);
end

%% Noise
% take the random crops and position them randomly on top of each other,
% then average
s = size(gfilts);
s(1) = s(1) + fwhm*2;
s(2) = s(2) + fwhm*2;
noise = nan(s);

for g=1:n
    for r=1:1 % 10 repeats
        % pick a position
        x = randi(size(noise,2)-size(crop{g},2));
        y = randi(size(noise,1)-size(crop{g},1));
        noise(y:(y+size(crop{g},1)-1),x:(x+size(crop{g},2)-1),:,g) = crop{g};
    end
end

noise = noise((fwhm+1):(s(1)-fwhm),(fwhm+1):(s(2)-fwhm),:,:);
%% Average
noise_ = nanmean(noise,4);
noise_ = noise_./max(noise_(:))*255;

%%
imagesc(uint8(noise_))

%% Histogram equalization
for rgb = 1:3
    i = I(:,:,rgb);
    h = hist(double(i(:)),0:1:255);
    I2(:,:,rgb) = histeq(I(:,:,rgb),h);
    N2(:,:,rgb) = histeq(uint8(noise_(:,:,rgb)),h);
end
figure; imagesc(uint8(I2));
figure; imagesc(uint8(N2));
%% Noise look
for i = 1:n
    imagesc(uint8(noise(:,:,:,i))); axis('image'); set(gca,'XTick',[],'YTick',[]);
%     imagesc(uint8(crop{i})); axis('image'); set(gca,'XTick',[],'YTick',[]);
    pause(0.1);
end

%% Generate combined final image (histogram equalized and combined according to object position)

x = 180;
y = 220;
sd = 60;
dist = hypot(xx-x,yy-y);
gauss = normpdf(dist,0,sd);
gauss = gauss ./ max(gauss(:));

final = repmat(1-gauss,1,1,3).*N2 + repmat(gauss,1,1,3).*double(I2);
figure;
imagesc(uint8(final));

%% More and more noise

figure(37);
for i = 0.5:0.01:1
    gauss = gauss*(1-i);
    imagesc(uint8(repmat(1-gauss,size(I,1),size(I,2),3).*N2 + repmat(gauss,size(I,1),size(I,2),3).*double(I2)));
    pause(0.1);
end