%This script is mostly copypasta from analysis.m to put together the bar
%graphs for the main images of squares figure in the paper

% Data format is: subj number, gender, age, image number, Y, X, reaction time, stop time (relative to new years 1970 in ms)
% Each new line is a new tap.

% Matlab file reading causes the new lines to run together, so the subject
% numbers and stop times in the subj variable below are formatted wierdly.

data_dir = '../subject_data/';

%% get all the data
Nsubj = 252;
subj = cell(Nsubj,1);
for k = 1:Nsubj
    if k>=100
        fname = [data_dir 'subject_0' num2str(k) '.txt'];
    elseif k>=10
        fname = [data_dir 'subject_00' num2str(k) '.txt'];
    else
        fname = [data_dir 'subject_000' num2str(k) '.txt'];
    end
    fid = fopen(fname);
    data = fread(fid, '*char')'; %read all contents into data as a char array (don't forget the `'` to make it a row rather than a column).
    fclose(fid);
    subj{k} = regexp(data, ',', 'split'); %This will return a cell array with the individual entries for each string you have between the commas.
    subj{k} = reshape(subj{k}(2:end),7,[]);
    
end

%% Organize onto the images (only first taps this time)
first_tap = true;

img_list = 1:80;
Nimg = length(img_list);
image_data = repmat(struct('Xdata',[],'Ydata',[],'RT',[],'subj',[],'tapnum',[],'stop_time',[],'start_time',[]),80,1);

for k = 1:Nsubj

    if first_tap;
        m = 2;
    else
        m = size(subj{k},2);
    end
    
    for l = 1:m
        im_idx = str2double(subj{k}{3,l});
        image_data(im_idx).subj(end+1) = k;
        image_data(im_idx).tapnum(end+1) = l; %NOTE: tapnum is messed up for the subject that had two data points drop out
        image_data(im_idx).Ydata(end+1)     = str2double(subj{k}{4,l});
        image_data(im_idx).Xdata(end+1)     = str2double(subj{k}{5,l});
        image_data(im_idx).RT(end+1)        = str2double(subj{k}{6,l});
        image_data(im_idx).stop_time(end+1) = str2double(subj{k}{7,l}(1:13));
        image_data(im_idx).start_time(end+1) = image_data(im_idx).stop_time(end) ...
                                               - image_data(im_idx).RT(end);
    end
end

% Data for all images together
X_all = [image_data(:).Xdata];
Y_all = [image_data(:).Ydata];
RT_all = [image_data(:).RT];
start_all = [image_data(:).start_time];
stop_all = [image_data(:).stop_time];
tapnum_all = [image_data(:).tapnum];

%% Compare images of all black (red) with gray (pink) (First taps)

% img_list = 1:5; %Black/gray comparison
% img_list = 21:25; %Red/pink comparison
img_list = [1:5 21:25];

offset = 5; %offset between spatially matching images (number of images in the group)

% Set up for the contingency table
%
%       insquare   elsewhere  total
% faint   a        c          K
% solid   b        d          -
% total   N        -          M
a = zeros(length(img_list),1); 
b = zeros(length(img_list),1); 
M = zeros(length(img_list),1);
K = zeros(length(img_list),1);
N = zeros(length(img_list),1);
p = ones(length(img_list),1);
insquare_faint_all = [];
insquare_solid_all = [];
for k = 1:length(img_list)
    %read image
    im = imread(['../images/squares' num2str(img_list(k)) '.png']);
    im_in = double(sum(im,3));
    
    %create mask for faint square
    mask = im_in>255 & im_in<255*3;
    mask(1005:1024,:) = 0;
    
    %Find taps in faint square
    X = image_data(img_list(k)).Xdata;
    Y = image_data(img_list(k)).Ydata;
    insquare_faint = zeros(length(X),1);
    for l = 1:length(X)
        insquare_faint(l) = mask(Y(l),X(l));
    end
    
    %Find taps in corresponding non-faint square
    X = image_data(img_list(k)+offset).Xdata;
    Y = image_data(img_list(k)+offset).Ydata;
    insquare_solid = zeros(length(X),1);
    for l = 1:length(X)
        insquare_solid(l) = mask(Y(l),X(l));
    end

    %plot
    figure(1);
    subplot(121);
    hold all;
    pad = 10;
    im(1005:1024,:) = 255; 
    im = [ones(pad,size(im,2)+2*pad,3); ones(size(im,1),pad,3) im ones(size(im,1),pad,3); 
          ones(pad,size(im,2)+2*pad,3)];
    imshow(imrotate(im,-90));
    scatter(1024-image_data(img_list(k)).Ydata+pad,image_data(img_list(k)).Xdata+pad,[],'g','filled');
    title([num2str(sum(insquare_faint)) ' out of ' num2str(length(image_data(img_list(k)).Xdata))]);
    
    subplot(122);
    hold all;
    im_offset = imread(['../images/squares' num2str(img_list(k)+offset) '.png']);
    im_offset(1005:1024,:) = 255;
    im_offset = [ones(pad,size(im_offset,2)+2*pad,3); 
                 ones(size(im_offset,1),pad,3) im_offset ones(size(im_offset,1),pad,3);
                 ones(pad,size(im_offset,2)+2*pad,3)];
    imshow(imrotate(im_offset,-90));
    scatter(1024-image_data(img_list(k)+offset).Ydata+pad,image_data(img_list(k)+offset).Xdata+pad,[],'g','filled');
    title([num2str(sum(insquare_solid)) ' out of ' num2str(length(image_data(img_list(k)+offset).Xdata))]);
    
    boldify
    % Set up for the contingency table
    %
    %       insquare   elsewhere  total
    % faint   a        c          K
    % solid   b        d          -
    % total   N        -          M
    a(k) = sum(insquare_faint);
    b(k) = sum(insquare_solid);
    M(k) = length(image_data(img_list(k)).Xdata)+length(image_data(img_list(k)+offset).Xdata);
    K(k) = length(image_data(img_list(k)).Xdata);
    N(k) = a(k)+b(k);
    
    if all( image_data(img_list(k)).subj == image_data(img_list(k)+offset).subj);
        [~,p(k)] = ttest(insquare_faint,insquare_solid);
        insquare_faint_all = [insquare_faint_all; insquare_faint];
        insquare_solid_all = [insquare_solid_all; insquare_solid];
    else
        p(k) = fexact( a(k),M(k),K(k),N(k), 'tail','r');
    end
    
    disp(['p = ' num2str(p(k))]);
%     pause;
    
end

if ~isempty(insquare_faint_all)
    [~, ptot] = ttest(insquare_faint_all,insquare_solid_all);
    disp(['p = ' num2str(ptot)]);
else
    ptot = fexact(sum(a),sum(M),sum(K),sum(N),'tail','r');
    disp(['p = ' num2str(ptot)]);
end

% This is super messy code-wise, but you have to run this using the
% firsttaps flag as both true and false (regenerating image_data each time)
% without closing the figure. 

%% Put together the figures for First tap data

[phat, pci] = binofit([sum(a(1:5)) sum(b(1:5)) sum(a(6:10)) sum(b(6:10))],[sum(K(1:5)) (sum(M(1:5))-sum(K(1:5))) sum(K(6:10)) (sum(M(6:10))-sum(K(6:10)))]);
err = sqrt(phat.*(1-phat))./sqrt([sum(K(1:5)) (sum(M(1:5))-sum(K(1:5))) sum(K(6:10)) (sum(M(6:10))-sum(K(6:10)))]);


% inx = (1:2)+2*(~first_tap);
% if ~first_tap
%     hold on;
% end
figure(2);
inx = 1:2;
bar(1,phat(inx(1)));
hold on;
bar(2,phat(inx(2)));
errorbar(1:2,phat(inx),err(inx),'k','linestyle','none','linewidth',2);
hold off;
title('Gray/Black First Tap','FontSize',20);
xlim([0.5 2.5]);
ylim([0 0.5]);
ylabel('Tap rate','FontSize',20)
set(gca,'FontSize',20)
set(gca,'XTick',[1 2])
set(gca,'XTicklabel','Gray Singleton|Non-singleton')

figure(3);
inx = 3:4;
bar(1,phat(inx(1)));
hold on;
bar(2,phat(inx(2)));
errorbar(1:2,phat(inx),err(inx),'k','linestyle','none','linewidth',2);
hold off;
title('Pink/Red First Tap','FontSize',20);
xlim([0.5 2.5]);
ylim([0 0.5]);
ylabel('Tap rate','FontSize',20)
set(gca,'FontSize',20)
set(gca,'XTick',[1 2])
set(gca,'XTicklabel','Pink Singleton|Non-singleton')
% bar(inx,[sum(a(1:5))./sum(K(1:5)) sum(b(1:5))./(sum(M(1:5))-sum(K(1:5))); sum(a(6:10))./sum(K(6:10)) sum(b(6:10))./(sum(M(6:10))-sum(K(6:10)))])
% [phat, pci] = binofit([sum(a(1:5)) sum(b(1:5)) sum(a(6:10)) sum(b(6:10))],[sum(K(1:5)) (sum(M(1:5))-sum(K(1:5))) sum(K(6:10)) (sum(M(6:10))-sum(K(6:10)))]);
% err = sqrt(phat.*(1-phat))./sqrt([sum(K(1:5)) (sum(M(1:5))-sum(K(1:5))) sum(K(6:10)) (sum(M(6:10))-sum(K(6:10)))]);
% x = [inx(1)-1/7 inx(1)+1/7 inx(2)-1/7 inx(2)+1/7];
% errorbar(x,phat,err,'k','linestyle','none');
% hold off;


%% Organize onto the images (All taps this time)
first_tap = false;

img_list = 1:80;
Nimg = length(img_list);
image_data = repmat(struct('Xdata',[],'Ydata',[],'RT',[],'subj',[],'tapnum',[],'stop_time',[],'start_time',[]),80,1);

for k = 1:Nsubj

    if first_tap;
        m = 2;
    else
        m = size(subj{k},2);
    end
    
    for l = 1:m
        im_idx = str2double(subj{k}{3,l});
        image_data(im_idx).subj(end+1) = k;
        image_data(im_idx).tapnum(end+1) = l; %NOTE: tapnum is messed up for the subject that had two data points drop out
        image_data(im_idx).Ydata(end+1)     = str2double(subj{k}{4,l});
        image_data(im_idx).Xdata(end+1)     = str2double(subj{k}{5,l});
        image_data(im_idx).RT(end+1)        = str2double(subj{k}{6,l});
        image_data(im_idx).stop_time(end+1) = str2double(subj{k}{7,l}(1:13));
        image_data(im_idx).start_time(end+1) = image_data(im_idx).stop_time(end) ...
                                               - image_data(im_idx).RT(end);
    end
end

% Data for all images together
X_all = [image_data(:).Xdata];
Y_all = [image_data(:).Ydata];
RT_all = [image_data(:).RT];
start_all = [image_data(:).start_time];
stop_all = [image_data(:).stop_time];
tapnum_all = [image_data(:).tapnum];

%% Compare images of all black (red) with gray (pink) (All taps)

% img_list = 1:5; %Black/gray comparison
% img_list = 21:25; %Red/pink comparison
img_list = [1:5 21:25];

offset = 5; %offset between spatially matching images (number of images in the group)

% Set up for the contingency table
%
%       insquare   elsewhere  total
% faint   a        c          K
% solid   b        d          -
% total   N        -          M
a = zeros(length(img_list),1); 
b = zeros(length(img_list),1); 
M = zeros(length(img_list),1);
K = zeros(length(img_list),1);
N = zeros(length(img_list),1);
p = ones(length(img_list),1);
insquare_faint_all = [];
insquare_solid_all = [];
for k = 1:length(img_list)
    %read image
    im = imread(['../images/squares' num2str(img_list(k)) '.png']);
    im_in = double(sum(im,3));
    
    %create mask for faint square
    mask = im_in>255 & im_in<255*3;
    mask(1005:1024,:) = 0;
    
    %Find taps in faint square
    X = image_data(img_list(k)).Xdata;
    Y = image_data(img_list(k)).Ydata;
    insquare_faint = zeros(length(X),1);
    for l = 1:length(X)
        insquare_faint(l) = mask(Y(l),X(l));
    end
    
    %Find taps in corresponding non-faint square
    X = image_data(img_list(k)+offset).Xdata;
    Y = image_data(img_list(k)+offset).Ydata;
    insquare_solid = zeros(length(X),1);
    for l = 1:length(X)
        insquare_solid(l) = mask(Y(l),X(l));
    end

    %plot
    figure(1);
    subplot(121);
    hold all;
    pad = 10;
    im(1005:1024,:) = 255; 
    im = [ones(pad,size(im,2)+2*pad,3); ones(size(im,1),pad,3) im ones(size(im,1),pad,3); 
          ones(pad,size(im,2)+2*pad,3)];
    imshow(imrotate(im,-90));
    scatter(1024-image_data(img_list(k)).Ydata+pad,image_data(img_list(k)).Xdata+pad,[],'g','filled');
    title([num2str(sum(insquare_faint)) ' out of ' num2str(length(image_data(img_list(k)).Xdata))]);
    
    subplot(122);
    hold all;
    im_offset = imread(['../images/squares' num2str(img_list(k)+offset) '.png']);
    im_offset(1005:1024,:) = 255;
    im_offset = [ones(pad,size(im_offset,2)+2*pad,3); 
                 ones(size(im_offset,1),pad,3) im_offset ones(size(im_offset,1),pad,3);
                 ones(pad,size(im_offset,2)+2*pad,3)];
    imshow(imrotate(im_offset,-90));
    scatter(1024-image_data(img_list(k)+offset).Ydata+pad,image_data(img_list(k)+offset).Xdata+pad,[],'g','filled');
    title([num2str(sum(insquare_solid)) ' out of ' num2str(length(image_data(img_list(k)+offset).Xdata))]);
    
    boldify
    % Set up for the contingency table
    %
    %       insquare   elsewhere  total
    % faint   a        c          K
    % solid   b        d          -
    % total   N        -          M
    a(k) = sum(insquare_faint);
    b(k) = sum(insquare_solid);
    M(k) = length(image_data(img_list(k)).Xdata)+length(image_data(img_list(k)+offset).Xdata);
    K(k) = length(image_data(img_list(k)).Xdata);
    N(k) = a(k)+b(k);
    
    if all( image_data(img_list(k)).subj == image_data(img_list(k)+offset).subj);
        [~,p(k)] = ttest(insquare_faint,insquare_solid);
        insquare_faint_all = [insquare_faint_all; insquare_faint];
        insquare_solid_all = [insquare_solid_all; insquare_solid];
    else
        p(k) = fexact( a(k),M(k),K(k),N(k), 'tail','r');
    end
    
    disp(['p = ' num2str(p(k))]);
%     pause;
    
end

if ~isempty(insquare_faint_all)
    [~, ptot] = ttest(insquare_faint_all,insquare_solid_all);
    disp(['p = ' num2str(ptot)]);
else
    ptot = fexact(sum(a),sum(M),sum(K),sum(N),'tail','r');
    disp(['p = ' num2str(ptot)]);
end

%% Put together the figures for All tap data

[phat, pci] = binofit([sum(a(1:5)) sum(b(1:5)) sum(a(6:10)) sum(b(6:10))],[sum(K(1:5)) (sum(M(1:5))-sum(K(1:5))) sum(K(6:10)) (sum(M(6:10))-sum(K(6:10)))]);
err = sqrt(phat.*(1-phat))./sqrt([sum(K(1:5)) (sum(M(1:5))-sum(K(1:5))) sum(K(6:10)) (sum(M(6:10))-sum(K(6:10)))]);


figure(4);
inx = 1:2;
bar(1,phat(inx(1)));
hold on;
bar(2,phat(inx(2)));
errorbar(1:2,phat(inx),err(inx),'k','linestyle','none','linewidth',2);
hold off;
title('Gray/Black All Taps','FontSize',20);
xlim([0.5 2.5]);
ylim([0 0.5]);
ylabel('Tap rate','FontSize',20)
set(gca,'FontSize',20)
set(gca,'XTick',[1 2])
set(gca,'XTicklabel','Gray Singleton|Non-singleton')

figure(5);
inx = 3:4;
bar(1,phat(inx(1)));
hold on;
bar(2,phat(inx(2)));
errorbar(1:2,phat(inx),err(inx),'k','linestyle','none','linewidth',2);
hold off;
title('Pink/Red All Taps','FontSize',20);
xlim([0.5 2.5]);
ylim([0 0.5]);
ylabel('Tap rate','FontSize',20)
set(gca,'FontSize',20)
set(gca,'XTick',[1 2])
set(gca,'XTicklabel','Pink Singleton|Non-singleton')