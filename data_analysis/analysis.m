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

%% Organize onto the images
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

%% look for double taps
thresh = 400;
sum(RT_all<thresh)
[sortstart, inx] = sort(start_all);
revsort = 1:length(inx); revsort(inx) = revsort; %reverse sorting indices
sortstop = stop_all(inx);
TR_sorted = [1E6 sortstart(2:end)-sortstop(1:end-1)];
TR_all = TR_sorted(revsort);
sum(TR_all<thresh)

sel = RT_all<thresh;
scatter(1024-Y_all(sel),X_all(sel)); xlim([0 1024]);ylim([0 768]); axis ij

% tapdiff_all = diff(sort([start_all stop_all]));
% sum(tapdiff_all<thresh)
%% Analyze taps by subject over images

mux_all = zeros(1,Nsubj);
muy_all = zeros(1,Nsubj);
muRT_all = zeros(1,Nsubj);
sig_all = zeros(2,2,Nsubj);
age_all = zeros(1,Nsubj);
ismale_all = zeros(1,Nsubj);

for k = 1:Nsubj
%     for l = 1:size(subj{k},2);
%         
%         im_idx = str2double(subj{k}{3,l});
%         image_data(im_idx).Ydata(end+1) = str2double(subj{k}{4,l});
%         image_data(im_idx).Xdata(end+1) = str2double(subj{k}{5,l});
%         image_data(im_idx).RT(end+1)    = str2double(subj{k}{6,l});
%     end
    s = subj{k};
    X = str2double(s(5,:));
    Y = str2double(s(4,:));
    RT = str2double(s(6,:));
    age = str2double(s{2,1}(1:2));
    age_all(k) = age;
    ismale = strcmp(s{1,1}(1),'m');
    ismale_all(k) = ismale;
    
    scatter(1024-Y,X);
    xlim([0 1024]);ylim([0 768]);
    axis ij; axis equal;
    title(['Subject ' num2str(k)]);
    
    mux = mean(X); 
    muy = mean(Y);
    muRT = mean(RT);
    hold all;
    scatter(1024-muy,mux,'r');
    var = [X-mux; Y-muy]*[X-mux; Y-muy]'/(length(X)-1);
    sig = sqrtm(var);
    sig_all(:,:,k) = sig;
    mux_all(k) = mux;
    muy_all(k) = muy;
    muRT_all(k) = muRT;
    hold off;
%     pause

end

% hist(squeeze(sig_all(1,1,:)),50:100:1000)
% hist(squeeze(sig_all(2,2,:)),50:100:1000)
% hist(squeeze(sig_all(2,1,:)),-350:100:1000)
xsigall = squeeze(sig_all(1,1,:));
agelist = unique(age_all);
for k = 1:length(agelist)
    xsig_meanage(k) = mean(xsigall(age_all==agelist(k)));
end
scatter(age_all,squeeze(sig_all(1,1,:)))
hold all;
scatter(agelist,xsig_meanage,'r');
hold off;



% scatter(age_all,squeeze(sig_all(2,2,:)))

scatter(age_all,squeeze(sig_all(2,1,:)))

scatter(1024-muy_all,mux_all,[],age_all)
axis ij; axis equal;
xlim([0 1024]);ylim([0 768]);

%% Generate Demographics table
age_all_idx = zeros(size(age_all));
for k = 1:length(agelist)
    age_all_idx(age_all==agelist(k)) = k;
end

demographics = accumarray([ismale_all'+1, age_all_idx'],1,[2 length(agelist)]);
age_tot = sum(demographics,1);
gend_tot = sum(demographics,2);

bar(demographics');
legend('Male','Female')
xlabel('Age Group');
ylabel('Number of Subjects');
% Manually go in to change tick values to 18-22, 23-30, 31-40, 41-50, 51+

%% Analyze taps independent of image or subject

figure(1)
subplot(121);
scatter(1024-Y_all,X_all,[],RT_all/1000,'filled');
title('Results for all images');
xlim([0 1024]);ylim([0 768]);
axis ij; axis equal;

subplot(122);
RT_all_clip = RT_all;RT_all_clip(RT_all_clip>3000) = 3000;
scatter(1024-Y_all,X_all,[],RT_all_clip/1000,'filled');
title('Results for all images (RT clipped)');
xlim([0 1024]);ylim([0 768]);
axis ij; axis equal;

figure(2);
hist(RT_all/1000,100)
xlabel('Reaction Times (sec)');
ylabel('Count');


%% Show all images with taps overlaid
load('fixation_points.mat');
for k = 1:length(img_list)
	figure(1);
    im = imread(['../images/squares' num2str(img_list(k)) '.png']);
    if img_list(k)<=30 %pad images that are less than 1024 pixels wide for display purposes
        im(1005:1024,:,:)=255;
    end
    im = imrotate(im,-90);
    f = imshow(im);
    hold all;
    if img_list(k)<=30
        scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],'MarkerEdgeColor','g','MarkerFaceColor','g');
    else
        
        %grid lines
        M = size(im,1);
        N = size(im,2);

        binsize = 64;
        for l = 1:binsize:M
            x = [1 N];
            y = [l l];
            plot(x,y,'Color','w','LineStyle','-');
            plot(x,y,'Color','k','LineStyle',':');
        end

        for l = 1:binsize:N
            x = [l l];
            y = [1 M];
            plot(x,y,'Color','w','LineStyle','-');
            plot(x,y,'Color','k','LineStyle',':');
        end
        
        scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],'MarkerEdgeColor','r','MarkerFaceColor','r');
        scatter(fix_points_x{k-30},fix_points_y{k-30},[],'MarkerEdgeColor','c','MarkerFaceColor','c');
    end

    
    title(num2str(k))
    hold off;
    boldify
    
    
%     if first_tap
%         saveas(f,['./images_252subj/firsttap_' num2str(img_list(k)) '.png']);
%     else
%         saveas(f,['./images_252subj/alltaps_' num2str(img_list(k)) '.png']);
%     end

	pause
end

%% Show all images with taps colored by RT
load('fixation_points.mat');
for k = 1:length(img_list)
	figure(1);
    im = imread(['../images/squares' num2str(img_list(k)) '.png']);
    if img_list(k)<=30 %pad images that are less than 1024 pixels wide for display purposes
        im(1005:1024,:,:)=255;
    end
    im = imrotate(im,-90);
    f = imshow(im);
    hold all;
    if img_list(k)<=30
        scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],image_data(img_list(k)).RT/1000);
    else
        scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],image_data(img_list(k)).RT/1000);
        scatter(fix_points_x{k-30},fix_points_y{k-30},[],'MarkerEdgeColor','g','MarkerFaceColor','b');
    end
    caxis([0 5]);
    title(num2str(k))
    hold off;
    boldify
	pause

    
%     if first_tap
%         saveas(f,['./images_252subj/firsttap_' num2str(img_list(k)) '.png']);
%     else
%         saveas(f,['./images_252subj/alltaps_' num2str(img_list(k)) '.png']);
%     end

end

%% Show taps and interest data 
% Requires getting the interest map data "int_data_norm" first

% img_list = 31:80;
% list = [1:25 76:100];
% 
% for k = 1:length(img_list)
% figure(1)
% subplot(121)
% im = imread(['../images/squares' num2str(img_list(k)) '.png']);
%     im = permute(im,[2 1 3]);
%     im = im(:,end:-1:1,:);
%     imshow(im);
%     hold all;
%     scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],'g');
%     title(num2str(k))
%     hold off;
% subplot(122)
% imagesc(int_data_norm{list(k)});colorbar; colormap(gray);axis equal;
%     hold all;
%     scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],'g');
%     hold off;
% pause
% end

%% Collapse all data for squares
X_sq = [image_data(1:30).Xdata];
Y_sq = [image_data(1:30).Ydata];
RT_sq = [image_data(1:30).RT];

figure(1)
subplot(121);
scatter(1024-Y_sq,X_sq,[],RT_sq/1000,'filled');
title('Results for square images');
xlim([0 1024]);ylim([0 768]);
axis ij; axis equal;

subplot(122);
RT_sq_clip = RT_sq;RT_sq_clip(RT_sq_clip>3000) = 3000;
scatter(1024-Y_sq,X_sq,[],RT_sq_clip/1000,'filled');
title('Results for square images (RT clipped)');
xlim([0 1024]);ylim([0 768]);
axis ij; axis equal;

figure(2);
hist(RT_sq/1000,50)
xlabel('Reaction Times (sec)');
ylabel('Count');





%% Compare images of all black (red) with gray (pink) 

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

figure(2);
inx = (1:2)+2*(~first_tap);
if ~first_tap
    hold on;
end
bar(inx,[sum(a(1:5))./sum(K(1:5)) sum(b(1:5))./(sum(M(1:5))-sum(K(1:5))); sum(a(6:10))./sum(K(6:10)) sum(b(6:10))./(sum(M(6:10))-sum(K(6:10)))])
hold on;
[phat, pci] = binofit([sum(a(1:5)) sum(b(1:5)) sum(a(6:10)) sum(b(6:10))],[sum(K(1:5)) (sum(M(1:5))-sum(K(1:5))) sum(K(6:10)) (sum(M(6:10))-sum(K(6:10)))]);
err = sqrt(phat.*(1-phat))./sqrt([sum(K(1:5)) (sum(M(1:5))-sum(K(1:5))) sum(K(6:10)) (sum(M(6:10))-sum(K(6:10)))]);
x = [inx(1)-1/7 inx(1)+1/7 inx(2)-1/7 inx(2)+1/7];
errorbar(x,phat,err,'k','linestyle','none');
hold off;
ylabel('Tap rate')

%% Analyze square tap rates over tap number

% img_list = 1:5; %Black/gray comparison
% img_list = 21:25; %Red/pink comparison
img_list = [1:5 11:20 21:25];

offset = 5; %offset between spatially matching images (number of images in the group)

a = zeros(length(img_list),6);
b = zeros(length(img_list),6);
L = zeros(length(img_list),6);
K = zeros(length(img_list),6);

insquare_faint_all = [];
insquare_solid_all = [];
group_faint_all = cell(1,2);
group_solid_all = cell(1,2);
for k = 1:length(img_list)
    %read image
    im = imread(['../images/squares' num2str(img_list(k)) '.png']);
    im_in = double(sum(im,3));
    
    %create mask for faint square
    mask = im_in>255 & im_in<255*3;
    mask(1005:1024,:) = 0;
    
    %set up images with black border for display purposes
    pad = 10;
    im(1005:1024,:) = 255; 
    im = [ones(pad,size(im,2)+2*pad,3); ones(size(im,1),pad,3) im ones(size(im,1),pad,3); 
          ones(pad,size(im,2)+2*pad,3)]; %#ok<AGROW>
    im_offset = imread(['../images/squares' num2str(img_list(k)+offset) '.png']);
    im_offset(1005:1024,:) = 255;
    im_offset = [ones(pad,size(im_offset,2)+2*pad,3); 
                 ones(size(im_offset,1),pad,3) im_offset ones(size(im_offset,1),pad,3);
                 ones(pad,size(im_offset,2)+2*pad,3)]; %#ok<AGROW>   
    
    for tapnum = 2:2:12
        %Find taps in faint square
        sel_faint = image_data(img_list(k)).tapnum==tapnum;
        X = image_data(img_list(k)).Xdata(sel_faint);
        Y = image_data(img_list(k)).Ydata(sel_faint);
        insquare_faint = zeros(length(X),1);
        for l = 1:length(X)
            insquare_faint(l) = mask(Y(l),X(l));
        end

        %Find taps in corresponding non-faint square
        sel_solid = image_data(img_list(k)+offset).tapnum==tapnum;
        X = image_data(img_list(k)+offset).Xdata(sel_solid);
        Y = image_data(img_list(k)+offset).Ydata(sel_solid);
        insquare_solid = zeros(length(X),1);
        for l = 1:length(X)
            insquare_solid(l) = mask(Y(l),X(l));
        end

        %plot
        figure(1);
        subplot(121);

        imshow(imrotate(im,-90));
        hold all;
        scatter(1024-image_data(img_list(k)).Ydata(sel_faint)+pad,image_data(img_list(k)).Xdata(sel_faint)+pad,[],'g','filled');
        title([num2str(sum(insquare_faint)) ' out of ' num2str(sum(sel_faint))]);

        subplot(122);
        im_offset = imread(['../images/squares' num2str(img_list(k)+offset) '.png']);
        im_offset(1005:1024,:) = 255;
        im_offset = [ones(pad,size(im_offset,2)+2*pad,3); 
                     ones(size(im_offset,1),pad,3) im_offset ones(size(im_offset,1),pad,3);
                     ones(pad,size(im_offset,2)+2*pad,3)]; %#ok<AGROW>
        imshow(imrotate(im_offset,-90));
        hold all;
        scatter(1024-image_data(img_list(k)+offset).Ydata(sel_solid)+pad,image_data(img_list(k)+offset).Xdata(sel_solid)+pad,[],'g','filled');
        title([num2str(sum(insquare_solid)) ' out of ' num2str(sum(sel_solid))]);

        boldify
        
        a(k,tapnum/2) = sum(insquare_faint);
        b(k,tapnum/2) = sum(insquare_solid);
        K(k,tapnum/2) = sum(sel_faint);
        L(k,tapnum/2) = sum(sel_solid);
        
        rate(k,tapnum/2) = mean(insquare_faint);
        err(k,tapnum/2) = 1.96*std(insquare_faint)/sqrt(length(insquare_faint));
        
        hold off;
        
        insquare_faint_all = [insquare_faint_all; insquare_faint];
        insquare_solid_all = [insquare_solid_all; insquare_solid];
        
        group_faint_all{1} = [group_faint_all{1}; repmat(k,[length(insquare_faint) 1])]; 
        group_faint_all{2} = [group_faint_all{2}; repmat(tapnum/2,[length(insquare_faint) 1])]; 
        
        group_solid_all{1} = [group_solid_all{1}; repmat(k,[length(insquare_solid) 1])]; 
        group_solid_all{2} = [group_solid_all{2}; repmat(tapnum/2,[length(insquare_solid) 1])];

%         pause
    end
    
end

% Set up for the contingency table
%
%                       insquare   elsewhere  total
% first presentation    a1         -          K1
% other presentation    a2         -          K2
% total                 N          -          M

a1 = sum(a(:,1:2),2);
a2 = sum(a(:,5:6),2);
K1 = sum(K(:,1:2),2);
K2 = sum(K(:,5:6),2);
N = a1+a2;
M = K1+K2;

b1 = b(:,1);
b2 = sum(b(:,2:end),2);
L1 = L(:,1);
L2 = sum(L(:,2:end),2);
Nb = b1+b2;
Mb = L1+L2;

for k = 1:length(a1)
    p(k) = fexact( a1(k),M(k),K1(k),N(k));
    disp(['pfaint = ' num2str(p(k))]);
    pb(k) = fexact( b1(k),Mb(k),L1(k),Nb(k));
    disp(['psolid = ' num2str(pb(k))]);
end
% group_faint_all = group_faint_all(2);
% [pfaint, t1] = anovan(insquare_faint_all,group_faint_all(1),'display','off')
% [psolid, t2] = anovan(insquare_solid_all,group_solid_all,'display','off')

% [tapnum, imgnum] = meshgrid(1:6,1:10);
% imgnum_binary = repmat(eye(10),6,1);
% 
% [B, dev, stats] = mnrfit([tapnum(:) imgnum_binary],[a(:) K(:)-a(:)])
%This is wrong too? I think so.

%% Reaction times over tap numbers
mean_RT_tapnum = zeros(12,1);
err_RT_tapnum = zeros(12,1);
% kurtosis_RT_tapnum = zeros(6,1);
% skewness_RT_tapnum = zeros(6,1);
for k = 1:12
    sel = tapnum_all==k;
%     figure(1);
%     subplot(4,3,k);
%     hist(RT_all(sel)/1000,.05:.2:10);
%     xlim([0 10]);
    mean_RT_tapnum(k) = mean(RT_all(sel));
    err_RT_tapnum(k) = 1.96*std(RT_all(sel))/sqrt(sum(sel));
%     kurtosis_RT_tapnum(k) = kurtosis(RT_all(sel));
%     skewness_RT_tapnum(k) = skewness(RT_all(sel));

%     xlabel('Time (Sec)');
%     title(['Tap Number ' num2str(k)]);
end
boldify
figure(2);
subplot(121);

hist(RT_all/1000,linspace(0,5,50))
xlim([0 5]);
xlabel('RT (sec)');
ylabel('Count');

subplot(122);
errorbar(1:2:12,mean_RT_tapnum(1:2:12)/1000,err_RT_tapnum(1:2:12)/1000)
hold all;
errorbar(2:2:12,mean_RT_tapnum(2:2:12)/1000,err_RT_tapnum(2:2:12)/1000)
hold off
xlabel('Presentation Number');
ylabel('Mean RT (Sec)');
legend('Natural Scenes','Images of Squares')
ylim([1.2 2.6001])
boldify

% figure(3);
% plot([skewness_RT_tapnum kurtosis_RT_tapnum/10]);
% 
% xlabel('Tap Number');
% ylabel('RT Kurtosis');
%% Compute singleton tap rates

img_list_start = [1 11 16 21];
grouplen = 5;
d = zeros(30,1);
a = zeros(30,1);
K = zeros(30,1);
d_rank = zeros(30,1);

for start = 1:length(img_list_start)
    
    for pic = img_list_start(start):img_list_start(start)+grouplen-1
        
        % Read image
        im = imread(['../images/squares' num2str(pic) '.png']);
        im_in = double(sum(im,3));
        
        % Create mask for signleton square
        [freq]  = hist(im_in(:),3);
        [~, inx]= min(freq);
        bins = unique(im_in(:));
        mask = im_in==bins(inx);
        mask(1005:1024,:) = 0;
        [row_sqr, col_sqr] = find(mask);
                
        % Euclidean distance to center
        d(pic) = sqrt(mean(col_sqr-size(mask,2)/2)^2 + mean(row_sqr-size(mask,1)/2)^2);
        
        % Find distance from center in terms of rank order 
        % (e.g. second closest square to the center)
        %find square centers
        sqr_tmp = conv2(im_in,ones(120,120),'same');
        Xmax = conv2(sign(diff(sqr_tmp,1,2)),[1 -1])==2;
        Ymax = conv2(sign(diff(sqr_tmp,1,1)),[1; -1])==2;
        ctr_mask = Xmax & Ymax; %object centers
        ctr_mask(1005:1024,:) = 0;
        [ctr_row, ctr_col] = find(ctr_mask);
        [sing_row, sing_col] = find(ctr_mask & mask);
        
        d_eu = sqrt((ctr_row-size(mask,1)/2).^2+(ctr_col-size(mask,2)/2).^2);
        [d_sort, inx] = sort(d_eu);
        d_rank(pic) = find(ctr_row(inx)==sing_row & ctr_col(inx) == sing_col);
        
        
        %Find taps in faint square
        X = image_data(pic).Xdata;
        Y = image_data(pic).Ydata;
        insquare_faint = zeros(length(X),1);
        for l = 1:length(X)
            insquare_faint(l) = mask(Y(l),X(l));
        end
        
        a(pic) = sum(insquare_faint);
        K(pic) = length(image_data(pic).Xdata);
    end
end

[rate, conf_int] = binofit(a,K); %find rate and 95% confidence interval
L = rate-conf_int(:,1);
U = conf_int(:,2)-rate;

%% Compute tap rates for non-singletons in control images

Nobj = 10;
pic_list = [6:10 26:30];
d_null = zeros(length(pic_list),Nobj);
d_rank_null = repmat(1:Nobj,[length(pic_list) 1]);
rate_null = zeros(length(pic_list),Nobj);
a_null = zeros(size(rate_null));
K_null = zeros(size(rate_null));

for p = 1:length(pic_list)
    pic = pic_list(p);
    im = imread(['../images/squares' num2str(pic) '.png']);
    im_in = double(sum(im,3));
    
    %find square centers
    sqr_tmp = conv2(im_in,ones(120,120),'same');
    Xmax = conv2(sign(diff(sqr_tmp,1,2)),[1 -1])==2;
    Ymax = conv2(sign(diff(sqr_tmp,1,1)),[1; -1])==2;
    ctr_mask = Xmax & Ymax;
    [ctr_row, ctr_col] = find(ctr_mask);
    
    d_null(p,:) = sqrt((ctr_row-size(mask,1)/2).^2+(ctr_col-size(mask,2)/2).^2);
    [d_sort, inx]= sort(d_null(p,:));
    d_rank_null(p,inx) = d_rank_null(p,:); %this line reverses the sort matrix into a rank matrix
    
    tmp = accumarray([ctr_row ctr_col],d_rank_null(p,:),size(ctr_mask));
    tmp = conv2(tmp,ones(20,20),'same');
    
%     figure(1);
%     subplot(121);
%     imshow(im);
%     subplot(122);
%     imagesc(tmp);
%     
%     pause
    
    X = image_data(pic).Xdata;
    Y = image_data(pic).Ydata;

    for k = 1:length(ctr_row)
        mask = zeros(size(im_in));
        mask(ctr_row(k)-59:ctr_row(k)+60,ctr_col(k)-59:ctr_col(k)+60)=1;
        mask(1005:1024,:) = 0;
        a_null(p,k) = sum(mask(sub2ind(size(mask),Y,X)));
%         imagesc(mask);
%         hold on;scatter(X,Y,'g');hold off;
%         title(num2str(rate_null(p,k)));
        
%         pause
        
    end
    K_null(p,:) = length(X);
    rate_null(p,:) = a_null(p,:)./K_null(p,:);
    
end


%% Compute non-singleton mean tap rates in bins based on distance from center
binwidth = 100;
bins_start = 0:binwidth:600-binwidth;
bins_end = binwidth:binwidth:600;
rate_null_mean = zeros(length(bins_start),1);
for k = 1:length(bins_start)
    sel_faint = d_null>=bins_start(k) & d_null<bins_end(k);
    rate_null_mean(k) = mean(rate_null(sel_faint));
end
bins_ctr = (bins_start+bins_end)/2;

rate_null_mean_rank = zeros(max(d_rank_null(:)),1);
for k = 1:max(d_rank_null(:))
    sel_faint = d_rank_null == k;
    rate_null_mean_rank(k) = mean(rate_null(sel_faint));
end

% std_err = 1.96*sqrt(rate.*(1-rate)./K);

%% Perform Pairwise Linear model comparisons
d_expand = cell(size(d));
rate_expand = cell(size(rate));
for k = 1:length(d)
   d_expand{k} = repmat(d(k),K(k),1);
   rate_expand{k} = [ones(a(k),1); zeros(K(k)-a(k),1)];
end

d_null_expand = cell(size(d_null));
rate_null_expand = cell(size(rate_null));
for k = 1:length(d_null(:))
    d_null_expand{k} = repmat(d_null(k),K_null(k),1);
    rate_null_expand{k} = [ones(a_null(k),1); zeros(K_null(k)-a_null(k),1)];
end

x = cell(5,1);
x{5} = vertcat(d_null_expand{:});
x{1} = vertcat(d_expand{1:5}); %gray/black
x{2} = vertcat(d_expand{11:15}); %black/gray
x{3} = vertcat(d_expand{16:20}); %blue/yellow
x{4} = vertcat(d_expand{21:25}); %pink/red
y = cell(5,1);
y{5} = vertcat(rate_null_expand{:});
y{1} = vertcat(rate_expand{1:5});
y{2} = vertcat(rate_expand{11:15});
y{3} = vertcat(rate_expand{16:20});
y{4} = vertcat(rate_expand{21:25});

%params for full linear model
x_all = vertcat(x{:});
y_all = vertcat(y{:});
ind1 = [ones(length(x{1}),1); zeros(length(x{2}),1); zeros(length(x{3}),1); zeros(length(x{4}),1); zeros(length(x{5}),1)];
ind2 = [zeros(length(x{1}),1); ones(length(x{2}),1); zeros(length(x{3}),1); zeros(length(x{4}),1); zeros(length(x{5}),1)];
ind3 = [zeros(length(x{1}),1); zeros(length(x{2}),1); ones(length(x{3}),1); zeros(length(x{4}),1); zeros(length(x{5}),1)];
ind4 = [zeros(length(x{1}),1); zeros(length(x{2}),1); zeros(length(x{3}),1); ones(length(x{4}),1); zeros(length(x{5}),1)];
ind5 = [zeros(length(x{1}),1); zeros(length(x{2}),1); zeros(length(x{3}),1); zeros(length(x{4}),1); ones(length(x{5}),1)];
mdl_full = LinearModel.fit([x_all ind1 ind2 ind3 ind4], y_all,'interactions');



%loop through pairwise models
p_intercept = ones(5,5);
p_slope = ones(5,5);
% slope = zeros(5,5);

for k = 1:4
    for l = k+1:5
        x_all = [x{k}; x{l}];
        y_all = [y{k}; y{l}];
        ind = [zeros(length(x{k}),1); ones(length(x{l}),1)];
        mdl = LinearModel.fit([x_all ind], y_all,'interactions');
        p = mdl.Coefficients.pValue;
        
%         if p(4)>0.1
%             p_slope(k,l) = p(4);
%             mdl = LinearModel.fit([x_all ind], y_all);
%             p = mdl.Coefficients.pValue;
% %             coef = mdl.Coefficients.Estimate;
% %             slope(k,l) = coef(2);
%         end
        if length(p)>2
            p_intercept(k,l) = p(3);
            if length(p)>3
                p_slope(k,l) = p(4);
            end
        end
        
%         pause
    end
end


% thresh = 0.05;
% p_intercept<thresh
% a = sort(p_intercept(:));
% a = a(1:10);
% sig  = a<thresh./(10./(1:10)');
% sum(sig)

%individual models
mdl_ind = cell(5,1);
intercept = zeros(5,1);
intercept_SE = zeros(5,1);
for k = 1:5
    mdl_ind{k} = LinearModel.fit(x{k}, y{k});
    coef = mdl_ind{k}.Coefficients.Estimate;
    intercept(k) = coef(1);
    intercept_SE(k) = mdl_ind{k}.Coefficients.SE(1);
%     intercept_CI(k,:) = coef_CI(1,:);
end

%%
list = [2 3 4 1 0]; %used for x axis placement on bar graphs
figure(1);
for l = 1:5
    bar(list(l),intercept(l));
    hold on;
end
errorbar(list,intercept,intercept_SE,'o')
hold off;
ylabel('Intercept value');
xlim([0.5 4.5])

x_all = vertcat(x{:});
y_all = vertcat(y{:});
mdl_final = LinearModel.fit([x_all ind1 ind2 ind3 ind4 ind5 x_all.*(ind5)], y_all,'intercept',false);
coef = mdl_final.Coefficients.Estimate;
coef_SE = mdl_final.Coefficients.SE;


figure(2)
for l = 1:5
    bar(list(l),coef(l+1));
    hold on;
    
end
errorbar(list,coef(2:6),coef_SE(2:6,1),'o');
xlim([-0.5 4.5])
hold off;
legend('gray/black','black/gray','blue/yellow','pink/red','no singleton','error bars')
% ind4_1 = [zeros(length(x{1}),1); ones(length(x{4}),1)];
% mdl2 = LinearModel.fit([[x{1}; x{4}] ind4_1],[y1; y4])
% 
% ind3_2 = [zeros(length(x{2}),1); ones(length(x{3}),1)];
% mdl3 = LinearModel.fit([[x{2}; x{3}] ind3_2],[y2; y3])
% 
% mdl4 = LinearModel.fit([x ind2+ind3], y)




%% Plot tap rates
figure(2);
% errorbar(d(1:5),rate(1:5),L(1:5),U(1:5),'o')
scatter(d(1:5),rate(1:5),'filled')
hold all;


% errorbar(d(11:15),rate(11:15),L(11:15),U(11:15),'o')
% errorbar(d(16:20),rate(16:20),L(16:20),U(16:20),'o')
% errorbar(d(21:25),rate(21:25),L(21:25),U(21:25),'o')
scatter(d(11:15),rate(11:15),'filled')
scatter(d(16:20),rate(16:20),'filled')
scatter(d(21:25),rate(21:25),'filled')
scatter(d_null(:),rate_null(:));


hold off;
ylabel('Tap rate');
xlabel('Distance from center (pixels)');


figure(2)
% X_nosing = [image_data([6:10 26:30]).Xdata];
% Y_nosing = [image_data([6:10 26:30]).Ydata];
% d_all = sqrt((X_nosing-size(mask,2)/2).^2 + (Y_nosing-size(mask,1)/2).^2);

hold all;
% [p, d_bins] = hist(d_all,60:120:max(d_all));
% plot(d_bins,p/sum(p));
dtmp1 = d_null(1:5,:);dtmp2 = d_null(6:10,:);
rtmp1 = rate_null(1:5,:); rtmp2 = rate_null(6:10,:);
% scatter(dtmp1(:),rtmp1(:),'x');
% scatter(dtmp2(:),rtmp2(:),'x');



coef = mdl_final.Coefficients.Estimate;
b = coef(1)+[0 coef(3) coef(6) coef(4) coef(5)];
figure(1);
for l = 1:length(b)
    bar(l,b(l));
    hold all;
end
hold off;
ylim([0 1]);
legend('No Singleton','Gray/Black','Pink/Red','Black/Gray','Blue/Yellow');
ylabel('Intercept value');

figure(2);
hold on
w = [d_null(:) ones(size(d_null(:)))]\rate_null(:);
plot([0 500],w(1)*[0 500]+w(2),'-g')
w = [d(1:5) ones(5,1)]\rate(1:5);
plot([0 500],w(1)*[0 500]+w(2),'-k')
w = [d(11:15) ones(5,1)]\rate(11:15);
plot([0 500],w(1)*[0 500]+w(2),'-k')
w = [d(16:20) ones(5,1)]\rate(16:20);
plot([0 500],w(1)*[0 500]+w(2),'-b')
w = [d(21:25) ones(5,1)]\rate(21:25);
plot([0 500],w(1)*[0 500]+w(2),'-r')


hold off
legend('Gray/Black','Black/Gray','Blue/Yellow','Pink/Red','Non Singletons');
ylim([0 1]);
xlim([0 500]);
boldify

figure(4);
scatter(d(1:5),rate(1:5),'filled')
hold all;
scatter(d(11:15),rate(11:15),'filled')
scatter(d(16:20),rate(16:20),'filled')
scatter(d(21:25),rate(21:25),'filled')
scatter(d_null(:),rate_null(:));


w = [coef(1)+coef(7) coef(6)];
plot([0 500],w(1)*[0 500]+w(2),'-g')
w = [coef(1) coef(2)];
plot([0 500],w(1)*[0 500]+w(2),'-k')
w = [coef(1) coef(3)];
plot([0 500],w(1)*[0 500]+w(2),'-k')
w = [coef(1) coef(4)];
plot([0 500],w(1)*[0 500]+w(2),'-b')
w = [coef(1) coef(5)];
plot([0 500],w(1)*[0 500]+w(2),'-r')
hold off
legend('Gray/Black','Black/Gray','Blue/Yellow','Pink/Red','No Singleton Fit');
ylim([0 1]);
ylabel('Singleton Tap rate');
xlabel('Distance from center (pixels)');


figure(3);
errorbar(d_rank(1:5),rate(1:5),L(1:5),U(1:5),'o')
% scatter(d_rank(1:5),rate(1:5),'filled')
hold all;


errorbar(d_rank(11:15),rate(11:15),L(11:15),U(11:15),'o')
errorbar(d_rank(16:20),rate(16:20),L(16:20),U(16:20),'o')
errorbar(d_rank(21:25),rate(21:25),L(21:25),U(21:25),'o')
% scatter(d_rank(11:15),rate(11:15),'filled')
% scatter(d_rank(16:20),rate(16:20),'filled')
% scatter(d_rank(21:25),rate(21:25),'filled')


hold off;
ylabel('Tap rate');
xlabel('Distance from center (Rank)');


figure(3)
% X_nosing = [image_data([6:10 26:30]).Xdata];
% Y_nosing = [image_data([6:10 26:30]).Ydata];
% d_all = sqrt((X_nosing-size(mask,2)/2).^2 + (Y_nosing-size(mask,1)/2).^2);

hold all;
% [p, d_bins] = hist(d_all,60:120:max(d_all));
% plot(d_bins,p/sum(p));
dtmp1 = d_rank_null(1:5,:);dtmp2 = d_rank_null(6:10,:);
rtmp1 = rate_null(1:5,:); rtmp2 = rate_null(6:10,:);
bar(1:max(d_rank_null(:)),rate_null_mean_rank);
scatter(dtmp1(:),rtmp1(:),'x');
scatter(dtmp2(:),rtmp2(:),'x');


% w = [d(1:5) ones(5,1)]\rate(1:5);
% plot([0 500],w(1)*[0 500]+w(2),'-k')
% w = [d(11:15) ones(5,1)]\rate(11:15);
% plot([0 500],w(1)*[0 500]+w(2),'-k')
% w = [d(16:20) ones(5,1)]\rate(16:20);
% plot([0 500],w(1)*[0 500]+w(2),'-b')
% w = [d(21:25) ones(5,1)]\rate(21:25);
% plot([0 500],w(1)*[0 500]+w(2),'-r')


hold off
legend('Gray singleton with Black distractors','Black/Gray','Blue/Yellow','Pink/Red','No singleton avg');
boldify


