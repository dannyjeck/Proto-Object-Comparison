% Data format is: subj number, gender, age, image number, Y, X, reaction time, stop time (relative to new years 1970 in ms)
% Each new line is a new tap.

% Matlab file reading causes the new lines to run together, so the subject
% numbers and stop times in the subj variable below are formatted wierdly.

dir = '../data_prelim2/';

%% get all the data
Nsubj = 138;
subj = cell(Nsubj,1);
for k = 1:Nsubj
    if k>=100
        fname = [dir 'subject_0' num2str(k) '.txt'];
    elseif k>=10
        fname = [dir 'subject_00' num2str(k) '.txt'];
    else
        fname = [dir 'subject_000' num2str(k) '.txt'];
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
image_data = repmat(struct('Xdata',[],'Ydata',[],'RT',[]),80,1);

for k = 1:Nsubj

    if first_tap;
        m = 2;
    else
        m = size(subj{k},2);
    end
    
    for l = 1:m
        im_idx = str2double(subj{k}{3,l});
        image_data(im_idx).Ydata(end+1) = str2double(subj{k}{4,l});
        image_data(im_idx).Xdata(end+1) = str2double(subj{k}{5,l});
        image_data(im_idx).RT(end+1)    = str2double(subj{k}{6,l});
    end
end

% Data for all images together
X_all = [image_data(:).Xdata];
Y_all = [image_data(:).Ydata];
RT_all = [image_data(:).RT];


%% Analyze taps by subject over images


%% Analyze taps independent of image or subject

scatter(1024-Y_all,X_all,[],RT_all,'filled');
xlim([0 1024]);ylim([0 768]);
axis ij; axis equal;

%% Show all images with taps overlaid

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
        scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],'MarkerEdgeColor','g','MarkerFaceColor','r');
    else
        scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],'MarkerEdgeColor','g','MarkerFaceColor','r');
    end
    title(num2str(k))
    hold off;
    
%     if first_tap
%         saveas(f,['./images_138subj/firsttap_' num2str(img_list(k)) '.png']);
%     else
%         saveas(f,['./images_138subj/alltaps_' num2str(img_list(k)) '.png']);
%     end

	pause
end

%% Show taps and interest data 
% Requires getting the interest map data "int_data_norm" first

img_list = 31:80;
list = [1:25 76:100];

for k = 1:length(img_list)
figure(1)
subplot(121)
im = imread(['../images/squares' num2str(img_list(k)) '.png']);
    im = permute(im,[2 1 3]);
    im = im(:,end:-1:1,:);
    imshow(im);
    hold all;
    scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],'g');
    title(num2str(k))
    hold off;
subplot(122)
imagesc(int_data_norm{list(k)});colorbar; colormap(gray);axis equal;
    hold all;
    scatter(1024-image_data(img_list(k)).Ydata,image_data(img_list(k)).Xdata,[],'g');
    hold off;
pause
end
%% Compare images of all black (red) with gray (pink) 
img_list = 1:5; %Black/gray comparison
% img_list = 21:25; %Red/pink comparison
offset = 5;

for k = 1:length(img_list)
    im = imread(['../images/squares' num2str(img_list(k)) '.png']);
    im_in = double(sum(im,3));
    
    mask = im_in>255 & im_in<255*3;
    mask(1005:1024,:) = 0;
    
    
    X = image_data(img_list(k)).Xdata;
    Y = image_data(img_list(k)).Ydata;
    insquare_faint = zeros(length(X),1);
    for l = 1:length(X)
        insquare_faint(l) = mask(Y(l),X(l));
    end
    
    X = image_data(img_list(k)+offset).Xdata;
    Y = image_data(img_list(k)+offset).Ydata;
    insquare_solid = zeros(length(X),1);
    for l = 1:length(X)
        insquare_solid(l) = mask(Y(l),X(l));
    end

    
    figure(1);
    subplot(121);
    hold all;
    imshow(['../images/squares' num2str(img_list(k)) '.png']);
    scatter(image_data(img_list(k)).Xdata,image_data(img_list(k)).Ydata,[],'g');
    title([num2str(sum(insquare_faint)) ' out of ' num2str(length(image_data(img_list(k)).Xdata))]);
    
    subplot(122);
    hold all;
    imshow(['../images/squares' num2str(img_list(k)+offset) '.png']);
    scatter(image_data(img_list(k)+offset).Xdata,image_data(img_list(k)+offset).Ydata,[],'g');
    title([num2str(sum(insquare_solid)) ' out of ' num2str(length(image_data(img_list(k)+offset).Xdata))]);
    
    %        R       NR
    %  pos   a        c       K
    %  neg   b        d       -
    % total  N        -       M
    a(k) = sum(insquare_faint);
    M(k) = length(image_data(img_list(k)).Xdata)+length(image_data(img_list(k)+offset).Xdata);
    K(k) = length(image_data(img_list(k)).Xdata);
    N(k) = sum(insquare_faint)+sum(insquare_solid);
    p(k) = fexact( a(k),M(k),K(k),N(k), 'tail','r');
    disp(['p = ' num2str(p(k))]);
    pause;
    
end
ptot = fexact( sum(a),sum(M),sum(K),sum(N), 'tail','r');
disp(['p = ' num2str(ptot)]);
