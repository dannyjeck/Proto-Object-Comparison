% After gathering data from 138 subjects, I realized that I had an uneven
% number of first taps between images I wanted to compare. This makes it
% more difficult to aggregate data across images to compare the groups
% overall.
%
% Annnd after gathering data from 196 subjects, I realized I had a sign
% wrong and need to double-correct now.
%
% By Danny Jeck

rng default 

faint_nsubj = [11 7 14 8 3 2 15 3 7 8];
solid_nsubj = [5 11 4 4 7 10 3  7 5 6];

offset = 5;

im_list = [1:5  21:25];
im_list((faint_nsubj-solid_nsubj)>0) = im_list((faint_nsubj-solid_nsubj)>0)+offset;
nsubj = abs(faint_nsubj-solid_nsubj);

list_out = [];
for k = 1:length(im_list(:));
    list_out(end+1:end+nsubj(k),1) = im_list(k)*ones(nsubj(k),1);
end

list_out = list_out(randperm(length(list_out)));

%I put list_out in the javascript for the test to correct for the
%imbalance.