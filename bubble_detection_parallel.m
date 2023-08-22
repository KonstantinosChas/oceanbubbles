%--------------------------------------------------------------------------
% 2D-image processing algorithm : 
% cluster version for HPC-UCL - Linux version
% Project: bubbles of breaking waves
% part 0: specify number of loaded files and pre-allocate 
% part 1: image pre-processing - Segmentation
% part 2: 3-step modified Hough Transform and Snakes (active contours);
% part 4 : analysis of results: bubble size distributions (load data to
% other script
%--------------------------------------------------------------------------
clc;
clear;
close;
warning off;
tic;
%--------------------------------------------------------------------------
% profile on -history;
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE THESE VALUES:
% name, NUM_in,NUM and file parameters should be adjusted 
T = 0;   % (0) windows or (1) linux
disp('did you change the parameters?');
disp('which phase? (peak,trough,posp2,minp2)');
% phase = input('','s');
phase = 'minp2';
name =['GW_175_',phase];   % name structure of images. images in folder are numbered after name (e.g. ..trough_0001)
name_0 = ['GW_175_',phase,'_000'];
name_10 = ['GW_175_',phase,'_00'];
name_100 = ['GW_175_',phase,'_0'];
name_1000 = ['GW_175_',phase,'_'];
% when running with no display
set(0,'DefaultFigureVisible','off'); % no graphics/plots
% choose which images are going to be loaded and analysed 
NUM_in = 1;   % choose first image index (check file names in folder)
NUM = 2;    %   choose last image index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
date = '191020';
wave = 'GW_175';
ROI = 'ROI1';
run = 'run1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('initial image index?');
% NUM_in =  input(''); 
% disp('last image index?');
% NUM =  input(''); 
% disp('number of peak images?');
% NUM_peak =  input(''); % number of iterations (images);
% disp('number of trough images?');
% NUM_trough = input('');
% disp('number of +p2 images?');
% NUM_posp2 = input('');
% disp('number of -p2 images?');
% NUM_minp2 = input('');
% NUM = NUM_minp2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% specify directories and paths 
if T ==0
    %%%windows 
    dir = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\images\',wave,'\',...
    date,'\',ROI,'\',run,'\'];  % directory of images to analyze
    dir_script = 'C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\matlab\';
    dir_target = 'C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\images\calibration\bubbles'; % directory of targets to use for calibration 
elseif T==1
    %%%linux
    dir = ['/mnt/DATA/users/konstantinos/',date,'/',wave,'/',stage,'/',ROI,'/',run,'/'];  % directory of images to analyze
    dir_script = '/home/konstantinos/apps/MATLAB/scripts/bubbles';
    dir_target = '/mnt/DATA/users/konstantinos/cam_calibration/bubbles/'; % directory of targets to use for calibration 
end
addpath(dir);
addpath(dir_script);
cd(dir_script);
cd(dir);
% allocate parameters + variables
% parameters for image processing 
scaling_factor = 1./4.1;   % scaling image factor - bubble resolution
imtype = 'jpg';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load images 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_name = cell(1,NUM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Load images and pre-processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor i = NUM_in : NUM % parallel computing
% for i = NUM_in:NUM
    if i>=1000
        im_name{i} = ([name_1000,num2str(i)]);
    elseif i>=100 && i<1000
        im_name{i} = ([name_100,num2str(i)]);
    elseif i<100 && i>=10
        im_name{i} = ([name_10,num2str(i)]);
    elseif i<10
        im_name{i} = ([name_0,num2str(i)]);
     end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor i = NUM_in : NUM % parallel computing
% for i = NUM_in : NUM
%     clearvars -except im_name NUM im_type
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     pack;
%     clear images bin_im h fd L wr rm im
%     clear fim em LIm g2 L2 f2
%     clear folderName
    images = struct();
    bin_im = struct();
    h = struct();
    fd = struct();
    L = struct();
    wr = struct();
    rm = struct();
    im = struct();
    fim = struct();
    em = struct();
    LIm = struct();
    g2 = struct();
    L2 = struct();
    f2 = struct();

    folderName        =   imread([im_name{i},'.',imtype]);
%     folderfieldname{i} =   im_name{i};
    images.(im_name{i}) = folderName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% image calibration 
%     figure(1); imshow(rgb);
    images.(im_name{i}) = camcal_server(20,[im_name{i},'.',imtype],dir,dir_target); % input: [size of chessboard square in mm, image]
%     disp(['image calibrated ',num2str(i)]);
%% pre-processing : convert to grayscale - binary - threshold - filters - 
    bin_im.(im_name{i}) = im2bw(images.(im_name{i}),graythresh(images.(im_name{i})));
%     figure(2); imshow(bin_im.(folderfieldname{i}));
%     [L1.(im_name{i}),num.(im_name{i})] = bwlabel(bin_im.(im_name{i}),4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.(im_name{i}) = fspecial('sobel'); % motion-blur filter 
    fd.(im_name{i}) = im2double(images.(im_name{i}));  % convert to double
    images.(im_name{i}) = sqrt(imfilter(fd.(im_name{i}),h.(im_name{i}),'replicate').^2 ...
        + imfilter(fd.(im_name{i}),h.(im_name{i})','replicate').^2);  % filter, specifying the replicate boundary option
    L.(im_name{i}) = watershed(images.(im_name{i})); % apply watershed transform
    wr.(im_name{i}) = L.(im_name{i}) == 0;
    L.(im_name{i}) = L.(im_name{i});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WATERSHED TRANSFORM AND SEGMENTATION
    rm.(im_name{i}) = imregionalmin(images.(im_name{i}));  % find local minima 
    im.(im_name{i}) = imextendedmin(images.(im_name{i}),2); % extended maxima transform 
    fim.(im_name{i}) = images.(im_name{i});
    fim.(im_name{i})(im.(im_name{i})(im.(im_name{i}))) = 175;
    fim.(im_name{i}) = im.(im_name{i});
    LIm.(im_name{i}) = watershed(bwdist(im.(im_name{i}))); % watershed of the distance transform of the internal marker image
    em.(im_name{i}) = LIm.(im_name{i}) == 0;
    LIm.(im_name{i}) = LIm.(im_name{i});
    g2.(im_name{i}) = imimposemin(images.(im_name{i}),im.(im_name{i}) | em.(im_name{i}));  % regional minima only occur in marked locations
    L2.(im_name{i})= watershed(g2.(im_name{i}));
    f2.(im_name{i}) = images.(im_name{i});
    f2.(im_name{i})(L2.(im_name{i}) ==0) = 255;
    L2.(im_name{i}) = L2.(im_name{i});
    f2.(im_name{i}) = f2.(im_name{i}); 
%     disp(['image segmented ',num2str(i)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOUGH TRANSFORM - detect circular objects, 
% i.e. bubbles in pre-processed images. 
% 3-stage Modified Hough Transform: function hough_detect
% stage 1: identify circles in segmented images allowing for large errors (high sensitivity)
% stage 2: use results from stage 1 as masks. then use active contours (Snakes) method 
% stage 3: perform final Hough transform with different settings on Snake processed images. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  best settings for active AND quiescent stage: 
%  edgethreshold: 0.40, 0.20, 0.20, 0.20 
%  sensitivity: 0.9,0.95,0.95,0.90
    rmin_xs = 1; % in pixels
    rmax_xs= 30;
    centers = [];
    radii = [];
    [centers.(im_name{i}), radii.(im_name{i})] = hough_detect(rmin_xs,rmax_xs,f2.(im_name{i}));

% REMEMBER TO CHANGE RMIN-RMAX IN ANALYSIS!!
% Save files in allocated directory.
    if T==0
    %%%windows
        dir_save = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\'...
        ,wave,'\',date,'\',ROI,'\',run];
    elseif T==1
    %%%linux
        dir_save =  ['/mnt/DATA/users/konstantinos/output/',date,'/',wave,'/',stage,'/',ROI,'/',run,'/']; 
    end
    cd(dir_save);
    filecodex = '-v7.3';
    savetofile(radii,filecodex,['radii.',im_name{i},'.mat'])
    savetofile(centers,filecodex,['centers.',im_name{i},'.mat'])
%     save(['radii.',im_name{i},'.mat'],'radii','-v7.3');
    disp(['bubble detection done ',num2str(i)]);
%     [M,X,C] = inmem;

end
toc;
% save(['radii_tot_',name,'.mat'],'-struct','radii','-v7.3');
% save(['centers_tot_',name,'.mat'],'-struct','centers','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = profile('info');
% save myprofiledata p;
% fin_time = toc;
disp('all images analyzed, goodbye!'); 
cd(dir);
% exit;


  
  