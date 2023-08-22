%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------BUBBLE SPACE/TIME ANALYSIS------------------------
% Bubble Size Distributions and bubble statistics
% REGIONS OF INTEREST (ROIold) ARE SUBDIVIDED INTO SMALLER ROI
% NEW ROI ARE ROIold/4 in x direction AND ARE RENAMED AND REPOSITIONED IN NEW 
% UNIFIED X-Y COORDINATE SYSTEM
% ROI1-ROI9 ARE NEW ROI
% LOAD ROI1-ROI6 AND SORT/RESHAPE/RENAME
% Compute BSD, VBSD, void fraction 
%% 1. time averaged statistics and 
%% 1a. instanteneous statistics
%% 2. ensemble averages of many experimental runs
% advice : do not to mix phase shifts in same Matlab execution (check NUM_case for execution sequence)
%% ----------------------BUBBLE SPACE/TIME ANALYSIS------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close;
T = 0;   % (0) windows or (1) linux
set(0,'DefaultFigureVisible','off'); % no graphics/plots
% disp('how many repetitions or cases?');
% choose what image data are going to be loaded and analysed 
NUM_in = 1;   % first image (check file names in folder)
NUM = 9999;   %   last image 
% NUM_case = input('');
NUM_case = 18;  % choose number of cases (repetitions, phase shifts, amps, spectra)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 0;  % choose whether to load results from previous runs L=0(NO)/1(YES)
avg = 1;  % option to find averages of many wave repetitions avg= 0(NO)/1(YES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ALWAYS CHECK: DATE,WAVE,NUMBER OF IMAGES,TIME INDICES IN SCRIPT');
date = '191020';
wave = 'GW_18';
stage = 'entire_video';    % maybe create txt parameter file with .exe
stage_save = 'entire_video';
rep = 'wavereps';
comments = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BSD parameters
image_length = 800; % in px
image_width = 1280; % in px
scaling_factor = 1./4.1; % pixel to mm conversion
image_length_mm = image_length.*scaling_factor; % in mm
image_width_mm = image_width.*scaling_factor; % in mm
per_bin = (2.4*10.^-1)/1;   % norm factor radii bin in mm
% per_bin = 1;  % norm factor not used
vol_ROI = 350.*50.*400;  % volume of each ROI in mm^3;
% vol_ROI_mm = vol_ROI.*10^9; % volume of each ROI in mm^3
% vol_cube = 450.*200.*400; % vol factor not used
vol_cube = 10^9;    % vol of 1m wide cube in mm^3
vol_factor = vol_cube./vol_ROI;
norm_bin_size = 1./per_bin;   % norm factor bin size
rmin = 1;               % min radius in px
rmax = 30;              % max radius in px
% alocate global variables
clear im_name
im_name = cell(1,NUM-NUM_in+1);
mean_radius_mm = zeros(1,NUM-NUM_in+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------IMPORTANT-----------------------------------
ROI0_idx = 0;ROI1_idx = 0; ROI2_idx = 0; ROI3_idx = 0; ROI4_idx = 0; ROI5_idx = 0; ROI6_idx = 0;
ROI0_mtx_idx = 0;ROI1_mtx_idx = 0; ROI2_mtx_idx = 0; ROI3_mtx_idx = 0; ROI4_mtx_idx = 0; 
ROI5_mtx_idx = 0;ROI6_mtx_idx = 0; ROI7_mtx_idx = 0; ROI8_mtx_idx = 0; ROI9_mtx_idx = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BSD_avg_ROI0 = struct();
BSD_avg_ROI1 = struct();
BSD_avg_ROI2 = struct();
BSD_avg_ROI3 = struct();
BSD_avg_ROI4 = struct();
BSD_avg_ROI5 = struct();
BSD_avg_ROI6 = struct();
BSD_avg_ROI7 = struct();
BSD_avg_ROI8 = struct();
BSD_avg_ROI9 = struct();
%% ------------------------------------------------------------------------ 
phase = cell(1,NUM_case);
%% ------------------------------------------------------------------------ 
%% ------------------------------------------------------------------------ 
phase(1:end) = {'minp2'};
%% ------------------------------------------------------------------------ 
%% ------------------------------------------------------------------------ 
ROI = {'ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI1','ROI2','ROI3','ROI4','ROI5','ROI6'};
run = {'run1','run1','run1','run1','run1','run1','run2','run2','run2','run2','run2','run2','run3','run3','run3','run3','run3','run3'};
% ROI = {'ROI0','ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI0','ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI0','ROI1','ROI2','ROI3','ROI4','ROI5','ROI6'};
% run = {'run1','run1','run1','run1','run1','run1','run1','run2','run2','run2','run2','run2','run2','run2','run3','run3','run3','run3','run3','run3','run3'};
% run = {'run1','run2','run3'};
%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if L (LOAD) = 0 BSD is computed 
if L == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for case_idx = 1:NUM_case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     disp('you have to specify run, ROI and phase shift');
    %     disp('which run? (e.g. run1,..');
    % %     run{case_idx} = 'run1';
    %     run{case_idx} = input('','s');
    %     disp('which ROI? (e.g. ROI1-ROI5)');
    %     ROI{case_idx} = input('','s');
    % %     disp('which case? (e.g. peak,trough,posp2,minp2)');
    %     phase{case_idx} = 'minp2';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear centers radii
        if T ==0
    %WINDOWS 
            dir = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\',wave,'\',...
            date,'\server\',stage,'\',ROI{case_idx},'\',run{case_idx},'\'];  % directory of images to analyze
            dir_script = 'C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\matlab\server';
            dir_save = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\'...
            ,wave,'\',date,'\server\',stage_save,'\ROI1-8\'];
            dir_save_rep = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\'...
            ,wave,'\',date,'\server\',stage_save,'\ROI1-8\',rep,'\'];
        elseif T==1
            % Linux directories
            
        end
        addpath(dir);
        addpath(dir_script);
        cd(dir_script);
        cd(dir);
        addpath(dir_save);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clear variables 
    %% FILE NAMES BASED ON IMAGE NAMES
    %       phase = {'peak'};
        name =[wave,'_',phase{case_idx}];   % names of images (e.g. ...peak_0001)
        name_0 = [wave,'_',phase{case_idx},'_000'];
        name_10 = [wave,'_',phase{case_idx},'_00'];
        name_100 = [wave,'_',phase{case_idx},'_0'];
        name_1000 = [wave,'_',phase{case_idx},'_'];
        for i = NUM_in:NUM
                if i>=1000
                    im_name{i-NUM_in+1} = ([name_1000,num2str(i)]);
                elseif i>=100 && i<1000
                    im_name{i-NUM_in+1} = ([name_100,num2str(i)]);
                elseif i<100 && i>=10
                    im_name{i-NUM_in+1} = ([name_10,num2str(i)]);
                elseif i<10
                    im_name{i-NUM_in+1} = ([name_0,num2str(i)]);
                end
%% this is only due to a mistake in file names in GW_18, it will be changed in future and be obsolete 
        wave_aux = 'GW_175';
        name_aux =[wave_aux,'_',phase{case_idx}];   
        name_aux_0 = [wave_aux,'_',phase{case_idx},'_000'];
        name_aux_10 = [wave_aux,'_',phase{case_idx},'_00'];
        name_aux_100 = [wave_aux,'_',phase{case_idx},'_0'];
        name_aux_1000 = [wave_aux,'_',phase{case_idx},'_'];
%% this is only due to a mistake in file names, it will be changed in future and be obsolete 
%% (BECAUSE INPUT_NAME#OUTPUT_NAME IN BUBBLE CODE)
                if i>=1000
                    im_name_aux{i-NUM_in+1} = ([name_aux_1000,num2str(i)]);
                elseif i>=100 && i<1000
                    im_name_aux{i-NUM_in+1} = ([name_aux_100,num2str(i)]);
                elseif i<100 && i>=10
                    im_name_aux{i-NUM_in+1} = ([name_aux_10,num2str(i)]);
                elseif i<10
                    im_name_aux{i-NUM_in+1} = ([name_aux_0,num2str(i)]);
                end
%% this is only due to a mistake in file names, it will be changed in future and be obsolete  

        end
        for i = NUM_in:NUM
            radii_in(i-NUM_in+1) = load(['radii.',im_name{i-NUM_in+1},'.mat']); %#ok<SAGROW>
            centers_in(i-NUM_in+1) = load(['centers.',im_name{i-NUM_in+1},'.mat']); %#ok<SAGROW> %  this is auxilliary
            radii.(im_name{i-NUM_in+1}) = radii_in(i-NUM_in+1).data.(im_name_aux{i-NUM_in+1});  % rename struct after loading
            centers_a.(im_name{i-NUM_in+1}) = centers_in(i-NUM_in+1).data.(im_name_aux{i-NUM_in+1}); %% this is auxilliary
            centers.(im_name{i-NUM_in+1}) = centers_in(i-NUM_in+1).data.(im_name_aux{i-NUM_in+1}); % the same
        end

        clear radii_in centers_in
        clear vol_b vol_b_xy radii_mm_xy centers_xy radii_mm_xy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % allocate BSD parameters
        for i = NUM_in:NUM
    % find max, min and means
            min_radius.(im_name{i-NUM_in+1}) = min(radii.(im_name{i-NUM_in+1}));
            max_radius.(im_name{i-NUM_in+1}) = max(radii.(im_name{i-NUM_in+1}));
            mean_radius.(im_name{i-NUM_in+1}) = mean(radii.(im_name{i-NUM_in+1}));
            tot_min_radius.(im_name{i-NUM_in+1}) = min(min_radius.(im_name{i-NUM_in+1}));
            tot_max_radius.(im_name{i-NUM_in+1}) = min(max_radius.(im_name{i-NUM_in+1}));
            tot_mean_radius.(im_name{i-NUM_in+1}) = mean(mean_radius.(im_name{i-NUM_in+1}));
    %--------------------------------------------------------------------------
    %% scale pixels to mm 
    %% choose precision:
            radii.(im_name{i-NUM_in+1}) = round(radii.(im_name{i-NUM_in+1}),2);
            radii_mm.(im_name{i-NUM_in+1}) = radii.(im_name{i-NUM_in+1}).*scaling_factor;
            radii_mm.(im_name{i-NUM_in+1}) = round(radii_mm.(im_name{i-NUM_in+1}),3);
    %--------------------------------------------------------------------------        
            mean_radius_mm(i-NUM_in+1) = mean_radius.(im_name{i-NUM_in+1}).*scaling_factor;
            max_radius_mm.(im_name{i-NUM_in+1}) = max_radius.(im_name{i-NUM_in+1}).*scaling_factor;
            min_radius_mm.(im_name{i-NUM_in+1}) = min_radius.(im_name{i-NUM_in+1}).*scaling_factor;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [radii_mm.(im_name{i-NUM_in+1}),sort_idx_radii.(im_name{i-NUM_in+1})] = sort(radii_mm.(im_name{i-NUM_in+1})); % sort radii in ascending order
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unified coordinate system for all ROI in two steps:
    % 1. local (0,0) for each ROI should be at bottom left corner of image for convenience
    % 2. global (0,0) is then adjusted to be the bottom left of ROI1_old
    % 3. global (0,0) is at oldROI1 x,y origin which corresponds in x=0 focus location in flume 
        for i = NUM_in:NUM
            centers_b.(im_name{i-NUM_in+1}) = centers.(im_name{i-NUM_in+1}); % auxilliary centers output
            if strcmp(ROI{case_idx},'ROI0')
                centers.(im_name{i-NUM_in+1})(:,2) = abs(centers.(im_name{i-NUM_in+1})(:,2) - image_length) - image_length./4; 
                centers_b.(im_name{i-NUM_in+1})(:,2) = abs(centers_b.(im_name{i-NUM_in+1})(:,2) - image_length) - image_length./4;                 
            elseif strcmp(ROI{case_idx},'ROI1')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% XY IN CENTERS (1280x800 pixels) IS YX IN FLUME!!!
                centers.(im_name{i-NUM_in+1})(:,2) = abs(centers.(im_name{i-NUM_in+1})(:,2) - image_length) + 0; % in pixels
                centers_b.(im_name{i-NUM_in+1})(:,2) = abs(centers_b.(im_name{i-NUM_in+1})(:,2) - image_length) + 0; % in pixels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% these are still the old ROI:
            elseif strcmp(ROI{case_idx},'ROI2')
                centers.(im_name{i-NUM_in+1})(:,2) = abs(centers.(im_name{i-NUM_in+1})(:,2) - image_length) + image_length./4; 
                centers_b.(im_name{i-NUM_in+1})(:,2) = abs(centers_b.(im_name{i-NUM_in+1})(:,2) - image_length) + image_length./4; 

            elseif strcmp(ROI{case_idx},'ROI3')
                centers.(im_name{i-NUM_in+1})(:,2) = abs(centers.(im_name{i-NUM_in+1})(:,2) - image_length) + 2.*image_length./4; 
                centers_b.(im_name{i-NUM_in+1})(:,2) = abs(centers_b.(im_name{i-NUM_in+1})(:,2) - image_length) + 2.*image_length./4; 

            elseif strcmp(ROI{case_idx},'ROI4')
                centers.(im_name{i-NUM_in+1})(:,2) = abs(centers.(im_name{i-NUM_in+1})(:,2) - image_length) + 3.*image_length./4; 
                centers_b.(im_name{i-NUM_in+1})(:,2) = abs(centers_b.(im_name{i-NUM_in+1})(:,2) - image_length) + 3.*image_length./4; 

            elseif strcmp(ROI{case_idx},'ROI5')
                centers.(im_name{i-NUM_in+1})(:,2) = abs(centers.(im_name{i-NUM_in+1})(:,2) - image_length) + 4.*image_length./4; 
                centers_b.(im_name{i-NUM_in+1})(:,2) = abs(centers_b.(im_name{i-NUM_in+1})(:,2) - image_length) + 4.*image_length./4; 
            elseif strcmp(ROI{case_idx},'ROI6')
                centers.(im_name{i-NUM_in+1})(:,2) = abs(centers.(im_name{i-NUM_in+1})(:,2) - image_length) + 5.*image_length./4; 
                centers_b.(im_name{i-NUM_in+1})(:,2) = abs(centers_b.(im_name{i-NUM_in+1})(:,2) - image_length) + 5.*image_length./4; 
            end
        end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     % sort
     %% XY IN CENTERS IS YX IN FLUME!!!
        for i = NUM_in:NUM
            centers_mm.(im_name{i-NUM_in+1}) = centers.(im_name{i-NUM_in+1}).*scaling_factor; % scale in mm
            centers_mm.(im_name{i-NUM_in+1}) = round(centers_mm.(im_name{i-NUM_in+1}),1);  % significant figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [centers_xy.(im_name{i-NUM_in+1})(:,2),sort_idx_centers.(im_name{i-NUM_in+1})] = sort(centers_mm.(im_name{i-NUM_in+1})(:,2)); % sort centers in increasing x-flume-coordinate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            centers_xy.(im_name{i-NUM_in+1})(:,1) = centers_mm.(im_name{i-NUM_in+1})(sort_idx_centers.(im_name{i-NUM_in+1})(:,1)); % same for corresponding y-flume values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            centers.(im_name{i-NUM_in+1})(:,1) = centers_mm.(im_name{i-NUM_in+1})(sort_idx_radii.(im_name{i-NUM_in+1}),1); % sort centers according to sorted radii y coord
            centers.(im_name{i-NUM_in+1})(:,2) = centers_mm.(im_name{i-NUM_in+1})(sort_idx_radii.(im_name{i-NUM_in+1}),2); % sort centers according to sorted radii x coord
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% distributions:
%% create bins and bin array with equal bin sizes (adjust rad_bin_size)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------------------------------------------------------
        rad_bin_size = round(scaling_factor/1,2); % bin size in mm
        for i = NUM_in:NUM
            rmin_mm = round(rmin.*scaling_factor,1);
            rmax_mm = round(rmax.*scaling_factor,1);
%--------------------------------------------------------------------------
            nbin_radii_hg = (rmax_mm - (rmin_mm-rmin_mm./5))./(rad_bin_size);
%--------------------------------------------------------------------------
            nbin_radii(i-NUM_in+1) = (max_radius_mm.(im_name{i-NUM_in+1}) - min_radius_mm.(im_name{i-NUM_in+1}))./rad_bin_size; %#ok<SAGROW>
        end
        BSDX_median = zeros(1,round(nbin_radii_hg));
        rad_bin(1) = round(rmin_mm,1);
        for k = 1:round(nbin_radii_hg)
            rad_bin(k+1) = rad_bin(k)+ rad_bin_size;
        end
        BSDX_median(1) = rad_bin(1) + rad_bin_size/2;
        for k = 1:round(nbin_radii_hg)-1
            BSDX_median(k+1) = BSDX_median(k)+rad_bin_size;
        end
    %     BSDX_median = BSDX_median(2:end);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = NUM_in:NUM
    %% calculate volume of each bubble
            vol_b.(im_name{i-NUM_in+1}) = (4./3).*pi.*(radii_mm.(im_name{i-NUM_in+1})(:).^3); % volume of each bubble 
            vol_b.(im_name{i-NUM_in+1}) = sort(vol_b.(im_name{i-NUM_in+1}));  % sort in increasing order
        end
        for i = NUM_in:NUM
    %% sort volume and radii to ascending x-direction order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vol_b_xy.(im_name{i-NUM_in+1})(:,1) = vol_b.(im_name{i-NUM_in+1})(sort_idx_centers.(im_name{i-NUM_in+1})(:,1));  % sort vol according to centers sorting (in increasing x-dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            radii_mm_xy.(im_name{i-NUM_in+1})(:,1) = radii_mm.(im_name{i-NUM_in+1})(sort_idx_centers.(im_name{i-NUM_in+1})(:,1));  % sort radii according to centers sorting
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------------------STATISTICS---------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % allocate parameters:
        BSD_ROI0 = cell(1,NUM-NUM_in+1);
        BSD_ROI1 = cell(1,NUM-NUM_in+1);BSD_ROI2 = cell(1,NUM-NUM_in+1);BSD_ROI3 = cell(1,NUM-NUM_in+1);BSD_ROI4 = cell(1,NUM-NUM_in+1);
        BSD_ROI5 = cell(1,NUM-NUM_in+1);BSD_ROI6 = cell(1,NUM-NUM_in+1);BSD_ROI7 = cell(1,NUM-NUM_in+1);BSD_ROI8 = cell(1,NUM-NUM_in+1);
        BSD_ROI9 = cell(1,NUM-NUM_in+1);
        
        aux_VBSD_ROI0 = cell(1,NUM-NUM_in+1);
        aux_VBSD_ROI1 = cell(1,NUM-NUM_in+1);aux_VBSD_ROI2 = cell(1,NUM-NUM_in+1);aux_VBSD_ROI3 = cell(1,NUM-NUM_in+1);aux_VBSD_ROI4 = cell(1,NUM-NUM_in+1);
        aux_VBSD_ROI5 = cell(1,NUM-NUM_in+1);aux_VBSD_ROI6 = cell(1,NUM-NUM_in+1);aux_VBSD_ROI7 = cell(1,NUM-NUM_in+1);aux_VBSD_ROI8 = cell(1,NUM-NUM_in+1);
        aux_VBSD_ROI9 = cell(1,NUM-NUM_in+1);
        
        count_vol_ROI0 = cell(1,NUM-NUM_in+1);
        count_vol_ROI1 = cell(1,NUM-NUM_in+1);count_vol_ROI2 = cell(1,NUM-NUM_in+1);count_vol_ROI3 = cell(1,NUM-NUM_in+1);count_vol_ROI4 = cell(1,NUM-NUM_in+1);
        count_vol_ROI5 = cell(1,NUM-NUM_in+1);count_vol_ROI6 = cell(1,NUM-NUM_in+1);count_vol_ROI7 = cell(1,NUM-NUM_in+1);count_vol_ROI8 = cell(1,NUM-NUM_in+1);
        count_vol_ROI9 = cell(1,NUM-NUM_in+1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------------------create new ROI------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

%% ------------------------------------------------------------------------
        if strcmp(ROI{case_idx},'ROI6')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear bsd_sums 
            clear aux_BSD_ROI6_2 aux_BSD_ROI7_2 aux_BSD_ROI8_2 aux_BSD_ROI9_2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = NUM_in:NUM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcualate single valued functions per unit length for ROI5
    % conditional BSD for old_ROI5/4 segments
                for k = 2:round(nbin_radii_hg)+1 % iteration to find number of bubbles for each bin
                    count_ROI6 = 0; % no bubbles initially in each bin/each image
                    count_ROI7 = 0;
                    count_ROI8 = 0;
                    count_ROI9 = 0;
    %for each bin:
                    for j = 1:numel(centers_xy.(im_name{i-NUM_in+1})(:,2))
    % oldROI6 spans between 300-500 mm. so in newROI this is between 5oldROI/4-8oldROI/4 
    %choose which bubbles fall into specific center coordinates:
    %% same as above for different center coordinates:
                        if centers_xy.(im_name{i-NUM_in+1})(j,2) <= 6.*image_length_mm/4
    %% newROI6
    % BSD per unit length
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI6 = count_ROI6 + 1;
                                BSD_ROI6{i-NUM_in+1}(1,k-1) =  count_ROI6;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI6{i-NUM_in+1}(1,k-1) =  count_ROI6;  % add 0 if radius is not found for the kth bin
                            end

                         elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 6.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 7.*image_length_mm/4
    % repeat, as above for:
    %% newROI7              
    % BSD per unit length
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI7 = count_ROI7 + 1;
                                BSD_ROI7{i-NUM_in+1}(1,k-1) =  count_ROI7;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI7{i-NUM_in+1}(1,k-1) =  count_ROI7;  % add 0 if radius is not found for the kth bin
                            end

                          elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 7.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 8.*image_length_mm/4
    %% newROI8
    % BSD per unit length
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI8 = count_ROI8 + 1;
                                BSD_ROI8{i-NUM_in+1}(1,k-1) =  count_ROI8;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI8{i-NUM_in+1}(1,k-1) =  count_ROI8;  % add 0 if radius is not found for the kth bin
                            end
                            
                          elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 8.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 9.*image_length_mm/4
    %% newROI9
    % BSD per unit length
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI9 = count_ROI9 + 1;
                                BSD_ROI9{i-NUM_in+1}(1,k-1) =  count_ROI9;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI9{i-NUM_in+1}(1,k-1) =  count_ROI9;  % add 0 if radius is not found for the kth bin
                            end
                        end
                    end
                end
            end
            for i=NUM_in:NUM
    % check for empty arrays and replace with zeros
                    if isempty(BSD_ROI6{i-NUM_in+1}) 
                        BSD_ROI6{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI7{i-NUM_in+1}) 
                        BSD_ROI7{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI8{i-NUM_in+1}) 
                        BSD_ROI8{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI9{i-NUM_in+1}) 
                        BSD_ROI9{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
            end
            for i=NUM_in:NUM
%% VBSD per unit length. Method:
% 0.  vol_b_xy is already sorted in ascending x-dir. Therefore, total number of bubbles N per ROI (sum BSD_ROI) can sort volumes per ROI.
% 1.  sort volume of bubbles per ROI (result: vol_b_xy_ROI). Sorting is in ascending order of volume.
% 2.  for each BSD radius bin the number of bubbles is known (count_vol_ROI). 
% 3.  therefore, the sorted volumes (vol_b_xy_ROI) can be sorted even further for each bin and the volume per bin can be found (aux_VBSD)
% 4.  this is done per image. aux_VBSD is for each image. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (think twice what this is)
                bsd_sums = [1,sum(BSD_ROI6{1,i-NUM_in+1}),sum(BSD_ROI6{1,i-NUM_in+1}+BSD_ROI7{1,i-NUM_in+1}),...
                sum(BSD_ROI6{1,i-NUM_in+1}+BSD_ROI7{1,i-NUM_in+1}+BSD_ROI8{1,i-NUM_in+1}),...
                sum(BSD_ROI6{1,i-NUM_in+1}+BSD_ROI7{1,i-NUM_in+1}+BSD_ROI8{1,i-NUM_in+1}+BSD_ROI9{1,i-NUM_in+1})]; % indices to calculate Volume for each ROI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ascending order of volume of bubbles per ROI:
                vol_b_xy_ROI6.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2)));
                vol_b_xy_ROI7.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3)));
                vol_b_xy_ROI8.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4)));
                vol_b_xy_ROI9.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5)));
    % 3-4. calculate VBSD for each new ROI/each bin/each image 
                count_vol_ROI6{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI5
                count_vol_ROI7{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI6
                count_vol_ROI8{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI7
                count_vol_ROI9{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI8

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  VBSD is as sum of volume of each bubble. Each bin radius is calculated separately:
                for k = 2:numel(BSDX_median)+1
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % find the indices of bubbles in BSD matrices that belong to same bins  
    % calculate the sum of volumes of bubbles for each bin radius:
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                    count_vol_ROI6{i-NUM_in+1}(1,k)  = count_vol_ROI6{i-NUM_in+1}(1,k-1) + BSD_ROI6{i-NUM_in+1}(1,k-1);   

                    aux_VBSD_ROI6{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI6.(im_name{i-NUM_in+1})...
                        (count_vol_ROI6{i-NUM_in+1}(1,k-1):count_vol_ROI6{i-NUM_in+1}(1,k)-1));  
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------                
                    count_vol_ROI7{i-NUM_in+1}(1,k)  = count_vol_ROI7{i-NUM_in+1}(1,k-1) + BSD_ROI7{i-NUM_in+1}(1,k-1);  

                    aux_VBSD_ROI7{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI7.(im_name{i-NUM_in+1})...
                    (count_vol_ROI7{i-NUM_in+1}(1,k-1):count_vol_ROI7{i-NUM_in+1}(1,k)-1));  
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------                
                    count_vol_ROI8{i-NUM_in+1}(1,k)  = count_vol_ROI8{i-NUM_in+1}(1,k-1) + BSD_ROI8{i-NUM_in+1}(1,k-1);  

                    aux_VBSD_ROI8{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI8.(im_name{i-NUM_in+1})...
                        (count_vol_ROI8{i-NUM_in+1}(1,k-1):count_vol_ROI8{i-NUM_in+1}(1,k)-1)); 
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                    count_vol_ROI9{i-NUM_in+1}(1,k)  = count_vol_ROI9{i-NUM_in+1}(1,k-1) + BSD_ROI9{i-NUM_in+1}(1,k-1);  

                    aux_VBSD_ROI9{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI9.(im_name{i-NUM_in+1})...
                        (count_vol_ROI9{i-NUM_in+1}(1,k-1):count_vol_ROI9{i-NUM_in+1}(1,k)-1)); 
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate instanteneous void fraction for ROI1-8
                aux_alpha_ins_ROI6{i-NUM_in+1} = sum(aux_VBSD_ROI6{i-NUM_in+1});
                aux_alpha_ins_ROI7{i-NUM_in+1} = sum(aux_VBSD_ROI7{i-NUM_in+1});
                aux_alpha_ins_ROI8{i-NUM_in+1} = sum(aux_VBSD_ROI8{i-NUM_in+1});
                aux_alpha_ins_ROI9{i-NUM_in+1} = sum(aux_VBSD_ROI9{i-NUM_in+1});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% for comparison, calculate tot vol_bubbles per ROI in another way:
    %         for i=NUM_in:NUM
    %             vol_b_tot_ROI6{i-NUM_in+1} = sum(vol_b_xy_ROI6.(im_name{i-NUM_in+1}));
    %             vol_b_tot_ROI7{i-NUM_in+1} = sum(vol_b_xy_ROI7.(im_name{i-NUM_in+1}));
    %             vol_b_tot_ROI8{i-NUM_in+1} = sum(vol_b_xy_ROI8.(im_name{i-NUM_in+1}));
    %             vol_b_tot_ROI9{i-NUM_in+1} = sum(vol_b_xy_ROI9.(im_name{i-NUM_in+1}));

    %         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     save(['BSD_avg_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI6_2','-v7.3');   
    %     save(['BSD_avg_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI7_2','-v7.3'); 
    %     save(['BSD_avg_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI8_2','-v7.3');  
    %     save(['BSD_avg_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI9_2','-v7.3');  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=NUM_in:NUM
    % store BSD results for each new ROI for each image 
                for k = 1:numel(BSDX_median)
                    aux_BSD_ROI6(i-NUM_in+1,k) = (BSD_ROI6{1,i-NUM_in+1}(k));  
                    aux_BSD_ROI7(i-NUM_in+1,k) = (BSD_ROI7{1,i-NUM_in+1}(k));  
                    aux_BSD_ROI8(i-NUM_in+1,k) = (BSD_ROI8{1,i-NUM_in+1}(k));  
                    aux_BSD_ROI9(i-NUM_in+1,k) = (BSD_ROI9{1,i-NUM_in+1}(k));  

                end
            end
            for k = 1:numel(BSDX_median)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% find time averaged distribution for each case. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                aux_BSD_ROI6_2(1,k) = mean(aux_BSD_ROI6(:,k));
                aux_BSD_ROI7_2(1,k) = mean(aux_BSD_ROI7(:,k));
                aux_BSD_ROI8_2(1,k) = mean(aux_BSD_ROI8(:,k));
                aux_BSD_ROI9_2(1,k) = mean(aux_BSD_ROI9(:,k));

            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            ROI6_idx = ROI6_idx+1; % ROI6 index IMPORTANT!! (represents oldROI)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% store in ROI1-8.ROI-5 type of struct
            BSD_avg_ROI6(ROI6_idx).ROI5= aux_BSD_ROI6_2(:)';   
            BSD_avg_ROI7(ROI6_idx).ROI5 = aux_BSD_ROI7_2(:)';   
            BSD_avg_ROI8(ROI6_idx).ROI5 = aux_BSD_ROI8_2(:)';
            BSD_avg_ROI9(ROI6_idx).ROI5 = aux_BSD_ROI9_2(:)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% time average for VBSD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_VBSD_ROI6_2(i-NUM_in+1,k) = (aux_VBSD_ROI6{1,i-NUM_in+1}(k)); 
                aux_VBSD_ROI7_2(i-NUM_in+1,k) = (aux_VBSD_ROI7{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI8_2(i-NUM_in+1,k) = (aux_VBSD_ROI8{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI9_2(i-NUM_in+1,k) = (aux_VBSD_ROI9{1,i-NUM_in+1}(k));  % place VBSD values of each image/rep. in a 2D matrix 
            end
        end
%% time averaged VBSD:
        for k = 1:numel(BSDX_median)
                VBSD_ROI6_avg_2(1,k) = mean(aux_VBSD_ROI6_2(:,k));
                VBSD_ROI7_avg_2(1,k) = mean(aux_VBSD_ROI7_2(:,k));
                VBSD_ROI8_avg_2(1,k) = mean(aux_VBSD_ROI8_2(:,k));
                VBSD_ROI9_avg_2(1,k) = mean(aux_VBSD_ROI9_2(:,k));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find time averages of specific time slots of videos:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% aux_VBSD_2 is equivalent to aux_BSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                VBSD_ROI6_act_1(1,k) = mean(aux_VBSD_ROI6_2(1:round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_1(1,k) = mean(aux_VBSD_ROI7_2(1:round(numel(aux_VBSD_ROI7_2(:,k))/4),k));
                VBSD_ROI8_act_1(1,k) = mean(aux_VBSD_ROI8_2(1:round(numel(aux_VBSD_ROI8_2(:,k))/4),k));
                VBSD_ROI9_act_1(1,k) = mean(aux_VBSD_ROI9_2(1:round(numel(aux_VBSD_ROI9_2(:,k))/4),k));


                VBSD_ROI6_act_2(1,k) = mean(aux_VBSD_ROI6_2(round(numel(aux_VBSD_ROI6_2(:,k))/4):2*round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_2(1,k) = mean(aux_VBSD_ROI7_2(round(numel(aux_VBSD_ROI7_2(:,k))/4):2*round(numel(aux_VBSD_ROI7_2(:,k))/4),k));
                VBSD_ROI8_act_2(1,k) = mean(aux_VBSD_ROI8_2(round(numel(aux_VBSD_ROI8_2(:,k))/4):2*round(numel(aux_VBSD_ROI8_2(:,k))/4),k));
                VBSD_ROI9_act_2(1,k) = mean(aux_VBSD_ROI9_2(round(numel(aux_VBSD_ROI9_2(:,k))/4):2*round(numel(aux_VBSD_ROI9_2(:,k))/4),k));


                VBSD_ROI6_act_3(1,k) = mean(aux_VBSD_ROI6_2(2*round(numel(aux_VBSD_ROI6_2(:,k))/4):3*round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_3(1,k) = mean(aux_VBSD_ROI7_2(2*round(numel(aux_VBSD_ROI7_2(:,k))/4):3*round(numel(aux_VBSD_ROI7_2(:,k))/4),k));
                VBSD_ROI8_act_3(1,k) = mean(aux_VBSD_ROI8_2(2*round(numel(aux_VBSD_ROI8_2(:,k))/4):3*round(numel(aux_VBSD_ROI8_2(:,k))/4),k));
                VBSD_ROI9_act_3(1,k) = mean(aux_VBSD_ROI9_2(2*round(numel(aux_VBSD_ROI9_2(:,k))/4):3*round(numel(aux_VBSD_ROI9_2(:,k))/4),k));


                VBSD_ROI6_act_4(1,k) = mean(aux_VBSD_ROI6_2(3*round(numel(aux_VBSD_ROI6_2(:,k))/4):4*(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_4(1,k) = mean(aux_VBSD_ROI7_2(3*round(numel(aux_VBSD_ROI7_2(:,k))/4):4*(numel(aux_VBSD_ROI7_2(:,k))/4),k));
                VBSD_ROI8_act_4(1,k) = mean(aux_VBSD_ROI8_2(3*round(numel(aux_VBSD_ROI8_2(:,k))/4):4*(numel(aux_VBSD_ROI8_2(:,k))/4),k));
                VBSD_ROI9_act_4(1,k) = mean(aux_VBSD_ROI9_2(3*round(numel(aux_VBSD_ROI9_2(:,k))/4):4*(numel(aux_VBSD_ROI9_2(:,k))/4),k));
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                BSD_ROI6_act_1(1,k) = mean(aux_BSD_ROI6(1:round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_1(1,k) = mean(aux_BSD_ROI7(1:round(numel(aux_BSD_ROI7(:,k))/4),k));
                BSD_ROI8_act_1(1,k) = mean(aux_BSD_ROI8(1:round(numel(aux_BSD_ROI8(:,k))/4),k));
                BSD_ROI9_act_1(1,k) = mean(aux_BSD_ROI9(1:round(numel(aux_BSD_ROI9(:,k))/4),k));



                BSD_ROI6_act_2(1,k) = mean(aux_BSD_ROI6(round(numel(aux_BSD_ROI6(:,k))/4):2*round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_2(1,k) = mean(aux_BSD_ROI7(round(numel(aux_BSD_ROI7(:,k))/4):2*round(numel(aux_BSD_ROI7(:,k))/4),k));
                BSD_ROI8_act_2(1,k) = mean(aux_BSD_ROI8(round(numel(aux_BSD_ROI8(:,k))/4):2*round(numel(aux_BSD_ROI8(:,k))/4),k));
                BSD_ROI9_act_2(1,k) = mean(aux_BSD_ROI9(round(numel(aux_BSD_ROI9(:,k))/4):2*round(numel(aux_BSD_ROI9(:,k))/4),k));


                BSD_ROI6_act_3(1,k) = mean(aux_BSD_ROI6(2*round(numel(aux_BSD_ROI6(:,k))/4):3*round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_3(1,k) = mean(aux_BSD_ROI7(2*round(numel(aux_BSD_ROI7(:,k))/4):3*round(numel(aux_BSD_ROI7(:,k))/4),k));
                BSD_ROI8_act_3(1,k) = mean(aux_BSD_ROI8(2*round(numel(aux_BSD_ROI8(:,k))/4):3*round(numel(aux_BSD_ROI8(:,k))/4),k));
                BSD_ROI9_act_3(1,k) = mean(aux_BSD_ROI9(2*round(numel(aux_BSD_ROI9(:,k))/4):3*round(numel(aux_BSD_ROI9(:,k))/4),k));


                BSD_ROI6_act_4(1,k) = mean(aux_BSD_ROI6(3*round(numel(aux_BSD_ROI6(:,k))/4):4*(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_4(1,k) = mean(aux_BSD_ROI7(3*round(numel(aux_BSD_ROI7(:,k))/4):4*(numel(aux_BSD_ROI7(:,k))/4),k));
                BSD_ROI8_act_4(1,k) = mean(aux_BSD_ROI8(3*round(numel(aux_BSD_ROI8(:,k))/4):4*(numel(aux_BSD_ROI8(:,k))/4),k));
                BSD_ROI9_act_4(1,k) = mean(aux_BSD_ROI9(3*round(numel(aux_BSD_ROI9(:,k))/4):4*(numel(aux_BSD_ROI9(:,k))/4),k));

        end

    %% store in VBSD_ROI(1-8).oldROI(1-5) type of struct
    %% which is used later to calculate ensenble averages (many runs)  
            VBSD_avg_ROI6(ROI6_idx).ROI6 = VBSD_ROI6_avg_2(:)';   
            VBSD_avg_ROI7(ROI6_idx).ROI6 = VBSD_ROI7_avg_2(:)';  
            VBSD_avg_ROI8(ROI6_idx).ROI6 = VBSD_ROI8_avg_2(:)'; 
            VBSD_avg_ROI9(ROI6_idx).ROI6 = VBSD_ROI9_avg_2(:)';   


            VBSD_act_1_ROI6(ROI6_idx).ROI6 = VBSD_ROI6_act_1(:)';   
            VBSD_act_1_ROI7(ROI6_idx).ROI6 = VBSD_ROI7_act_1(:)'; 
            VBSD_act_1_ROI8(ROI6_idx).ROI6 = VBSD_ROI8_act_1(:)';  
            VBSD_act_1_ROI9(ROI6_idx).ROI6 = VBSD_ROI9_act_1(:)';  


            VBSD_act_2_ROI6(ROI6_idx).ROI6 = VBSD_ROI6_act_2(:)';   
            VBSD_act_2_ROI7(ROI6_idx).ROI6 = VBSD_ROI7_act_2(:)'; 
            VBSD_act_2_ROI8(ROI6_idx).ROI6 = VBSD_ROI8_act_2(:)';  
            VBSD_act_2_ROI9(ROI6_idx).ROI6 = VBSD_ROI9_act_2(:)';  


            VBSD_act_3_ROI6(ROI6_idx).ROI6 = VBSD_ROI6_act_3(:)';   
            VBSD_act_3_ROI7(ROI6_idx).ROI6 = VBSD_ROI7_act_3(:)'; 
            VBSD_act_3_ROI8(ROI6_idx).ROI6 = VBSD_ROI8_act_3(:)'; 
            VBSD_act_3_ROI9(ROI6_idx).ROI6 = VBSD_ROI9_act_3(:)';  


            VBSD_act_4_ROI6(ROI6_idx).ROI6 = VBSD_ROI6_act_4(:)';   
            VBSD_act_4_ROI7(ROI6_idx).ROI6 = VBSD_ROI7_act_4(:)'; 
            VBSD_act_4_ROI8(ROI6_idx).ROI6 = VBSD_ROI8_act_4(:)'; 
            VBSD_act_4_ROI9(ROI6_idx).ROI6 = VBSD_ROI9_act_4(:)';  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            BSD_act_1_ROI6(ROI6_idx).ROI6 = BSD_ROI6_act_1(:)';   
            BSD_act_1_ROI7(ROI6_idx).ROI6 = BSD_ROI7_act_1(:)'; 
            BSD_act_1_ROI8(ROI6_idx).ROI6 = BSD_ROI8_act_1(:)';  
            BSD_act_1_ROI9(ROI6_idx).ROI6 = BSD_ROI9_act_1(:)';  

            
            BSD_act_2_ROI6(ROI6_idx).ROI6 = BSD_ROI6_act_2(:)';   
            BSD_act_2_ROI7(ROI6_idx).ROI6 = BSD_ROI7_act_2(:)'; 
            BSD_act_2_ROI8(ROI6_idx).ROI6 = BSD_ROI8_act_2(:)';  
            BSD_act_2_ROI9(ROI6_idx).ROI6 = BSD_ROI9_act_2(:)';  

            
            BSD_act_3_ROI6(ROI6_idx).ROI6 = BSD_ROI6_act_3(:)';   
            BSD_act_3_ROI7(ROI6_idx).ROI6 = BSD_ROI7_act_3(:)'; 
            BSD_act_3_ROI8(ROI6_idx).ROI6 = BSD_ROI8_act_3(:)';  
            BSD_act_3_ROI9(ROI6_idx).ROI6 = BSD_ROI9_act_3(:)';  

            
            BSD_act_4_ROI6(ROI6_idx).ROI6 = BSD_ROI6_act_4(:)';   
            BSD_act_4_ROI7(ROI6_idx).ROI6 = BSD_ROI7_act_4(:)'; 
            BSD_act_4_ROI8(ROI6_idx).ROI6 = BSD_ROI8_act_4(:)';  
            BSD_act_4_ROI9(ROI6_idx).ROI6 = BSD_ROI9_act_4(:)';  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% store results in matrices for further analysis 
        if exist('aux_BSD_ROI6_2','var') == 1
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
            ROI6_mtx_idx = ROI6_mtx_idx + 1;  % IMPORTANT!!
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
            for k = 1:numel(BSDX_median)
    % place BSD in matrix to use for ensemble avrgs
                BSD_avg_ROI6_mtx(ROI6_mtx_idx,k) = aux_BSD_ROI6_2(1,k);
                VBSD_avg_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_avg_2(1,k);

                VBSD_act_1_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_1(1,k);
                VBSD_act_2_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_2(1,k);
                VBSD_act_3_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_3(1,k);
                VBSD_act_4_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_4(1,k);

                BSD_act_1_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_1(1,k);
                BSD_act_2_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_2(1,k);
                BSD_act_3_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_3(1,k);
                BSD_act_4_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_4(1,k);
            end
            for i = NUM_in : NUM
                alpha_ins_ROI6(ROI6_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI6{i-NUM_in+1};
            end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI6_mtx(ROI6_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI6_2(i-NUM_in+1,k); 
                end
            end
        end
        if exist('aux_BSD_ROI7_2','var') == 1
            ROI7_mtx_idx = ROI7_mtx_idx + 1;  % IMPORTANT!!
            for k = 1:numel(BSDX_median)
    % place BSD in matrix to use for ensemble avrgs
                BSD_avg_ROI7_mtx(ROI7_mtx_idx,k) = aux_BSD_ROI7_2(1,k);
                VBSD_avg_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_avg_2(1,k);

                VBSD_act_1_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_1(1,k);
                VBSD_act_2_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_2(1,k);
                VBSD_act_3_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_3(1,k);
                VBSD_act_4_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_4(1,k);

                BSD_act_1_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_1(1,k);
                BSD_act_2_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_2(1,k);
                BSD_act_3_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_3(1,k);
                BSD_act_4_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_4(1,k);
            end
            for i = NUM_in : NUM
                alpha_ins_ROI7(ROI7_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI7{i-NUM_in+1};
            end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI7_mtx(ROI7_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI7_2(i-NUM_in+1,k); 
                end
            end
        end
        if exist('aux_BSD_ROI8_2','var') == 1
            ROI8_mtx_idx = ROI8_mtx_idx + 1;  % IMPORTANT!!
            for k = 1:numel(BSDX_median)
                BSD_avg_ROI8_mtx(ROI8_mtx_idx,k) = aux_BSD_ROI8_2(1,k);
                VBSD_avg_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_avg_2(1,k);
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                VBSD_act_1_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_act_1(1,k);
                VBSD_act_2_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_act_2(1,k);
                VBSD_act_3_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_act_3(1,k);
                VBSD_act_4_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_act_4(1,k);

                BSD_act_1_ROI8_mtx(ROI8_mtx_idx,k) = BSD_ROI8_act_1(1,k);
                BSD_act_2_ROI8_mtx(ROI8_mtx_idx,k) = BSD_ROI8_act_2(1,k);
                BSD_act_3_ROI8_mtx(ROI8_mtx_idx,k) = BSD_ROI8_act_3(1,k);
                BSD_act_4_ROI8_mtx(ROI8_mtx_idx,k) = BSD_ROI8_act_4(1,k);
            end
            for i = NUM_in : NUM
                alpha_ins_ROI8(ROI8_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI8{i-NUM_in+1};
            end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI8_mtx(ROI8_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI8_2(i-NUM_in+1,k); 
                end
            end
        end  
        if exist('aux_BSD_ROI9_2','var') == 1
            ROI9_mtx_idx = ROI9_mtx_idx + 1;  % IMPORTANT!!
            for k = 1:numel(BSDX_median)
                BSD_avg_ROI9_mtx(ROI9_mtx_idx,k) = aux_BSD_ROI9_2(1,k);
                VBSD_avg_ROI9_mtx(ROI9_mtx_idx,k) = VBSD_ROI9_avg_2(1,k);
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                VBSD_act_1_ROI9_mtx(ROI9_mtx_idx,k) = VBSD_ROI9_act_1(1,k);
                VBSD_act_2_ROI9_mtx(ROI9_mtx_idx,k) = VBSD_ROI9_act_2(1,k);
                VBSD_act_3_ROI9_mtx(ROI9_mtx_idx,k) = VBSD_ROI9_act_3(1,k);
                VBSD_act_4_ROI9_mtx(ROI9_mtx_idx,k) = VBSD_ROI9_act_4(1,k);

                BSD_act_1_ROI9_mtx(ROI9_mtx_idx,k) = BSD_ROI9_act_1(1,k);
                BSD_act_2_ROI9_mtx(ROI9_mtx_idx,k) = BSD_ROI9_act_2(1,k);
                BSD_act_3_ROI9_mtx(ROI9_mtx_idx,k) = BSD_ROI9_act_3(1,k);
                BSD_act_4_ROI9_mtx(ROI9_mtx_idx,k) = BSD_ROI9_act_4(1,k);
            end
            for i = NUM_in : NUM
                alpha_ins_ROI9(ROI9_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI9{i-NUM_in+1};
            end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI9_mtx(ROI9_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI9_2(i-NUM_in+1,k); 
                end
            end
        end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % save results 
            cd(dir_save)
            save('BSDX_median.mat','BSDX_median','-v7.3');
            save(['BSD_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI6','-v7.3');   
            save(['BSD_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI7','-v7.3'); 
            save(['BSD_ins_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI8','-v7.3');  
            save(['BSD_ins_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI9','-v7.3');  

    % save results
            save(['VBSD_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI6_2','-v7.3');   
            save(['VBSD_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI7_2','-v7.3'); 
            save(['VBSD_ins_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI8_2','-v7.3'); 
            save(['VBSD_ins_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI9_2','-v7.3'); 

            
    %     save(['VBSD_avg_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI9_avg_2','-v7.3');   

    % 
    %     save(['VBSD_act_1_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_1','-v7.3');
    %     save(['VBSD_act_1_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI9_act_1','-v7.3');   

    %     
    %     save(['VBSD_act_2_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_2','-v7.3');
    %     save(['VBSD_act_2_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI9_act_2','-v7.3');   

    %     
    %     save(['VBSD_act_3_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_3','-v7.3');
    %     save(['VBSD_act_3_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI9_act_3','-v7.3');   

    %         
    %     save(['VBSD_act_4_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI9_act_4','-v7.3');   

    
            save(['alpha_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI6','-v7.3');   
            save(['alpha_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI7','-v7.3');   
            save(['alpha_ins_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI8','-v7.3');   
            save(['alpha_ins_ROI9_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI9','-v7.3');   

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(ROI{case_idx},'ROI5')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear bsd_sums aux_BSD_ROI5_2 aux_BSD_ROI6_2 aux_BSD_ROI7_2 aux_BSD_ROI8_2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = NUM_in:NUM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcualate single valued functions per unit length for ROI5
    % conditional BSD for old_ROI5/4 segments
                for k = 2:round(nbin_radii_hg)+1 % iteration to find number of bubbles for each bin
                    count_ROI5 = 0; % no bubbles initially in each bin/each image
                    count_ROI6 = 0;
                    count_ROI7 = 0;
                    count_ROI8 = 0;
    %for each bin:
                    for j = 1:numel(centers_xy.(im_name{i-NUM_in+1})(:,2))
    % oldROI5 spans between 200-400 mm. so in newROI this is between 5oldROI/4-8oldROI/4 
    %choose which bubbles fall into specific center coordinates:
                        if centers_xy.(im_name{i-NUM_in+1})(j,2) <= 5.*image_length_mm/4 
    %% new ROI5:
    %check if radius belongs in the jth rad_bin
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble RADIUS falls within the Kth bin size
    %if yes, count +1
                                count_ROI5 = count_ROI5 + 1;
                                BSD_ROI5{i-NUM_in+1}(1,k-1) =  count_ROI5;  % COUNT +1 for each detected radius in the kth bin
    %                             BSD_idx_ROI5{i-NUM_in+1}(j,k-1) = j; % store indices of radii per bin 
                            else
    %if not, do not 
                                BSD_ROI5{i-NUM_in+1}(1,k-1) =  count_ROI5;  % add 0 if radius is not found for the kth bin
                            end
    %same as above for different center coordinates:
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 5.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 6.*image_length_mm/4
    %newROI6
    % BSD per unit length
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI6 = count_ROI6 + 1;
                                BSD_ROI6{i-NUM_in+1}(1,k-1) =  count_ROI6;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI6{i-NUM_in+1}(1,k-1) =  count_ROI6;  % add 0 if radius is not found for the kth bin
                            end

                         elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 6.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 7.*image_length_mm/4
    % repeat, as above for:
     %newROI7              
    % BSD per unit length
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI7 = count_ROI7 + 1;
                                BSD_ROI7{i-NUM_in+1}(1,k-1) =  count_ROI7;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI7{i-NUM_in+1}(1,k-1) =  count_ROI7;  % add 0 if radius is not found for the kth bin
                            end

                          elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 7.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 8.*image_length_mm/4
    %newROI8
    % BSD per unit length
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI8 = count_ROI8 + 1;
                                BSD_ROI8{i-NUM_in+1}(1,k-1) =  count_ROI8;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI8{i-NUM_in+1}(1,k-1) =  count_ROI8;  % add 0 if radius is not found for the kth bin
                            end
                        end
                    end
                end
            end
            for i=NUM_in:NUM
    % check for empty arrays and replace with zeros
                    if isempty(BSD_ROI5{i-NUM_in+1}) 
                        BSD_ROI5{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI6{i-NUM_in+1}) 
                        BSD_ROI6{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI7{i-NUM_in+1}) 
                        BSD_ROI7{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI8{i-NUM_in+1}) 
                        BSD_ROI8{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
            end
            for i=NUM_in:NUM
%% VBSD per unit length. Method:
% 0.  vol_b_xy is already sorted in ascending x-dir. Therefore, total number of bubbles N per ROI (sum BSD_ROI) can sort volumes per ROI.
% 1.  sort volume of bubbles per ROI (result: vol_b_xy_ROI). Sorting is in ascending order of volume.
% 2.  for each BSD radius bin the number of bubbles is known (count_vol_ROI). 
% 3.  therefore, the sorted volumes (vol_b_xy_ROI) can be sorted even further for each bin and the volume per bin can be found (aux_VBSD)
% 4.  this is done per image. aux_VBSD is for each image. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (think twice what this is)
                bsd_sums = [1,sum(BSD_ROI5{1,i-NUM_in+1}),sum(BSD_ROI5{1,i-NUM_in+1}+BSD_ROI6{1,i-NUM_in+1}),...
                sum(BSD_ROI5{1,i-NUM_in+1}+BSD_ROI6{1,i-NUM_in+1}+BSD_ROI7{1,i-NUM_in+1}),...
                sum(BSD_ROI5{1,i-NUM_in+1}+BSD_ROI6{1,i-NUM_in+1}+BSD_ROI7{1,i-NUM_in+1}+BSD_ROI8{1,i-NUM_in+1})]; % indices to calculate Volume for each ROI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ascending order of volume of bubbles per ROI:
                vol_b_xy_ROI5.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2)));
                vol_b_xy_ROI6.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3)));
                vol_b_xy_ROI7.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4)));
                vol_b_xy_ROI8.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5)));
    % 3-4. calculate VBSD for each new ROI/each bin/each image 
                count_vol_ROI5{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI5
                count_vol_ROI6{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI6
                count_vol_ROI7{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI7
                count_vol_ROI8{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI8

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  VBSD is as sum of volume of each bubble. Each bin radius is calculated separately:
                for k = 2:numel(BSDX_median)+1
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % find the indices of bubbles in BSD matrices that belong to same bins  
                    count_vol_ROI5{i-NUM_in+1}(1,k)  = count_vol_ROI5{i-NUM_in+1}(1,k-1) + BSD_ROI5{i-NUM_in+1}(1,k-1); 
    % calculate the sum of volumes of bubbles for each bin radius:
                    aux_VBSD_ROI5{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI5.(im_name{i-NUM_in+1})...
                        (count_vol_ROI5{i-NUM_in+1}(1,k-1):count_vol_ROI5{i-NUM_in+1}(1,k)-1)); 
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                    count_vol_ROI6{i-NUM_in+1}(1,k)  = count_vol_ROI6{i-NUM_in+1}(1,k-1) + BSD_ROI6{i-NUM_in+1}(1,k-1);   

                    aux_VBSD_ROI6{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI6.(im_name{i-NUM_in+1})...
                        (count_vol_ROI6{i-NUM_in+1}(1,k-1):count_vol_ROI6{i-NUM_in+1}(1,k)-1));  
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------                
                    count_vol_ROI7{i-NUM_in+1}(1,k)  = count_vol_ROI7{i-NUM_in+1}(1,k-1) + BSD_ROI7{i-NUM_in+1}(1,k-1);  

                    aux_VBSD_ROI7{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI7.(im_name{i-NUM_in+1})...
                    (count_vol_ROI7{i-NUM_in+1}(1,k-1):count_vol_ROI7{i-NUM_in+1}(1,k)-1));  
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------                
                    count_vol_ROI8{i-NUM_in+1}(1,k)  = count_vol_ROI8{i-NUM_in+1}(1,k-1) + BSD_ROI8{i-NUM_in+1}(1,k-1);  

                    aux_VBSD_ROI8{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI8.(im_name{i-NUM_in+1})...
                        (count_vol_ROI8{i-NUM_in+1}(1,k-1):count_vol_ROI8{i-NUM_in+1}(1,k)-1)); 
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate instanteneous void fraction for ROI1-8
                aux_alpha_ins_ROI5{i-NUM_in+1} = sum(aux_VBSD_ROI5{i-NUM_in+1});
                aux_alpha_ins_ROI6{i-NUM_in+1} = sum(aux_VBSD_ROI6{i-NUM_in+1});
                aux_alpha_ins_ROI7{i-NUM_in+1} = sum(aux_VBSD_ROI7{i-NUM_in+1});
                aux_alpha_ins_ROI8{i-NUM_in+1} = sum(aux_VBSD_ROI8{i-NUM_in+1});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% for comparison, calculate tot vol_bubbles per ROI in another way:
    %         for i=NUM_in:NUM
    %             vol_b_tot_ROI5{i-NUM_in+1} = sum(vol_b_xy_ROI5.(im_name{i-NUM_in+1}));
    %             vol_b_tot_ROI6{i-NUM_in+1} = sum(vol_b_xy_ROI6.(im_name{i-NUM_in+1}));
    %             vol_b_tot_ROI7{i-NUM_in+1} = sum(vol_b_xy_ROI7.(im_name{i-NUM_in+1}));
    %             vol_b_tot_ROI8{i-NUM_in+1} = sum(vol_b_xy_ROI8.(im_name{i-NUM_in+1}));
    %         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     save(['BSD_avg_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI5_2','-v7.3');   
    %     save(['BSD_avg_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI6_2','-v7.3');   
    %     save(['BSD_avg_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI7_2','-v7.3'); 
    %     save(['BSD_avg_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI8_2','-v7.3');  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=NUM_in:NUM
    % store BSD results for each new ROI for each image 
                for k = 1:numel(BSDX_median)
                    aux_BSD_ROI5(i-NUM_in+1,k) = (BSD_ROI5{1,i-NUM_in+1}(k));   
                    aux_BSD_ROI6(i-NUM_in+1,k) = (BSD_ROI6{1,i-NUM_in+1}(k));  
                    aux_BSD_ROI7(i-NUM_in+1,k) = (BSD_ROI7{1,i-NUM_in+1}(k));  
                    aux_BSD_ROI8(i-NUM_in+1,k) = (BSD_ROI8{1,i-NUM_in+1}(k));  
                end
            end
            for k = 1:numel(BSDX_median)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% find time averaged distribution for each case. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                aux_BSD_ROI5_2(1,k) = mean(aux_BSD_ROI5(:,k));
                aux_BSD_ROI6_2(1,k) = mean(aux_BSD_ROI6(:,k));
                aux_BSD_ROI7_2(1,k) = mean(aux_BSD_ROI7(:,k));
                aux_BSD_ROI8_2(1,k) = mean(aux_BSD_ROI8(:,k));
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            ROI5_idx = ROI5_idx+1; % ROI5 index  IMPORTANT!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% store in ROI1-8.ROI-5 type of struct
            BSD_avg_ROI5(ROI5_idx).ROI5 = aux_BSD_ROI5_2(:)';   
            BSD_avg_ROI6(ROI5_idx).ROI5= aux_BSD_ROI6_2(:)';   
            BSD_avg_ROI7(ROI5_idx).ROI5 = aux_BSD_ROI7_2(:)';   
            BSD_avg_ROI8(ROI5_idx).ROI5 = aux_BSD_ROI8_2(:)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% time average for VBSD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_VBSD_ROI5_2(i-NUM_in+1,k) = (aux_VBSD_ROI5{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI6_2(i-NUM_in+1,k) = (aux_VBSD_ROI6{1,i-NUM_in+1}(k)); 
                aux_VBSD_ROI7_2(i-NUM_in+1,k) = (aux_VBSD_ROI7{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI8_2(i-NUM_in+1,k) = (aux_VBSD_ROI8{1,i-NUM_in+1}(k));  % place VBSD values of each image/rep. in a 2D matrix 
            end
        end
%% time averaged VBSD:
        for k = 1:numel(BSDX_median)
                VBSD_ROI5_avg_2(1,k) = mean(aux_VBSD_ROI5_2(:,k));
                VBSD_ROI6_avg_2(1,k) = mean(aux_VBSD_ROI6_2(:,k));
                VBSD_ROI7_avg_2(1,k) = mean(aux_VBSD_ROI7_2(:,k));
                VBSD_ROI8_avg_2(1,k) = mean(aux_VBSD_ROI8_2(:,k));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% find time averages of specific time slots of videos:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% aux_VBSD_2 is equivalent to aux_BSD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                VBSD_ROI5_act_1(1,k) = mean(aux_VBSD_ROI5_2(1:round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_1(1,k) = mean(aux_VBSD_ROI6_2(1:round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_1(1,k) = mean(aux_VBSD_ROI7_2(1:round(numel(aux_VBSD_ROI7_2(:,k))/4),k));
                VBSD_ROI8_act_1(1,k) = mean(aux_VBSD_ROI8_2(1:round(numel(aux_VBSD_ROI8_2(:,k))/4),k));


                VBSD_ROI5_act_2(1,k) = mean(aux_VBSD_ROI5_2(round(numel(aux_VBSD_ROI5_2(:,k))/4):2*round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_2(1,k) = mean(aux_VBSD_ROI6_2(round(numel(aux_VBSD_ROI6_2(:,k))/4):2*round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_2(1,k) = mean(aux_VBSD_ROI7_2(round(numel(aux_VBSD_ROI7_2(:,k))/4):2*round(numel(aux_VBSD_ROI7_2(:,k))/4),k));
                VBSD_ROI8_act_2(1,k) = mean(aux_VBSD_ROI8_2(round(numel(aux_VBSD_ROI8_2(:,k))/4):2*round(numel(aux_VBSD_ROI8_2(:,k))/4),k));


                VBSD_ROI5_act_3(1,k) = mean(aux_VBSD_ROI5_2(2*round(numel(aux_VBSD_ROI5_2(:,k))/4):3*round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_3(1,k) = mean(aux_VBSD_ROI6_2(2*round(numel(aux_VBSD_ROI6_2(:,k))/4):3*round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_3(1,k) = mean(aux_VBSD_ROI7_2(2*round(numel(aux_VBSD_ROI7_2(:,k))/4):3*round(numel(aux_VBSD_ROI7_2(:,k))/4),k));
                VBSD_ROI8_act_3(1,k) = mean(aux_VBSD_ROI8_2(2*round(numel(aux_VBSD_ROI8_2(:,k))/4):3*round(numel(aux_VBSD_ROI8_2(:,k))/4),k));


                VBSD_ROI5_act_4(1,k) = mean(aux_VBSD_ROI5_2(3*round(numel(aux_VBSD_ROI5_2(:,k))/4):4*(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_4(1,k) = mean(aux_VBSD_ROI6_2(3*round(numel(aux_VBSD_ROI6_2(:,k))/4):4*(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_4(1,k) = mean(aux_VBSD_ROI7_2(3*round(numel(aux_VBSD_ROI7_2(:,k))/4):4*(numel(aux_VBSD_ROI7_2(:,k))/4),k));
                VBSD_ROI8_act_4(1,k) = mean(aux_VBSD_ROI8_2(3*round(numel(aux_VBSD_ROI8_2(:,k))/4):4*(numel(aux_VBSD_ROI8_2(:,k))/4),k));
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                BSD_ROI5_act_1(1,k) = mean(aux_BSD_ROI5(1:round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_1(1,k) = mean(aux_BSD_ROI6(1:round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_1(1,k) = mean(aux_BSD_ROI7(1:round(numel(aux_BSD_ROI7(:,k))/4),k));
                BSD_ROI8_act_1(1,k) = mean(aux_BSD_ROI8(1:round(numel(aux_BSD_ROI8(:,k))/4),k));


                BSD_ROI5_act_2(1,k) = mean(aux_BSD_ROI5(round(numel(aux_BSD_ROI5(:,k))/4):2*round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_2(1,k) = mean(aux_BSD_ROI6(round(numel(aux_BSD_ROI6(:,k))/4):2*round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_2(1,k) = mean(aux_BSD_ROI7(round(numel(aux_BSD_ROI7(:,k))/4):2*round(numel(aux_BSD_ROI7(:,k))/4),k));
                BSD_ROI8_act_2(1,k) = mean(aux_BSD_ROI8(round(numel(aux_BSD_ROI8(:,k))/4):2*round(numel(aux_BSD_ROI8(:,k))/4),k));


                BSD_ROI5_act_3(1,k) = mean(aux_BSD_ROI5(2*round(numel(aux_BSD_ROI5(:,k))/4):3*round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_3(1,k) = mean(aux_BSD_ROI6(2*round(numel(aux_BSD_ROI6(:,k))/4):3*round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_3(1,k) = mean(aux_BSD_ROI7(2*round(numel(aux_BSD_ROI7(:,k))/4):3*round(numel(aux_BSD_ROI7(:,k))/4),k));
                BSD_ROI8_act_3(1,k) = mean(aux_BSD_ROI8(2*round(numel(aux_BSD_ROI8(:,k))/4):3*round(numel(aux_BSD_ROI8(:,k))/4),k));


                BSD_ROI5_act_4(1,k) = mean(aux_BSD_ROI5(3*round(numel(aux_BSD_ROI5(:,k))/4):4*(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_4(1,k) = mean(aux_BSD_ROI6(3*round(numel(aux_BSD_ROI6(:,k))/4):4*(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_4(1,k) = mean(aux_BSD_ROI7(3*round(numel(aux_BSD_ROI7(:,k))/4):4*(numel(aux_BSD_ROI7(:,k))/4),k));
                BSD_ROI8_act_4(1,k) = mean(aux_BSD_ROI8(3*round(numel(aux_BSD_ROI8(:,k))/4):4*(numel(aux_BSD_ROI8(:,k))/4),k));

        end

    %% store in VBSD_ROI(1-8).oldROI(1-5) type of struct
    %% which is used later to calculate ensenble averages (many runs)  
            VBSD_avg_ROI5(ROI5_idx).ROI5 = VBSD_ROI5_avg_2(:)';   
            VBSD_avg_ROI6(ROI5_idx).ROI5 = VBSD_ROI6_avg_2(:)';   
            VBSD_avg_ROI7(ROI5_idx).ROI5 = VBSD_ROI7_avg_2(:)';  
            VBSD_avg_ROI8(ROI5_idx).ROI5 = VBSD_ROI8_avg_2(:)';   

            VBSD_act_1_ROI5(ROI5_idx).ROI5 = VBSD_ROI5_act_1(:)';   
            VBSD_act_1_ROI6(ROI5_idx).ROI5 = VBSD_ROI6_act_1(:)';   
            VBSD_act_1_ROI7(ROI5_idx).ROI5 = VBSD_ROI7_act_1(:)'; 
            VBSD_act_1_ROI8(ROI5_idx).ROI5 = VBSD_ROI8_act_1(:)';  

            VBSD_act_2_ROI5(ROI5_idx).ROI5 = VBSD_ROI5_act_2(:)';   
            VBSD_act_2_ROI6(ROI5_idx).ROI5 = VBSD_ROI6_act_2(:)';   
            VBSD_act_2_ROI7(ROI5_idx).ROI5 = VBSD_ROI7_act_2(:)'; 
            VBSD_act_2_ROI8(ROI5_idx).ROI5 = VBSD_ROI8_act_2(:)';  

            VBSD_act_3_ROI5(ROI5_idx).ROI5 = VBSD_ROI5_act_3(:)';   
            VBSD_act_3_ROI6(ROI5_idx).ROI5 = VBSD_ROI6_act_3(:)';   
            VBSD_act_3_ROI7(ROI5_idx).ROI5 = VBSD_ROI7_act_3(:)'; 
            VBSD_act_3_ROI8(ROI5_idx).ROI5 = VBSD_ROI8_act_3(:)';  

            VBSD_act_4_ROI5(ROI5_idx).ROI5 = VBSD_ROI5_act_4(:)';   
            VBSD_act_4_ROI6(ROI5_idx).ROI5 = VBSD_ROI6_act_4(:)';   
            VBSD_act_4_ROI7(ROI5_idx).ROI5 = VBSD_ROI7_act_4(:)'; 
            VBSD_act_4_ROI8(ROI5_idx).ROI5 = VBSD_ROI8_act_4(:)';  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BSD_act_1_ROI5(ROI5_idx).ROI5 = BSD_ROI5_act_1(:)';   
            BSD_act_1_ROI6(ROI5_idx).ROI5 = BSD_ROI6_act_1(:)';   
            BSD_act_1_ROI7(ROI5_idx).ROI5 = BSD_ROI7_act_1(:)'; 
            BSD_act_1_ROI8(ROI5_idx).ROI5 = BSD_ROI8_act_1(:)';  

            BSD_act_2_ROI5(ROI5_idx).ROI5 = BSD_ROI5_act_2(:)';   
            BSD_act_2_ROI6(ROI5_idx).ROI5 = BSD_ROI6_act_2(:)';   
            BSD_act_2_ROI7(ROI5_idx).ROI5 = BSD_ROI7_act_2(:)'; 
            BSD_act_2_ROI8(ROI5_idx).ROI5 = BSD_ROI8_act_2(:)';  

            BSD_act_3_ROI5(ROI5_idx).ROI5 = BSD_ROI5_act_3(:)';   
            BSD_act_3_ROI6(ROI5_idx).ROI5 = BSD_ROI6_act_3(:)';   
            BSD_act_3_ROI7(ROI5_idx).ROI5 = BSD_ROI7_act_3(:)'; 
            BSD_act_3_ROI8(ROI5_idx).ROI5 = BSD_ROI8_act_3(:)';  

            BSD_act_4_ROI5(ROI5_idx).ROI5 = BSD_ROI5_act_4(:)';   
            BSD_act_4_ROI6(ROI5_idx).ROI5 = BSD_ROI6_act_4(:)';   
            BSD_act_4_ROI7(ROI5_idx).ROI5 = BSD_ROI7_act_4(:)'; 
            BSD_act_4_ROI8(ROI5_idx).ROI5 = BSD_ROI8_act_4(:)';  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% store results in matrices for further analysis 
        if exist('aux_BSD_ROI5_2','var') == 1
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
            ROI5_mtx_idx = ROI5_mtx_idx + 1;  % IMPORTANT!!
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
            for k = 1:numel(BSDX_median)
    %--------------------------------------------------------------------------
    %% place time averaged BSD/VBSD in matrix to use for ensemble avrgs 
    %% each ROI has different number of rows due to the way oldROI overlap and give more or less exp. repetitions.
    %--------------------------------------------------------------------------            
                BSD_avg_ROI5_mtx(ROI5_mtx_idx,k) = aux_BSD_ROI5_2(1,k);
                VBSD_avg_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_avg_2(1,k);

                VBSD_act_1_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_1(1,k);
                VBSD_act_2_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_2(1,k);
                VBSD_act_3_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_3(1,k);
                VBSD_act_4_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_4(1,k);

                BSD_act_1_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_1(1,k);
                BSD_act_2_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_2(1,k);
                BSD_act_3_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_3(1,k);
                BSD_act_4_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_4(1,k);
            end
            for i = NUM_in : NUM
                alpha_ins_ROI5(ROI5_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI5{i-NUM_in+1};
            end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI5_mtx(ROI5_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI5_2(i-NUM_in+1,k); 
                end
            end
        end

        if exist('aux_BSD_ROI6_2','var') == 1
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
            ROI6_mtx_idx = ROI6_mtx_idx + 1;  % IMPORTANT!!
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
            for k = 1:numel(BSDX_median)
    % place BSD in matrix to use for ensemble avrgs
                BSD_avg_ROI6_mtx(ROI6_mtx_idx,k) = aux_BSD_ROI6_2(1,k);
                VBSD_avg_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_avg_2(1,k);

                VBSD_act_1_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_1(1,k);
                VBSD_act_2_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_2(1,k);
                VBSD_act_3_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_3(1,k);
                VBSD_act_4_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_4(1,k);

                BSD_act_1_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_1(1,k);
                BSD_act_2_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_2(1,k);
                BSD_act_3_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_3(1,k);
                BSD_act_4_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_4(1,k);
            end
            for i = NUM_in : NUM
                alpha_ins_ROI6(ROI6_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI6{i-NUM_in+1};
            end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI6_mtx(ROI6_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI6_2(i-NUM_in+1,k); 
                end
            end
        end
        if exist('aux_BSD_ROI7_2','var') == 1
            ROI7_mtx_idx = ROI7_mtx_idx + 1;  % IMPORTANT!!
            for k = 1:numel(BSDX_median)
    % place BSD in matrix to use for ensemble avrgs
                BSD_avg_ROI7_mtx(ROI7_mtx_idx,k) = aux_BSD_ROI7_2(1,k);
                VBSD_avg_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_avg_2(1,k);

                VBSD_act_1_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_1(1,k);
                VBSD_act_2_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_2(1,k);
                VBSD_act_3_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_3(1,k);
                VBSD_act_4_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_4(1,k);

                BSD_act_1_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_1(1,k);
                BSD_act_2_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_2(1,k);
                BSD_act_3_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_3(1,k);
                BSD_act_4_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_4(1,k);
            end
            for i = NUM_in : NUM
                alpha_ins_ROI7(ROI7_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI7{i-NUM_in+1};
            end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI7_mtx(ROI7_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI7_2(i-NUM_in+1,k); 
                end
            end
        end
        if exist('aux_BSD_ROI8_2','var') == 1
            ROI8_mtx_idx = ROI8_mtx_idx + 1;  % IMPORTANT!!
            for k = 1:numel(BSDX_median)
                BSD_avg_ROI8_mtx(ROI8_mtx_idx,k) = aux_BSD_ROI8_2(1,k);
                VBSD_avg_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_avg_2(1,k);
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                VBSD_act_1_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_act_1(1,k);
                VBSD_act_2_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_act_2(1,k);
                VBSD_act_3_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_act_3(1,k);
                VBSD_act_4_ROI8_mtx(ROI8_mtx_idx,k) = VBSD_ROI8_act_4(1,k);

                BSD_act_1_ROI8_mtx(ROI8_mtx_idx,k) = BSD_ROI8_act_1(1,k);
                BSD_act_2_ROI8_mtx(ROI8_mtx_idx,k) = BSD_ROI8_act_2(1,k);
                BSD_act_3_ROI8_mtx(ROI8_mtx_idx,k) = BSD_ROI8_act_3(1,k);
                BSD_act_4_ROI8_mtx(ROI8_mtx_idx,k) = BSD_ROI8_act_4(1,k);
            end
            for i = NUM_in : NUM
                alpha_ins_ROI8(ROI8_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI8{i-NUM_in+1};
            end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI8_mtx(ROI8_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI8_2(i-NUM_in+1,k); 
                end
            end
        end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % save results 
            cd(dir_save)
            save('BSDX_median.mat','BSDX_median','-v7.3');
            save(['BSD_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI5','-v7.3');   
            save(['BSD_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI6','-v7.3');   
            save(['BSD_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI7','-v7.3'); 
            save(['BSD_ins_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI8','-v7.3');  
    % save results
            save(['VBSD_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI5_2','-v7.3');   
            save(['VBSD_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI6_2','-v7.3');   
            save(['VBSD_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI7_2','-v7.3'); 
            save(['VBSD_ins_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI8_2','-v7.3'); 

    %     save(['VBSD_avg_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_avg_2','-v7.3');   
    % 
    %     save(['VBSD_act_1_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_1','-v7.3');
    %     save(['VBSD_act_1_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_act_1','-v7.3');   
    %     
    %     save(['VBSD_act_2_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_2','-v7.3');
    %     save(['VBSD_act_2_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_act_2','-v7.3');   
    %     
    %     save(['VBSD_act_3_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_3','-v7.3');
    %     save(['VBSD_act_3_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_act_3','-v7.3');   
    %         
    %     save(['VBSD_act_4_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI8_act_4','-v7.3');   

            save(['alpha_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI5','-v7.3');   
            save(['alpha_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI6','-v7.3');   
            save(['alpha_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI7','-v7.3');   
            save(['alpha_ins_ROI8_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI8','-v7.3');   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(ROI{case_idx},'ROI4')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear bsd_sums aux_BSD_ROI4_2 aux_BSD_ROI5_2 aux_BSD_ROI6_2 aux_BSD_ROI7_2 centers_ROI4 radii_ROI4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = NUM_in:NUM
    % single valued functions per unit length for ROI4/4
                for k = 2:round(nbin_radii_hg)+1
                    count_ROI4 = 0; % no bubbles initially in each bin/each image
                    count_ROI5 = 0;
                    count_ROI6 = 0;
                    count_ROI7 = 0;
                    for j = 1:numel(centers_xy.(im_name{i-NUM_in+1})(:,2))
                        if centers_xy.(im_name{i-NUM_in+1})(j,2) <= 4.*image_length_mm/4  
    %new:ROI4
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI4 = count_ROI4 + 1;
                                BSD_ROI4{i-NUM_in+1}(1,k-1) =  count_ROI4;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI4{i-NUM_in+1}(1,k-1) =  count_ROI4;  % add 0 if it doesn't find radius for the kth bin
                            end

                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 4.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 5.*image_length_mm/4
    %new:ROI5
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) 
                                count_ROI5 = count_ROI5 + 1;
                                BSD_ROI5{i-NUM_in+1}(1,k-1) =  count_ROI5; 
                            else
                                BSD_ROI5{i-NUM_in+1}(1,k-1) =  count_ROI5; 
                            end

                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 5.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 6.*image_length_mm/4
    %new:ROI6
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) 
                                count_ROI6 = count_ROI6 + 1;
                                BSD_ROI6{i-NUM_in+1}(1,k-1) =  count_ROI6;  
                            else
                                BSD_ROI6{i-NUM_in+1}(1,k-1) =  count_ROI6;  
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 6.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 7.*image_length_mm/4
    %new:ROI7
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k)
                                count_ROI7 = count_ROI7 + 1;
                                BSD_ROI7{i-NUM_in+1}(1,k-1) =  count_ROI7;  
                            else
                                BSD_ROI7{i-NUM_in+1}(1,k-1) =  count_ROI7;  
                            end
                        end
                    end
                end
            end
            for i=NUM_in:NUM
                    if isempty(BSD_ROI4{i-NUM_in+1}) 
                        BSD_ROI4{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI5{i-NUM_in+1}) 
                        BSD_ROI5{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI6{i-NUM_in+1}) 
                        BSD_ROI6{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI7{i-NUM_in+1}) 
                        BSD_ROI7{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
            end
            clear radii_ROI4 centers_ROI4
            for i=NUM_in:NUM
%% for methodology description, see details in oldROI5 part.
    % VBSD per unit length
    % 1. sort radii for each new ROI/each image according to the BSD_ROI4-7 findings
                bsd_sums = [1,sum(BSD_ROI4{1,i-NUM_in+1}),sum(BSD_ROI4{1,i-NUM_in+1}+BSD_ROI5{1,i-NUM_in+1}),...
                sum(BSD_ROI4{1,i-NUM_in+1}+BSD_ROI5{1,i-NUM_in+1}+BSD_ROI6{1,i-NUM_in+1}),...
                sum(BSD_ROI4{1,i-NUM_in+1}+BSD_ROI5{1,i-NUM_in+1}+BSD_ROI6{1,i-NUM_in+1}+BSD_ROI7{1,i-NUM_in+1})]; % indices to calculate Volume for each ROI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1 find total bubble numbers per BSD of new ROI
    % bsd_total_sums.(im_name{i-NUM_in+1}) = sum(BSD_ROI4{1,i-NUM_in+1});
    % 1.2 find radii and centers corresponding per BSD in new ROI
    
    
    %% WHAT IS THAT??? IS IT FOR EUGENY??
                radii_ROI4.(im_name{i-NUM_in+1})  = radii_mm_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2));
                centers_ROI4.(im_name{i-NUM_in+1})(:,1)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2),1);
                centers_ROI4.(im_name{i-NUM_in+1})(:,2)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2),2);
                centers_ROI4.(im_name{i-NUM_in+1})       = round(centers_ROI4.(im_name{i-NUM_in+1}),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort 
                vol_b_xy_ROI4.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2)));
                vol_b_xy_ROI5.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3)));
                vol_b_xy_ROI6.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4)));
                vol_b_xy_ROI7.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5)));
    % 2. calculate VBSD for each new ROI/each bin/each image 
                count_vol_ROI4{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI4
                count_vol_ROI5{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI5
                count_vol_ROI6{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI6
                count_vol_ROI7{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI7
                for k = 2:numel(BSDX_median)+1
    % calculate volume as sum of volumes of bubbles of each bin found in BSD_ROI
                    count_vol_ROI4{i-NUM_in+1}(1,k)  = count_vol_ROI4{i-NUM_in+1}(1,k-1) + BSD_ROI4{i-NUM_in+1}(1,k-1);  % find indices of bubbles that are in same bins in BSD 
                    aux_VBSD_ROI4{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI4.(im_name{i-NUM_in+1})...
                        (count_vol_ROI4{i-NUM_in+1}(1,k-1):count_vol_ROI4{i-NUM_in+1}(1,k)-1)); % calculate the sum of volumes of bubbles in same bin

                    count_vol_ROI5{i-NUM_in+1}(1,k)  = count_vol_ROI5{i-NUM_in+1}(1,k-1) + BSD_ROI5{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI5{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI5.(im_name{i-NUM_in+1})...
                        (count_vol_ROI5{i-NUM_in+1}(1,k-1):count_vol_ROI5{i-NUM_in+1}(1,k)-1)); 

                    count_vol_ROI6{i-NUM_in+1}(1,k)  = count_vol_ROI6{i-NUM_in+1}(1,k-1) + BSD_ROI6{i-NUM_in+1}(1,k-1);  
                    aux_VBSD_ROI6{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI6.(im_name{i-NUM_in+1})...
                        (count_vol_ROI6{i-NUM_in+1}(1,k-1):count_vol_ROI6{i-NUM_in+1}(1,k)-1)); 

                    count_vol_ROI7{i-NUM_in+1}(1,k)  = count_vol_ROI7{i-NUM_in+1}(1,k-1) + BSD_ROI7{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI7{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI7.(im_name{i-NUM_in+1})...
                        (count_vol_ROI7{i-NUM_in+1}(1,k-1):count_vol_ROI7{i-NUM_in+1}(1,k)-1)); 
                end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate instanteneous void fraction for ROI1-8
                aux_alpha_ins_ROI4{i-NUM_in+1} = sum(aux_VBSD_ROI4{i-NUM_in+1});
                aux_alpha_ins_ROI5{i-NUM_in+1} = sum(aux_VBSD_ROI5{i-NUM_in+1});
                aux_alpha_ins_ROI6{i-NUM_in+1} = sum(aux_VBSD_ROI6{i-NUM_in+1});
                aux_alpha_ins_ROI7{i-NUM_in+1} = sum(aux_VBSD_ROI7{i-NUM_in+1});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
    % % save radii and centers to create txt files     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % indices of images to use in order to aling videos in time (see bubble_time_analysis.m for more info) 
            if strcmp(phase{case_idx},'peak')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'trough')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'posp2')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'minp2')
                   time_idx_in = 1;
                   time_idx_end = NUM;
            end
            % align instant data in time 
            for i = time_idx_in: time_idx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                radii_ins_ROI4.(im_name{i+1-time_idx_in}) = radii_ROI4.(im_name{i});
                centers_ins_ROI4.(im_name{i+1-time_idx_in}) = centers_ROI4.(im_name{i});
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cd('C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\GW_175\191020\server\entire_video\ROI1-8\text files');
            for i = 1:(time_idx_end-time_idx_in)
            fid=fopen(['centers_ins_ROI4.',(im_name{i}),'.old',ROI{case_idx},'.',run{case_idx},'.txt'],'wt');
            header1 = 'y(mm)';
            header2 = 'x(mm)';
            fprintf(fid, [ ' ' header1 ' ' header2 '\r\n']);
            fprintf(fid, '%4.1f %4.1f\r\n',centers_ins_ROI4.(im_name{i})');
            fclose(fid);
            end
            for i = 1:(time_idx_end-time_idx_in)
            fid=fopen(['radii_ins_ROI4.',(im_name{i}),'.old',ROI{case_idx},'.',run{case_idx},'.txt'],'wt');
            header1 = 'r(mm)';
            fprintf(fid, [ ' ' header1 '\r\n']);
            fprintf(fid, '%4.2f\r\n',radii_ins_ROI4.(im_name{i}));
            fclose(fid);
            end
            save(['radii_ins_ROI4_',phase{case_idx},'.old',ROI{case_idx},'.',run{case_idx},'.mat'],'-struct','radii_ins_ROI4');
            save(['centers_ins_ROI4_',phase{case_idx},'.old',ROI{case_idx},'.',run{case_idx},'.mat'],'-struct','centers_ins_ROI4');
            cd(dir_save)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fid=fopen('BSD_ins_ROI4.txt','wt');
            header = 'radius bin -->';
            formatspec = '5.0f';
            fprintf(fid, [ header '\r\n']);
            fprintf(fid, '%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\r\n', BSDX_median);
            fprintf(fid, '%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\r\n', aux_BSD_ROI5');
            fclose(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     save(['BSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4_2','-v7.3');   % name structure of images (e.g. ...peak_0001)
    %     save(['BSD_avg_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI5_2','-v7.3');   
    %     save(['BSD_avg_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI6_2','-v7.3');   
    %     save(['BSD_avg_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI7_2','-v7.3');  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% time averaged BSD 
            for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_BSD_ROI4(i-NUM_in+1,k) = (BSD_ROI4{1,i-NUM_in+1}(k));  % place BSD values of each image/rep. in a 2D matrix 
                aux_BSD_ROI5(i-NUM_in+1,k) = (BSD_ROI5{1,i-NUM_in+1}(k));  
                aux_BSD_ROI6(i-NUM_in+1,k) = (BSD_ROI6{1,i-NUM_in+1}(k)); 
                aux_BSD_ROI7(i-NUM_in+1,k) = (BSD_ROI7{1,i-NUM_in+1}(k));  
            end
            end
            for k = 1:numel(BSDX_median)
    % find mean per radius bin for each case and place them in cell matrix
                aux_BSD_ROI4_2(1,k) = mean(aux_BSD_ROI4(:,k));
                aux_BSD_ROI5_2(1,k) = mean(aux_BSD_ROI5(:,k));
                aux_BSD_ROI6_2(1,k) = mean(aux_BSD_ROI6(:,k));
                aux_BSD_ROI7_2(1,k) = mean(aux_BSD_ROI7(:,k));
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ROI4_idx = ROI4_idx + 1; % new ROI4 index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store in BSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages 
            BSD_avg_ROI4(ROI4_idx).ROI4 = aux_BSD_ROI4_2(:)';  
            BSD_avg_ROI5(ROI4_idx).ROI4 = aux_BSD_ROI5_2(:)';   
            BSD_avg_ROI6(ROI4_idx).ROI4 = aux_BSD_ROI6_2(:)';   
            BSD_avg_ROI7(ROI4_idx).ROI4 = aux_BSD_ROI7_2(:)'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged VBSD
            for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_VBSD_ROI4_2(i-NUM_in+1,k) = (aux_VBSD_ROI4{1,i-NUM_in+1}(k));  % place VBSD values of each image/rep. in a 2D matrix 
                aux_VBSD_ROI5_2(i-NUM_in+1,k) = (aux_VBSD_ROI5{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI6_2(i-NUM_in+1,k) = (aux_VBSD_ROI6{1,i-NUM_in+1}(k)); 
                aux_VBSD_ROI7_2(i-NUM_in+1,k) = (aux_VBSD_ROI7{1,i-NUM_in+1}(k));  
            end
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:numel(BSDX_median)
    % find mean per radius bin for each case and place them in cell matrix
                VBSD_ROI4_avg_2(1,k) = mean(aux_VBSD_ROI4_2(:,k));
                VBSD_ROI5_avg_2(1,k) = mean(aux_VBSD_ROI5_2(:,k));
                VBSD_ROI6_avg_2(1,k) = mean(aux_VBSD_ROI6_2(:,k));
                VBSD_ROI7_avg_2(1,k) = mean(aux_VBSD_ROI7_2(:,k));
    % find time averages for specific time slots of videos
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% aux_VBSD_2 is equivalent to aux_BSD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                VBSD_ROI4_act_1(1,k) = mean(aux_VBSD_ROI4_2(1:round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_1(1,k) = mean(aux_VBSD_ROI5_2(1:round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_1(1,k) = mean(aux_VBSD_ROI6_2(1:round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_1(1,k) = mean(aux_VBSD_ROI7_2(1:round(numel(aux_VBSD_ROI7_2(:,k))/4),k));

                VBSD_ROI4_act_2(1,k) = mean(aux_VBSD_ROI4_2(round(numel(aux_VBSD_ROI4_2(:,k))/4):2*round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_2(1,k) = mean(aux_VBSD_ROI5_2(round(numel(aux_VBSD_ROI5_2(:,k))/4):2*round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_2(1,k) = mean(aux_VBSD_ROI6_2(round(numel(aux_VBSD_ROI6_2(:,k))/4):2*round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_2(1,k) = mean(aux_VBSD_ROI7_2(round(numel(aux_VBSD_ROI7_2(:,k))/4):2*round(numel(aux_VBSD_ROI7_2(:,k))/4),k));

                VBSD_ROI4_act_3(1,k) = mean(aux_VBSD_ROI4_2(2*round(numel(aux_VBSD_ROI4_2(:,k))/4):3*round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_3(1,k) = mean(aux_VBSD_ROI5_2(2*round(numel(aux_VBSD_ROI5_2(:,k))/4):3*round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_3(1,k) = mean(aux_VBSD_ROI6_2(2*round(numel(aux_VBSD_ROI6_2(:,k))/4):3*round(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_3(1,k) = mean(aux_VBSD_ROI7_2(2*round(numel(aux_VBSD_ROI7_2(:,k))/4):3*round(numel(aux_VBSD_ROI7_2(:,k))/4),k));

                VBSD_ROI4_act_4(1,k) = mean(aux_VBSD_ROI4_2(3*round(numel(aux_VBSD_ROI4_2(:,k))/4):4*(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_4(1,k) = mean(aux_VBSD_ROI5_2(3*round(numel(aux_VBSD_ROI5_2(:,k))/4):4*(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_4(1,k) = mean(aux_VBSD_ROI6_2(3*round(numel(aux_VBSD_ROI6_2(:,k))/4):4*(numel(aux_VBSD_ROI6_2(:,k))/4),k));
                VBSD_ROI7_act_4(1,k) = mean(aux_VBSD_ROI7_2(3*round(numel(aux_VBSD_ROI7_2(:,k))/4):4*(numel(aux_VBSD_ROI7_2(:,k))/4),k));
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                BSD_ROI4_act_1(1,k) = mean(aux_BSD_ROI4(1:round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_1(1,k) = mean(aux_BSD_ROI5(1:round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_1(1,k) = mean(aux_BSD_ROI6(1:round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_1(1,k) = mean(aux_BSD_ROI7(1:round(numel(aux_BSD_ROI7(:,k))/4),k));


                BSD_ROI4_act_2(1,k) = mean(aux_BSD_ROI4(round(numel(aux_BSD_ROI4(:,k))/4):2*round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_2(1,k) = mean(aux_BSD_ROI5(round(numel(aux_BSD_ROI5(:,k))/4):2*round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_2(1,k) = mean(aux_BSD_ROI6(round(numel(aux_BSD_ROI6(:,k))/4):2*round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_2(1,k) = mean(aux_BSD_ROI7(round(numel(aux_BSD_ROI7(:,k))/4):2*round(numel(aux_BSD_ROI7(:,k))/4),k));


                BSD_ROI4_act_3(1,k) = mean(aux_BSD_ROI4(2*round(numel(aux_BSD_ROI4(:,k))/4):3*round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_3(1,k) = mean(aux_BSD_ROI5(2*round(numel(aux_BSD_ROI5(:,k))/4):3*round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_3(1,k) = mean(aux_BSD_ROI6(2*round(numel(aux_BSD_ROI6(:,k))/4):3*round(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_3(1,k) = mean(aux_BSD_ROI7(2*round(numel(aux_BSD_ROI7(:,k))/4):3*round(numel(aux_BSD_ROI7(:,k))/4),k));


                BSD_ROI4_act_4(1,k) = mean(aux_BSD_ROI4(3*round(numel(aux_BSD_ROI4(:,k))/4):4*(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_4(1,k) = mean(aux_BSD_ROI5(3*round(numel(aux_BSD_ROI5(:,k))/4):4*(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_4(1,k) = mean(aux_BSD_ROI6(3*round(numel(aux_BSD_ROI6(:,k))/4):4*(numel(aux_BSD_ROI6(:,k))/4),k));
                BSD_ROI7_act_4(1,k) = mean(aux_BSD_ROI7(3*round(numel(aux_BSD_ROI7(:,k))/4):4*(numel(aux_BSD_ROI7(:,k))/4),k));
            end

    % store in VBSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages     
            VBSD_avg_ROI4(ROI4_idx).ROI4 = VBSD_ROI4_avg_2(:)';  
            VBSD_avg_ROI5(ROI4_idx).ROI4 = VBSD_ROI5_avg_2(:)';   
            VBSD_avg_ROI6(ROI4_idx).ROI4 = VBSD_ROI6_avg_2(:)';   
            VBSD_avg_ROI7(ROI4_idx).ROI4 = VBSD_ROI7_avg_2(:)';  

            VBSD_act_1_ROI4(ROI4_idx).ROI4 = VBSD_ROI4_act_1(:)';  
            VBSD_act_1_ROI5(ROI4_idx).ROI4 = VBSD_ROI5_act_1(:)';   
            VBSD_act_1_ROI6(ROI4_idx).ROI4 = VBSD_ROI6_act_1(:)';   
            VBSD_act_1_ROI7(ROI4_idx).ROI4 = VBSD_ROI7_act_1(:)'; 

            VBSD_act_2_ROI4(ROI4_idx).ROI4 = VBSD_ROI4_act_2(:)';  
            VBSD_act_2_ROI5(ROI4_idx).ROI4 = VBSD_ROI5_act_2(:)';   
            VBSD_act_2_ROI6(ROI4_idx).ROI4 = VBSD_ROI6_act_2(:)';   
            VBSD_act_2_ROI7(ROI4_idx).ROI4 = VBSD_ROI7_act_2(:)'; 

            VBSD_act_3_ROI4(ROI4_idx).ROI4 = VBSD_ROI4_act_3(:)';  
            VBSD_act_3_ROI5(ROI4_idx).ROI4 = VBSD_ROI5_act_3(:)';   
            VBSD_act_3_ROI6(ROI4_idx).ROI4 = VBSD_ROI6_act_3(:)';   
            VBSD_act_3_ROI7(ROI4_idx).ROI4 = VBSD_ROI7_act_3(:)'; 

            VBSD_act_4_ROI4(ROI4_idx).ROI4 = VBSD_ROI4_act_4(:)';  
            VBSD_act_4_ROI5(ROI4_idx).ROI4 = VBSD_ROI5_act_4(:)';   
            VBSD_act_4_ROI6(ROI4_idx).ROI4 = VBSD_ROI6_act_4(:)';   
            VBSD_act_4_ROI7(ROI4_idx).ROI4 = VBSD_ROI7_act_4(:)'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            BSD_act_1_ROI4(ROI4_idx).ROI4 = BSD_ROI4_act_1(:)';  
            BSD_act_1_ROI5(ROI4_idx).ROI4 = BSD_ROI5_act_1(:)';   
            BSD_act_1_ROI6(ROI4_idx).ROI4 = BSD_ROI6_act_1(:)';   
            BSD_act_1_ROI7(ROI4_idx).ROI4 = BSD_ROI7_act_1(:)'; 

            BSD_act_2_ROI4(ROI4_idx).ROI4 = BSD_ROI4_act_2(:)';  
            BSD_act_2_ROI5(ROI4_idx).ROI4 = BSD_ROI5_act_2(:)';   
            BSD_act_2_ROI6(ROI4_idx).ROI4 = BSD_ROI6_act_2(:)';   
            BSD_act_2_ROI7(ROI4_idx).ROI4 = BSD_ROI7_act_2(:)'; 

            BSD_act_3_ROI4(ROI4_idx).ROI4 = BSD_ROI4_act_3(:)';  
            BSD_act_3_ROI5(ROI4_idx).ROI4 = BSD_ROI5_act_3(:)';   
            BSD_act_3_ROI6(ROI4_idx).ROI4 = BSD_ROI6_act_3(:)';   
            BSD_act_3_ROI7(ROI4_idx).ROI4 = BSD_ROI7_act_3(:)'; 

            BSD_act_4_ROI4(ROI4_idx).ROI4 = BSD_ROI4_act_4(:)';  
            BSD_act_4_ROI5(ROI4_idx).ROI4 = BSD_ROI5_act_4(:)';   
            BSD_act_4_ROI6(ROI4_idx).ROI4 = BSD_ROI6_act_4(:)';   
            BSD_act_4_ROI7(ROI4_idx).ROI4 = BSD_ROI7_act_4(:)'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store results in matrices for further analysis 
            if exist('aux_BSD_ROI4_2','var') == 1
                ROI4_mtx_idx = ROI4_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI4_mtx(ROI4_mtx_idx,k) = aux_BSD_ROI4_2(1,k);
                    VBSD_avg_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_avg_2(1,k);
                    VBSD_act_1_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_1(1,k);
                    VBSD_act_2_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_2(1,k);
                    VBSD_act_3_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_3(1,k);
                    VBSD_act_4_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_4(1,k);

                    BSD_act_1_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_1(1,k);
                    BSD_act_2_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_2(1,k);
                    BSD_act_3_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_3(1,k);
                    BSD_act_4_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI4(ROI4_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI4{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI4_mtx(ROI4_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI4_2(i-NUM_in+1,k); 
                end
            end
            end  
            if exist('aux_BSD_ROI5_2','var') == 1
                ROI5_mtx_idx = ROI5_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI5_mtx(ROI5_mtx_idx,k) = aux_BSD_ROI5_2(1,k);
                    
                    VBSD_avg_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_avg_2(1,k);
                    VBSD_act_1_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_1(1,k);
                    VBSD_act_2_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_2(1,k);
                    VBSD_act_3_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_3(1,k);
                    VBSD_act_4_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_4(1,k);

                    BSD_act_1_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_1(1,k);
                    BSD_act_2_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_2(1,k);
                    BSD_act_3_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_3(1,k);
                    BSD_act_4_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI5(ROI5_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI5{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI5_mtx(ROI5_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI5_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI6_2','var') == 1
                ROI6_mtx_idx = ROI6_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI6_mtx(ROI6_mtx_idx,k) = aux_BSD_ROI6_2(1,k);
                    
                    VBSD_avg_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_avg_2(1,k);
                    VBSD_act_1_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_1(1,k);
                    VBSD_act_2_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_2(1,k);
                    VBSD_act_3_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_3(1,k);
                    VBSD_act_4_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_4(1,k);

                    BSD_act_1_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_1(1,k);
                    BSD_act_2_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_2(1,k);
                    BSD_act_3_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_3(1,k);
                    BSD_act_4_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI6(ROI6_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI6{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI6_mtx(ROI6_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI6_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI7_2','var') == 1
                ROI7_mtx_idx = ROI7_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI7_mtx(ROI7_mtx_idx,k) = aux_BSD_ROI7_2(1,k);
                    
                    VBSD_avg_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_avg_2(1,k);
                    VBSD_act_1_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_1(1,k);
                    VBSD_act_2_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_2(1,k);
                    VBSD_act_3_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_3(1,k);
                    VBSD_act_4_ROI7_mtx(ROI7_mtx_idx,k) = VBSD_ROI7_act_4(1,k);

                    BSD_act_1_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_1(1,k);
                    BSD_act_2_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_2(1,k);
                    BSD_act_3_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_3(1,k);
                    BSD_act_4_ROI7_mtx(ROI7_mtx_idx,k) = BSD_ROI7_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI7(ROI7_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI7{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI7_mtx(ROI7_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI7_2(i-NUM_in+1,k); 
                end
            end
            end
     % save results
            cd(dir_save);
            save('BSDX_median.mat','BSDX_median','-v7.3');
            save(['BSD_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4','-v7.3');  
            save(['BSD_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI5','-v7.3');   
            save(['BSD_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI6','-v7.3');   
            save(['BSD_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI7','-v7.3'); 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save results
            save(['VBSD_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI4_2','-v7.3'); 
            save(['VBSD_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI5_2','-v7.3');   
            save(['VBSD_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI6_2','-v7.3');   
            save(['VBSD_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI7_2','-v7.3'); 

    %     save(['VBSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_avg_2','-v7.3');   % name structure of images (e.g. ...peak_0001)
    %     save(['VBSD_avg_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_avg_2','-v7.3');     
    % 
    %     save(['VBSD_act_1_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_1','-v7.3');
    %     
    %     save(['VBSD_act_2_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_2','-v7.3');
    %     
    %     save(['VBSD_act_3_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_3','-v7.3');
    %         
    %     save(['VBSD_act_4_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI7_act_4','-v7.3');

            save(['alpha_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI4','-v7.3');  
            save(['alpha_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI5','-v7.3');   
            save(['alpha_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI6','-v7.3');   
            save(['alpha_ins_ROI7_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI7','-v7.3');   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(ROI{case_idx},'ROI3')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear bsd_sums aux_BSD_ROI3_2 aux_BSD_ROI4_2 aux_BSD_ROI5_2  aux_BSD_ROI6_2 
            clear centers_ROI4 radii_ROI4 bsd_sums_ROI4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = NUM_in:NUM
%% for methodology description, see details in oldROI5 part.
    % single valued functions per unit length for ROI4/4
                for k = 2:round(nbin_radii_hg)+1
                    count_ROI3 = 0;
                    count_ROI4 = 0; % no bubbles initially in each bin/each image
                    count_ROI5 = 0;
                    count_ROI6 = 0;
                    for j = 1:numel(centers_xy.(im_name{i-NUM_in+1})(:,2))
                        if centers_xy.(im_name{i-NUM_in+1})(j,2)  <= 3.*image_length_mm/4 
    %new:ROI3
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k)
                                count_ROI3 = count_ROI3 + 1;
                                BSD_ROI3{i-NUM_in+1}(1,k-1) =  count_ROI3;  
                            else
                                BSD_ROI3{i-NUM_in+1}(1,k-1) =  count_ROI3;  
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 3.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 4.*image_length_mm/4
    %new:ROI4
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI4 = count_ROI4 + 1;
                                BSD_ROI4{i-NUM_in+1}(1,k-1) =  count_ROI4;  % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI4{i-NUM_in+1}(1,k-1) =  count_ROI4;  % add 0 if it doesn't find radius for the kth bin
                            end

                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 4.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 5.*image_length_mm/4
    %new:ROI5
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) 
                                count_ROI5 = count_ROI5 + 1;
                                BSD_ROI5{i-NUM_in+1}(1,k-1) =  count_ROI5; 
                            else
                                BSD_ROI5{i-NUM_in+1}(1,k-1) =  count_ROI5; 
                            end

                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 5.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 6.*image_length_mm/4
    %new:ROI6
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) 
                                count_ROI6 = count_ROI6 + 1;
                                BSD_ROI6{i-NUM_in+1}(1,k-1) =  count_ROI6;  
                            else
                                BSD_ROI6{i-NUM_in+1}(1,k-1) =  count_ROI6;  
                            end
                        else
                        end
                    end
                end
            end
            for i=NUM_in:NUM
                    if isempty(BSD_ROI3{i-NUM_in+1}) 
                        BSD_ROI3{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI4{i-NUM_in+1}) 
                        BSD_ROI4{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI5{i-NUM_in+1}) 
                        BSD_ROI5{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI6{i-NUM_in+1}) 
                        BSD_ROI6{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
            end
            clear centers_ROI4 radii_ROI4
            for i=NUM_in:NUM
    % VBSD per unit length
%% for methodology description, see details in oldROI5 part.
    % 1. sort radii for each new ROI/each image according to the BSD_ROI4-7 findings
                bsd_sums = [1,sum(BSD_ROI3{1,i-NUM_in+1}),sum(BSD_ROI3{1,i-NUM_in+1}+BSD_ROI4{1,i-NUM_in+1}),...
                sum(BSD_ROI3{1,i-NUM_in+1}+BSD_ROI4{1,i-NUM_in+1}+BSD_ROI5{1,i-NUM_in+1}),...
                sum(BSD_ROI3{1,i-NUM_in+1}+BSD_ROI4{1,i-NUM_in+1}+BSD_ROI5{1,i-NUM_in+1}+BSD_ROI6{1,i-NUM_in+1})]; % indices to calculate Volume for each ROI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1 find radii and centers corresponding per BSD in new ROI : TEST
                radii_ROI4.(im_name{i-NUM_in+1})  = radii_mm_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3));
                centers_ROI4.(im_name{i-NUM_in+1})(:,1)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3),1);
                centers_ROI4.(im_name{i-NUM_in+1})(:,2)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3),2);
                centers_ROI4.(im_name{i-NUM_in+1})       = round(centers_ROI4.(im_name{i-NUM_in+1}),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort 
                vol_b_xy_ROI3.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2)));
                vol_b_xy_ROI4.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3)));
                vol_b_xy_ROI5.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4)));
                vol_b_xy_ROI6.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5)));
    % 2. calculate VBSD for each new ROI/each bin/each image 
                count_vol_ROI3{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI3
                count_vol_ROI4{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI4
                count_vol_ROI5{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI5
                count_vol_ROI6{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI6
                for k = 2:numel(BSDX_median)+1
    % calculate volume as sum of volumes of bubbles of each bin found in BSD_ROI
                    count_vol_ROI3{i-NUM_in+1}(1,k)  = count_vol_ROI3{i-NUM_in+1}(1,k-1) + BSD_ROI3{i-NUM_in+1}(1,k-1); % index of kth bubble bin in ROI
                    aux_VBSD_ROI3{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI3.(im_name{i-NUM_in+1})...  % find number of bubbles per bin in specific ROI
                        (count_vol_ROI3{i-NUM_in+1}(1,k-1):count_vol_ROI3{i-NUM_in+1}(1,k)-1)); 

                    count_vol_ROI4{i-NUM_in+1}(1,k)  = count_vol_ROI4{i-NUM_in+1}(1,k-1) + BSD_ROI4{i-NUM_in+1}(1,k-1);  % find indices of bubbles that are in same bins in BSD 
                    aux_VBSD_ROI4{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI4.(im_name{i-NUM_in+1})...
                        (count_vol_ROI4{i-NUM_in+1}(1,k-1):count_vol_ROI4{i-NUM_in+1}(1,k)-1)); % calculate the sum of volumes of bubbles in same bin

                    count_vol_ROI5{i-NUM_in+1}(1,k)  = count_vol_ROI5{i-NUM_in+1}(1,k-1) + BSD_ROI5{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI5{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI5.(im_name{i-NUM_in+1})...
                        (count_vol_ROI5{i-NUM_in+1}(1,k-1):count_vol_ROI5{i-NUM_in+1}(1,k)-1)); 

                    count_vol_ROI6{i-NUM_in+1}(1,k)  = count_vol_ROI6{i-NUM_in+1}(1,k-1) + BSD_ROI6{i-NUM_in+1}(1,k-1);  
                    aux_VBSD_ROI6{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI6.(im_name{i-NUM_in+1})...
                        (count_vol_ROI6{i-NUM_in+1}(1,k-1):count_vol_ROI6{i-NUM_in+1}(1,k)-1)); 
                end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate instanteneous void fraction for ROI1-8
                aux_alpha_ins_ROI3{i-NUM_in+1} = sum(aux_VBSD_ROI3{i-NUM_in+1});
                aux_alpha_ins_ROI4{i-NUM_in+1} = sum(aux_VBSD_ROI4{i-NUM_in+1});
                aux_alpha_ins_ROI5{i-NUM_in+1} = sum(aux_VBSD_ROI5{i-NUM_in+1});
                aux_alpha_ins_ROI6{i-NUM_in+1} = sum(aux_VBSD_ROI6{i-NUM_in+1});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
    % save txt files   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort indices of images, in order to align videos in time (see bubble_time_analysis.m for more info) 
            if strcmp(phase{case_idx},'peak')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'trough')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'posp2')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'minp2')
                   time_idx_in = 1;
                   time_idx_end = NUM;
            end
            % align instant data in time 
            for i = time_idx_in: time_idx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                radii_ins_ROI4.(im_name{i+1-time_idx_in}) = radii_ROI4.(im_name{i});
                centers_ins_ROI4.(im_name{i+1-time_idx_in}) = centers_ROI4.(im_name{i});
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cd('C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\GW_175\191020\server\entire_video\ROI1-8\text files');
            for i = 1:(time_idx_end-time_idx_in)
            fid=fopen(['centers_ins_ROI4.',(im_name{i}),'.old',ROI{case_idx},'.',run{case_idx},'.txt'],'wt');
            header1 = 'y(mm)';
            header2 = 'x(mm)';
            fprintf(fid, [ ' ' header1 ' ' header2 '\r\n']);
            fprintf(fid, '%4.1f %4.1f\r\n',centers_ROI4.(im_name{i})');
            fclose(fid);
            end
            for i = 1:(time_idx_end-time_idx_in)
            fid=fopen(['radii_ins_ROI4.',(im_name{i}),'.old',ROI{case_idx},'.',run{case_idx},'.txt'],'wt');
            header1 = 'r(mm)';
            fprintf(fid, [ ' ' header1 '\r\n']);
            fprintf(fid, '%4.2f\r\n',radii_ROI4.(im_name{i}));
            fclose(fid);
            end
            save(['radii_ins_ROI4_',phase{case_idx},'.old',ROI{case_idx},'.',run{case_idx},'.mat'],'-struct','radii_ROI4');
            save(['centers_ins_ROI4_',phase{case_idx},'.old',ROI{case_idx},'.',run{case_idx},'.mat'],'-struct','centers_ROI4');
            cd(dir_save)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %     save(['BSD_avg_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI3_2','-v7.3');  
    %     save(['BSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4_2','-v7.3');   % name structure of images (e.g. ...peak_0001)
    %     save(['BSD_avg_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI5_2','-v7.3');   
    %     save(['BSD_avg_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI6_2','-v7.3');   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged BSD 
            for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_BSD_ROI3(i-NUM_in+1,k) = (BSD_ROI3{1,i-NUM_in+1}(k));  
                aux_BSD_ROI4(i-NUM_in+1,k) = (BSD_ROI4{1,i-NUM_in+1}(k));  % place BSD values of each image/rep. in a 2D matrix 
                aux_BSD_ROI5(i-NUM_in+1,k) = (BSD_ROI5{1,i-NUM_in+1}(k));  
                aux_BSD_ROI6(i-NUM_in+1,k) = (BSD_ROI6{1,i-NUM_in+1}(k)); 
            end
            end
            for k = 1:numel(BSDX_median)
    % find mean per radius bin for each case and place them in cell matrix
                aux_BSD_ROI3_2(1,k) = mean(aux_BSD_ROI3(:,k));
                aux_BSD_ROI4_2(1,k) = mean(aux_BSD_ROI4(:,k));
                aux_BSD_ROI5_2(1,k) = mean(aux_BSD_ROI5(:,k));
                aux_BSD_ROI6_2(1,k) = mean(aux_BSD_ROI6(:,k));
            end
            ROI3_idx = ROI3_idx + 1; % new ROI3 index
    % store in BSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages 
            BSD_avg_ROI3(ROI3_idx).ROI3 = aux_BSD_ROI3_2(:)';  
            BSD_avg_ROI4(ROI3_idx).ROI3 = aux_BSD_ROI4_2(:)';  
            BSD_avg_ROI5(ROI3_idx).ROI3 = aux_BSD_ROI5_2(:)';   
            BSD_avg_ROI6(ROI3_idx).ROI3 = aux_BSD_ROI6_2(:)'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged VBSD
            for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_VBSD_ROI3_2(i-NUM_in+1,k) = (aux_VBSD_ROI3{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI4_2(i-NUM_in+1,k) = (aux_VBSD_ROI4{1,i-NUM_in+1}(k));  % place VBSD values of each image/rep. in a 2D matrix 
                aux_VBSD_ROI5_2(i-NUM_in+1,k) = (aux_VBSD_ROI5{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI6_2(i-NUM_in+1,k) = (aux_VBSD_ROI6{1,i-NUM_in+1}(k)); 
            end
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:numel(BSDX_median)
    % find mean per radius bin for each case and place them in cell matrix
                VBSD_ROI3_avg_2(1,k) = mean(aux_VBSD_ROI3_2(:,k));
                VBSD_ROI4_avg_2(1,k) = mean(aux_VBSD_ROI4_2(:,k));
                VBSD_ROI5_avg_2(1,k) = mean(aux_VBSD_ROI5_2(:,k));
                VBSD_ROI6_avg_2(1,k) = mean(aux_VBSD_ROI6_2(:,k));
    % find time averages for segments of videos
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% aux_VBSD_2 is equivalent to aux_BSD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                VBSD_ROI3_act_1(1,k) = mean(aux_VBSD_ROI3_2(1:round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_1(1,k) = mean(aux_VBSD_ROI4_2(1:round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_1(1,k) = mean(aux_VBSD_ROI5_2(1:round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_1(1,k) = mean(aux_VBSD_ROI6_2(1:round(numel(aux_VBSD_ROI6_2(:,k))/4),k));

                VBSD_ROI3_act_2(1,k) = mean(aux_VBSD_ROI3_2(round(numel(aux_VBSD_ROI3_2(:,k))/4):2*round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_2(1,k) = mean(aux_VBSD_ROI4_2(round(numel(aux_VBSD_ROI4_2(:,k))/4):2*round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_2(1,k) = mean(aux_VBSD_ROI5_2(round(numel(aux_VBSD_ROI5_2(:,k))/4):2*round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_2(1,k) = mean(aux_VBSD_ROI6_2(round(numel(aux_VBSD_ROI6_2(:,k))/4):2*round(numel(aux_VBSD_ROI6_2(:,k))/4),k));

                VBSD_ROI3_act_3(1,k) = mean(aux_VBSD_ROI3_2(2*round(numel(aux_VBSD_ROI3_2(:,k))/4):3*round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_3(1,k) = mean(aux_VBSD_ROI4_2(2*round(numel(aux_VBSD_ROI4_2(:,k))/4):3*round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_3(1,k) = mean(aux_VBSD_ROI5_2(2*round(numel(aux_VBSD_ROI5_2(:,k))/4):3*round(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_3(1,k) = mean(aux_VBSD_ROI6_2(2*round(numel(aux_VBSD_ROI6_2(:,k))/4):3*round(numel(aux_VBSD_ROI6_2(:,k))/4),k));

                VBSD_ROI3_act_4(1,k) = mean(aux_VBSD_ROI3_2(3*round(numel(aux_VBSD_ROI3_2(:,k))/4):4*(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_4(1,k) = mean(aux_VBSD_ROI4_2(3*round(numel(aux_VBSD_ROI4_2(:,k))/4):4*(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_4(1,k) = mean(aux_VBSD_ROI5_2(3*round(numel(aux_VBSD_ROI5_2(:,k))/4):4*(numel(aux_VBSD_ROI5_2(:,k))/4),k));
                VBSD_ROI6_act_4(1,k) = mean(aux_VBSD_ROI6_2(3*round(numel(aux_VBSD_ROI6_2(:,k))/4):4*(numel(aux_VBSD_ROI6_2(:,k))/4),k));
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                BSD_ROI3_act_1(1,k) = mean(aux_BSD_ROI3(1:round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_1(1,k) = mean(aux_BSD_ROI4(1:round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_1(1,k) = mean(aux_BSD_ROI5(1:round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_1(1,k) = mean(aux_BSD_ROI6(1:round(numel(aux_BSD_ROI6(:,k))/4),k));


                BSD_ROI3_act_2(1,k) = mean(aux_BSD_ROI3(round(numel(aux_BSD_ROI3(:,k))/4):2*round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_2(1,k) = mean(aux_BSD_ROI4(round(numel(aux_BSD_ROI4(:,k))/4):2*round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_2(1,k) = mean(aux_BSD_ROI5(round(numel(aux_BSD_ROI5(:,k))/4):2*round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_2(1,k) = mean(aux_BSD_ROI6(round(numel(aux_BSD_ROI6(:,k))/4):2*round(numel(aux_BSD_ROI6(:,k))/4),k));


                BSD_ROI3_act_3(1,k) = mean(aux_BSD_ROI3(2*round(numel(aux_BSD_ROI3(:,k))/4):3*round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_3(1,k) = mean(aux_BSD_ROI4(2*round(numel(aux_BSD_ROI4(:,k))/4):3*round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_3(1,k) = mean(aux_BSD_ROI5(2*round(numel(aux_BSD_ROI5(:,k))/4):3*round(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_3(1,k) = mean(aux_BSD_ROI6(2*round(numel(aux_BSD_ROI6(:,k))/4):3*round(numel(aux_BSD_ROI6(:,k))/4),k));


                BSD_ROI3_act_4(1,k) = mean(aux_BSD_ROI3(3*round(numel(aux_BSD_ROI3(:,k))/4):4*(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_4(1,k) = mean(aux_BSD_ROI4(3*round(numel(aux_BSD_ROI4(:,k))/4):4*(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_4(1,k) = mean(aux_BSD_ROI5(3*round(numel(aux_BSD_ROI5(:,k))/4):4*(numel(aux_BSD_ROI5(:,k))/4),k));
                BSD_ROI6_act_4(1,k) = mean(aux_BSD_ROI6(3*round(numel(aux_BSD_ROI6(:,k))/4):4*(numel(aux_BSD_ROI6(:,k))/4),k));
            end

    % store in VBSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages  
            VBSD_avg_ROI3(ROI3_idx).ROI3 = VBSD_ROI3_avg_2(:)';  
            VBSD_avg_ROI4(ROI3_idx).ROI3 = VBSD_ROI4_avg_2(:)';  
            VBSD_avg_ROI5(ROI3_idx).ROI3 = VBSD_ROI5_avg_2(:)';   
            VBSD_avg_ROI6(ROI3_idx).ROI3 = VBSD_ROI6_avg_2(:)';   

            VBSD_act_1_ROI3(ROI3_idx).ROI3 = VBSD_ROI3_act_1(:)'; 
            VBSD_act_1_ROI4(ROI3_idx).ROI3 = VBSD_ROI4_act_1(:)';  
            VBSD_act_1_ROI5(ROI3_idx).ROI3 = VBSD_ROI5_act_1(:)';   
            VBSD_act_1_ROI6(ROI3_idx).ROI3 = VBSD_ROI6_act_1(:)';   

            VBSD_act_2_ROI3(ROI3_idx).ROI3 = VBSD_ROI3_act_2(:)'; 
            VBSD_act_2_ROI4(ROI3_idx).ROI3 = VBSD_ROI4_act_2(:)';  
            VBSD_act_2_ROI5(ROI3_idx).ROI3 = VBSD_ROI5_act_2(:)';   
            VBSD_act_2_ROI6(ROI3_idx).ROI3 = VBSD_ROI6_act_2(:)';   

            VBSD_act_3_ROI3(ROI3_idx).ROI3 = VBSD_ROI3_act_3(:)'; 
            VBSD_act_3_ROI4(ROI3_idx).ROI3 = VBSD_ROI4_act_3(:)';  
            VBSD_act_3_ROI5(ROI3_idx).ROI3 = VBSD_ROI5_act_3(:)';   
            VBSD_act_3_ROI6(ROI3_idx).ROI3 = VBSD_ROI6_act_3(:)';   

            VBSD_act_4_ROI3(ROI3_idx).ROI3 = VBSD_ROI3_act_4(:)'; 
            VBSD_act_4_ROI4(ROI3_idx).ROI3 = VBSD_ROI4_act_4(:)';  
            VBSD_act_4_ROI5(ROI3_idx).ROI3 = VBSD_ROI5_act_4(:)';   
            VBSD_act_4_ROI6(ROI3_idx).ROI3 = VBSD_ROI6_act_4(:)';   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BSD_act_1_ROI3(ROI3_idx).ROI3 = BSD_ROI3_act_1(:)'; 
            BSD_act_1_ROI4(ROI3_idx).ROI3 = BSD_ROI4_act_1(:)';  
            BSD_act_1_ROI5(ROI3_idx).ROI3 = BSD_ROI5_act_1(:)';   
            BSD_act_1_ROI6(ROI3_idx).ROI3 = BSD_ROI6_act_1(:)';   

            BSD_act_2_ROI3(ROI3_idx).ROI3 = BSD_ROI3_act_2(:)'; 
            BSD_act_2_ROI4(ROI3_idx).ROI3 = BSD_ROI4_act_2(:)';  
            BSD_act_2_ROI5(ROI3_idx).ROI3 = BSD_ROI5_act_2(:)';   
            BSD_act_2_ROI6(ROI3_idx).ROI3 = BSD_ROI6_act_2(:)';   

            BSD_act_3_ROI3(ROI3_idx).ROI3 = BSD_ROI3_act_3(:)'; 
            BSD_act_3_ROI4(ROI3_idx).ROI3 = BSD_ROI4_act_3(:)';  
            BSD_act_3_ROI5(ROI3_idx).ROI3 = BSD_ROI5_act_3(:)';   
            BSD_act_3_ROI6(ROI3_idx).ROI3 = BSD_ROI6_act_3(:)';   

            BSD_act_4_ROI3(ROI3_idx).ROI3 = BSD_ROI3_act_4(:)'; 
            BSD_act_4_ROI4(ROI3_idx).ROI3 = BSD_ROI4_act_4(:)';  
            BSD_act_4_ROI5(ROI3_idx).ROI3 = BSD_ROI5_act_4(:)';   
            BSD_act_4_ROI6(ROI3_idx).ROI3 = BSD_ROI6_act_4(:)';   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store results in matrices for further analysis 
            if exist('aux_BSD_ROI3_2','var') == 1
                ROI3_mtx_idx = ROI3_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI3_mtx(ROI3_mtx_idx,k) = aux_BSD_ROI3_2(1,k);
                    
                    VBSD_avg_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_avg_2(1,k);
                    VBSD_act_1_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_1(1,k);
                    VBSD_act_2_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_2(1,k);
                    VBSD_act_3_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_3(1,k);
                    VBSD_act_4_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_4(1,k);

                    BSD_act_1_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_1(1,k);
                    BSD_act_2_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_2(1,k);
                    BSD_act_3_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_3(1,k);
                    BSD_act_4_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI3(ROI3_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI3{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI3_mtx(ROI3_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI3_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI4_2','var') == 1
                ROI4_mtx_idx = ROI4_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI4_mtx(ROI4_mtx_idx,k) = aux_BSD_ROI4_2(1,k);
                    
                    VBSD_avg_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_avg_2(1,k);
                    VBSD_act_1_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_1(1,k);
                    VBSD_act_2_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_2(1,k);
                    VBSD_act_3_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_3(1,k);
                    VBSD_act_4_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_4(1,k);

                    BSD_act_1_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_1(1,k);
                    BSD_act_2_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_2(1,k);
                    BSD_act_3_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_3(1,k);
                    BSD_act_4_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI4(ROI4_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI4{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI4_mtx(ROI4_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI4_2(i-NUM_in+1,k); 
                end
            end
            end   
            if exist('aux_BSD_ROI5_2','var') == 1
                ROI5_mtx_idx = ROI5_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI5_mtx(ROI5_mtx_idx,k) = aux_BSD_ROI5_2(1,k);
                    
                    VBSD_avg_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_avg_2(1,k);
                    VBSD_act_1_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_1(1,k);
                    VBSD_act_2_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_2(1,k);
                    VBSD_act_3_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_3(1,k);
                    VBSD_act_4_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_4(1,k);

                    BSD_act_1_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_1(1,k);
                    BSD_act_2_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_2(1,k);
                    BSD_act_3_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_3(1,k);
                    BSD_act_4_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI5(ROI5_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI5{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI5_mtx(ROI5_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI5_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI6_2','var') == 1
                ROI6_mtx_idx = ROI6_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI6_mtx(ROI6_mtx_idx,k) = aux_BSD_ROI6_2(1,k);
                    
                    VBSD_avg_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_avg_2(1,k);
                    VBSD_act_1_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_1(1,k);
                    VBSD_act_2_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_2(1,k);
                    VBSD_act_3_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_3(1,k);
                    VBSD_act_4_ROI6_mtx(ROI6_mtx_idx ,k) = VBSD_ROI6_act_4(1,k);

                    BSD_act_1_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_1(1,k);
                    BSD_act_2_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_2(1,k);
                    BSD_act_3_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_3(1,k);
                    BSD_act_4_ROI6_mtx(ROI6_mtx_idx ,k) = BSD_ROI6_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI6(ROI6_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI6{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI6_mtx(ROI6_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI6_2(i-NUM_in+1,k); 
                end
            end
            end
     % save results
            cd(dir_save)
            save('BSDX_median.mat','BSDX_median','-v7.3');
            save(['BSD_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI3','-v7.3'); 
            save(['BSD_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4','-v7.3');  
            save(['BSD_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI5','-v7.3');   
            save(['BSD_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI6','-v7.3');   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % save results
            save(['VBSD_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI3_2','-v7.3'); 
            save(['VBSD_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI4_2','-v7.3'); 
            save(['VBSD_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI5_2','-v7.3');   
            save(['VBSD_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI6_2','-v7.3');   

    %     save(['VBSD_avg_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_avg_2','-v7.3');     
    %     save(['VBSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_avg_2','-v7.3');   % name structure of images (e.g. ...peak_0001)
    %     save(['VBSD_avg_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_avg_2','-v7.3');   
    % 
    %     save(['VBSD_act_1_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_1','-v7.3');
    %     save(['VBSD_act_1_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_1','-v7.3');   
    %     
    %     save(['VBSD_act_2_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_2','-v7.3');
    %     save(['VBSD_act_2_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_2','-v7.3');   
    % 
    %     save(['VBSD_act_3_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_3','-v7.3');
    %     save(['VBSD_act_3_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_3','-v7.3');   
    %    
    %     save(['VBSD_act_4_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_4','-v7.3');
    %     save(['VBSD_act_4_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI6_act_4','-v7.3');   

            save(['alpha_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI3','-v7.3'); 
            save(['alpha_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI4','-v7.3');  
            save(['alpha_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI5','-v7.3');   
            save(['alpha_ins_ROI6_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI6','-v7.3');   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(ROI{case_idx},'ROI2')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear bsd_sums aux_BSD_ROI2_2 aux_BSD_ROI3_2 aux_BSD_ROI4_2 aux_BSD_ROI5_2  
            clear centers_ROI4 radii_ROI4 bsd_sums_ROI4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = NUM_in:NUM
%% for methodology description, see details in oldROI5 part.
                for k = 2:round(nbin_radii_hg)+1
                    count_ROI2 = 0;
                    count_ROI3 = 0;
                    count_ROI4 = 0; % no bubbles initially in each bin/each image
                    count_ROI5 = 0;
                    for j = 1:numel(centers_xy.(im_name{i-NUM_in+1})(:,2))
                        if centers_xy.(im_name{i-NUM_in+1})(j,2)  <= 2.*image_length_mm/4 
    %new:ROI2
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) 
                                count_ROI2 = count_ROI2 + 1;
                                BSD_ROI2{i-NUM_in+1}(1,k-1) =  count_ROI2;   % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI2{i-NUM_in+1}(1,k-1) =  count_ROI2;  % add 0 if it doesn't find radius for the kth bin
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 2.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 3.*image_length_mm/4
    %new:ROI3
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k)
                                count_ROI3 = count_ROI3 + 1;
                                BSD_ROI3{i-NUM_in+1}(1,k-1) =  count_ROI3;  
                            else
                                BSD_ROI3{i-NUM_in+1}(1,k-1) =  count_ROI3;  
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 3.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 4.*image_length_mm/4
    %new:ROI4
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI4 = count_ROI4 + 1;
                                BSD_ROI4{i-NUM_in+1}(1,k-1) =  count_ROI4; 
                            else
                                BSD_ROI4{i-NUM_in+1}(1,k-1) =  count_ROI4;  
                            end

                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 4.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 5.*image_length_mm/4
    %new:ROI5
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) 
                                count_ROI5 = count_ROI5 + 1;
                                BSD_ROI5{i-NUM_in+1}(1,k-1) =  count_ROI5; 
                            else
                                BSD_ROI5{i-NUM_in+1}(1,k-1) =  count_ROI5; 
                            end
                        else
                        end
                    end
                end
            end
            for i=NUM_in:NUM
                    if isempty(BSD_ROI2{i-NUM_in+1}) 
                        BSD_ROI2{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI3{i-NUM_in+1}) 
                        BSD_ROI3{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI4{i-NUM_in+1}) 
                        BSD_ROI4{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI5{i-NUM_in+1}) 
                        BSD_ROI5{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
            end
            clear centers_ROI4 radii_ROI4
            for i=NUM_in:NUM
%% for methodology description, see details in oldROI5 part.
    % VBSD per unit length in 2 steps
    % 1. sort radii for each new ROI/each image according to the BSD_ROI4-7 findings
                bsd_sums = [1,sum(BSD_ROI2{1,i-NUM_in+1}),sum(BSD_ROI2{1,i-NUM_in+1}+BSD_ROI3{1,i-NUM_in+1}),...
                sum(BSD_ROI2{1,i-NUM_in+1}+BSD_ROI3{1,i-NUM_in+1}+BSD_ROI4{1,i-NUM_in+1}),...
                sum(BSD_ROI2{1,i-NUM_in+1}+BSD_ROI3{1,i-NUM_in+1}+BSD_ROI4{1,i-NUM_in+1}+BSD_ROI5{1,i-NUM_in+1})]; % indices to calculate Volume for each ROI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1 find total bubble numbers per BSD of new ROI
    %             bsd_toal_sums.(im_name{i-NUM_in+1}) = sum(BSD_ROI4{1,i-NUM_in+1});
    % 1.2 find radii and centers corresponding per BSD in new ROI
                radii_ROI4.(im_name{i-NUM_in+1})  = radii_mm_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4));
                centers_ROI4.(im_name{i-NUM_in+1})(:,1)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4),1);
                centers_ROI4.(im_name{i-NUM_in+1})(:,2)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4),2);
                centers_ROI4.(im_name{i-NUM_in+1})       = round(centers_ROI4.(im_name{i-NUM_in+1}),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort 
                vol_b_xy_ROI2.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2)));
                vol_b_xy_ROI3.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3)));
                vol_b_xy_ROI4.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4)));
                vol_b_xy_ROI5.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5)));
    % 2. calculate VBSD for each new ROI/each bin/each image 
                count_vol_ROI2{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI2
                count_vol_ROI3{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI3
                count_vol_ROI4{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI4
                count_vol_ROI5{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI5
                for k = 2:numel(BSDX_median)+1
    % calculate volume as sum of volumes of bubbles of each bin found in BSD_ROI
                    count_vol_ROI2{i-NUM_in+1}(1,k)  = count_vol_ROI2{i-NUM_in+1}(1,k-1) + BSD_ROI2{i-NUM_in+1}(1,k-1);  
                    aux_VBSD_ROI2{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI2.(im_name{i-NUM_in+1})...           % find number of bubbles per bin in specific ROI
                        (count_vol_ROI2{i-NUM_in+1}(1,k-1):count_vol_ROI2{i-NUM_in+1}(1,k)-1)); 

                    count_vol_ROI3{i-NUM_in+1}(1,k)  = count_vol_ROI3{i-NUM_in+1}(1,k-1) + BSD_ROI3{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI3{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI3.(im_name{i-NUM_in+1})...
                        (count_vol_ROI3{i-NUM_in+1}(1,k-1):count_vol_ROI3{i-NUM_in+1}(1,k)-1)); 

                    count_vol_ROI4{i-NUM_in+1}(1,k)  = count_vol_ROI4{i-NUM_in+1}(1,k-1) + BSD_ROI4{i-NUM_in+1}(1,k-1);  % find indices of bubbles that are in same bins in BSD 
                    aux_VBSD_ROI4{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI4.(im_name{i-NUM_in+1})...
                        (count_vol_ROI4{i-NUM_in+1}(1,k-1):count_vol_ROI4{i-NUM_in+1}(1,k)-1)); % calculate the sum of volumes of bubbles in same bin

                    count_vol_ROI5{i-NUM_in+1}(1,k)  = count_vol_ROI5{i-NUM_in+1}(1,k-1) + BSD_ROI5{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI5{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI5.(im_name{i-NUM_in+1})...
                        (count_vol_ROI5{i-NUM_in+1}(1,k-1):count_vol_ROI5{i-NUM_in+1}(1,k)-1)); 
                end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate instanteneous void fraction for ROI1-8
                aux_alpha_ins_ROI2{i-NUM_in+1} = sum(aux_VBSD_ROI2{i-NUM_in+1});
                aux_alpha_ins_ROI3{i-NUM_in+1} = sum(aux_VBSD_ROI3{i-NUM_in+1});
                aux_alpha_ins_ROI4{i-NUM_in+1} = sum(aux_VBSD_ROI4{i-NUM_in+1});
                aux_alpha_ins_ROI5{i-NUM_in+1} = sum(aux_VBSD_ROI5{i-NUM_in+1});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
    % save txt files    
    % indices of images to use in order to aling videos in time (see bubble_time_analysis.m for more info) 
            if strcmp(phase{case_idx},'peak')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'trough')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'posp2')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'minp2')
                   time_idx_in = 1;
                   time_idx_end = NUM;
            end
            % align instant data in time 
            for i = time_idx_in: time_idx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                radii_ins_ROI4.(im_name{i+1-time_idx_in}) = radii_ROI4.(im_name{i});
                centers_ins_ROI4.(im_name{i+1-time_idx_in}) = centers_ROI4.(im_name{i});
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cd('C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\GW_175\191020\server\entire_video\ROI1-8\text files');
            for i = 1: (time_idx_end-time_idx_in)
            fid=fopen(['centers_ins_ROI4.',(im_name{i}),'.old',ROI{case_idx},'.',run{case_idx},'.txt'],'wt');
            header1 = 'y(mm)';
            header2 = 'x(mm)';
            fprintf(fid, [ ' ' header1 ' ' header2 '\r\n']);
            fprintf(fid, '%4.1f %4.1f\r\n',centers_ROI4.(im_name{i})');
            fclose(fid);
            end
            for i = 1: (time_idx_end-time_idx_in)
            fid=fopen(['radii_ins_ROI4.',(im_name{i}),'.old',ROI{case_idx},'.',run{case_idx},'.txt'],'wt');
            header1 = 'r(mm)';
            fprintf(fid, [ ' ' header1 '\r\n']);
            fprintf(fid, '%4.2f\r\n',radii_ROI4.(im_name{i}));
            fclose(fid);
            end
            save(['radii_ins_ROI4_',phase{case_idx},'.old',ROI{case_idx},'.',run{case_idx},'.mat'],'-struct','radii_ROI4');
            save(['centers_ins_ROI4_',phase{case_idx},'.old',ROI{case_idx},'.',run{case_idx},'.mat'],'-struct','centers_ROI4');
            cd(dir_save)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %     save(['BSD_avg_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI2_2','-v7.3');   
    %     save(['BSD_avg_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI3_2','-v7.3');  
    %     save(['BSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4_2','-v7.3');   
    %     save(['BSD_avg_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI5_2','-v7.3');   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged BSD 
            for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_BSD_ROI2(i-NUM_in+1,k) = (BSD_ROI2{1,i-NUM_in+1}(k)); % place BSD values of each image/rep. in a 2D matrix 
                aux_BSD_ROI3(i-NUM_in+1,k) = (BSD_ROI3{1,i-NUM_in+1}(k));  
                aux_BSD_ROI4(i-NUM_in+1,k) = (BSD_ROI4{1,i-NUM_in+1}(k));  
                aux_BSD_ROI5(i-NUM_in+1,k) = (BSD_ROI5{1,i-NUM_in+1}(k));  
            end
            end
            for k = 1:numel(BSDX_median)
    % find mean per radius bin for each case and place them in cell matrix
                aux_BSD_ROI2_2(1,k) = mean(aux_BSD_ROI2(:,k));
                aux_BSD_ROI3_2(1,k) = mean(aux_BSD_ROI3(:,k));
                aux_BSD_ROI4_2(1,k) = mean(aux_BSD_ROI4(:,k));
                aux_BSD_ROI5_2(1,k) = mean(aux_BSD_ROI5(:,k));
            end
            ROI2_idx = ROI2_idx + 1; % new ROI2 index
    % store in BSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages 
            BSD_avg_ROI2(ROI2_idx).ROI2 = aux_BSD_ROI2_2(:)';   
            BSD_avg_ROI3(ROI2_idx).ROI2 = aux_BSD_ROI3_2(:)';  
            BSD_avg_ROI4(ROI2_idx).ROI2 = aux_BSD_ROI4_2(:)';  
            BSD_avg_ROI5(ROI2_idx).ROI2 = aux_BSD_ROI5_2(:)';  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged VBSD
            for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_VBSD_ROI2_2(i-NUM_in+1,k) = (aux_VBSD_ROI2{1,i-NUM_in+1}(k)); 
                aux_VBSD_ROI3_2(i-NUM_in+1,k) = (aux_VBSD_ROI3{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI4_2(i-NUM_in+1,k) = (aux_VBSD_ROI4{1,i-NUM_in+1}(k));  % place VBSD values of each image/rep. in a 2D matrix 
                aux_VBSD_ROI5_2(i-NUM_in+1,k) = (aux_VBSD_ROI5{1,i-NUM_in+1}(k));  
            end
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:numel(BSDX_median)
    % find mean per radius bin for each case and place them in cell matrix
                VBSD_ROI2_avg_2(1,k) = mean(aux_VBSD_ROI2_2(:,k));
                VBSD_ROI3_avg_2(1,k) = mean(aux_VBSD_ROI3_2(:,k));
                VBSD_ROI4_avg_2(1,k) = mean(aux_VBSD_ROI4_2(:,k));
                VBSD_ROI5_avg_2(1,k) = mean(aux_VBSD_ROI5_2(:,k));
% find time averages for segments of videos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% aux_VBSD_2 is equivalent to aux_BSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                VBSD_ROI2_act_1(1,k) = mean(aux_VBSD_ROI2_2(1:round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_1(1,k) = mean(aux_VBSD_ROI3_2(1:round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_1(1,k) = mean(aux_VBSD_ROI4_2(1:round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_1(1,k) = mean(aux_VBSD_ROI5_2(1:round(numel(aux_VBSD_ROI5_2(:,k))/4),k));

                VBSD_ROI2_act_2(1,k) = mean(aux_VBSD_ROI2_2(round(numel(aux_VBSD_ROI2_2(:,k))/4):2*round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_2(1,k) = mean(aux_VBSD_ROI3_2(round(numel(aux_VBSD_ROI3_2(:,k))/4):2*round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_2(1,k) = mean(aux_VBSD_ROI4_2(round(numel(aux_VBSD_ROI4_2(:,k))/4):2*round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_2(1,k) = mean(aux_VBSD_ROI5_2(round(numel(aux_VBSD_ROI5_2(:,k))/4):2*round(numel(aux_VBSD_ROI5_2(:,k))/4),k));

                VBSD_ROI2_act_3(1,k) = mean(aux_VBSD_ROI2_2(2*round(numel(aux_VBSD_ROI2_2(:,k))/4):3*round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_3(1,k) = mean(aux_VBSD_ROI3_2(2*round(numel(aux_VBSD_ROI3_2(:,k))/4):3*round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_3(1,k) = mean(aux_VBSD_ROI4_2(2*round(numel(aux_VBSD_ROI4_2(:,k))/4):3*round(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_3(1,k) = mean(aux_VBSD_ROI5_2(2*round(numel(aux_VBSD_ROI5_2(:,k))/4):3*round(numel(aux_VBSD_ROI5_2(:,k))/4),k));

                VBSD_ROI2_act_4(1,k) = mean(aux_VBSD_ROI2_2(3*round(numel(aux_VBSD_ROI2_2(:,k))/4):4*(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_4(1,k) = mean(aux_VBSD_ROI3_2(3*round(numel(aux_VBSD_ROI3_2(:,k))/4):4*(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_4(1,k) = mean(aux_VBSD_ROI4_2(3*round(numel(aux_VBSD_ROI4_2(:,k))/4):4*(numel(aux_VBSD_ROI4_2(:,k))/4),k));
                VBSD_ROI5_act_4(1,k) = mean(aux_VBSD_ROI5_2(3*round(numel(aux_VBSD_ROI5_2(:,k))/4):4*(numel(aux_VBSD_ROI5_2(:,k))/4),k));
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                BSD_ROI2_act_1(1,k) = mean(aux_BSD_ROI2(1:round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_1(1,k) = mean(aux_BSD_ROI3(1:round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_1(1,k) = mean(aux_BSD_ROI4(1:round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_1(1,k) = mean(aux_BSD_ROI5(1:round(numel(aux_BSD_ROI5(:,k))/4),k));


                BSD_ROI2_act_2(1,k) = mean(aux_BSD_ROI2(round(numel(aux_BSD_ROI2(:,k))/4):2*round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_2(1,k) = mean(aux_BSD_ROI3(round(numel(aux_BSD_ROI3(:,k))/4):2*round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_2(1,k) = mean(aux_BSD_ROI4(round(numel(aux_BSD_ROI4(:,k))/4):2*round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_2(1,k) = mean(aux_BSD_ROI5(round(numel(aux_BSD_ROI5(:,k))/4):2*round(numel(aux_BSD_ROI5(:,k))/4),k));


                BSD_ROI2_act_3(1,k) = mean(aux_BSD_ROI2(2*round(numel(aux_BSD_ROI2(:,k))/4):3*round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_3(1,k) = mean(aux_BSD_ROI3(2*round(numel(aux_BSD_ROI3(:,k))/4):3*round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_3(1,k) = mean(aux_BSD_ROI4(2*round(numel(aux_BSD_ROI4(:,k))/4):3*round(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_3(1,k) = mean(aux_BSD_ROI5(2*round(numel(aux_BSD_ROI5(:,k))/4):3*round(numel(aux_BSD_ROI5(:,k))/4),k));


                BSD_ROI2_act_4(1,k) = mean(aux_BSD_ROI2(3*round(numel(aux_BSD_ROI2(:,k))/4):4*(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_4(1,k) = mean(aux_BSD_ROI3(3*round(numel(aux_BSD_ROI3(:,k))/4):4*(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_4(1,k) = mean(aux_BSD_ROI4(3*round(numel(aux_BSD_ROI4(:,k))/4):4*(numel(aux_BSD_ROI4(:,k))/4),k));
                BSD_ROI5_act_4(1,k) = mean(aux_BSD_ROI5(3*round(numel(aux_BSD_ROI5(:,k))/4):4*(numel(aux_BSD_ROI5(:,k))/4),k));
            end

    % store in VBSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages  
            VBSD_avg_ROI2(ROI2_idx).ROI2 = VBSD_ROI2_avg_2(:)';   
            VBSD_avg_ROI3(ROI2_idx).ROI2 = VBSD_ROI3_avg_2(:)';  
            VBSD_avg_ROI4(ROI2_idx).ROI2 = VBSD_ROI4_avg_2(:)';  
            VBSD_avg_ROI5(ROI2_idx).ROI2 = VBSD_ROI5_avg_2(:)';   

            VBSD_act_1_ROI2(ROI2_idx).ROI2 = VBSD_ROI2_act_1(:)'; 
            VBSD_act_1_ROI3(ROI2_idx).ROI2 = VBSD_ROI3_act_1(:)';   
            VBSD_act_1_ROI4(ROI2_idx).ROI2 = VBSD_ROI4_act_1(:)';  
            VBSD_act_1_ROI5(ROI2_idx).ROI2 = VBSD_ROI5_act_1(:)';   

            VBSD_act_2_ROI2(ROI2_idx).ROI2 = VBSD_ROI2_act_2(:)';   
            VBSD_act_2_ROI3(ROI2_idx).ROI2 = VBSD_ROI3_act_2(:)'; 
            VBSD_act_2_ROI4(ROI2_idx).ROI2 = VBSD_ROI4_act_2(:)';  
            VBSD_act_2_ROI5(ROI2_idx).ROI2 = VBSD_ROI5_act_2(:)';   

            VBSD_act_3_ROI2(ROI2_idx).ROI2 = VBSD_ROI2_act_3(:)';   
            VBSD_act_3_ROI3(ROI2_idx).ROI2 = VBSD_ROI3_act_3(:)'; 
            VBSD_act_3_ROI4(ROI2_idx).ROI2 = VBSD_ROI4_act_3(:)';  
            VBSD_act_3_ROI5(ROI2_idx).ROI2 = VBSD_ROI5_act_3(:)';   

            VBSD_act_4_ROI2(ROI2_idx).ROI2 = VBSD_ROI2_act_4(:)';   
            VBSD_act_4_ROI3(ROI2_idx).ROI2 = VBSD_ROI3_act_4(:)'; 
            VBSD_act_4_ROI4(ROI2_idx).ROI2 = VBSD_ROI4_act_4(:)';  
            VBSD_act_4_ROI5(ROI2_idx).ROI2 = VBSD_ROI5_act_4(:)';   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BSD_act_1_ROI2(ROI2_idx).ROI2 = BSD_ROI2_act_1(:)'; 
            BSD_act_1_ROI3(ROI2_idx).ROI2 = BSD_ROI3_act_1(:)';   
            BSD_act_1_ROI4(ROI2_idx).ROI2 = BSD_ROI4_act_1(:)';  
            BSD_act_1_ROI5(ROI2_idx).ROI2 = BSD_ROI5_act_1(:)';   

            BSD_act_2_ROI2(ROI2_idx).ROI2 = BSD_ROI2_act_2(:)';   
            BSD_act_2_ROI3(ROI2_idx).ROI2 = BSD_ROI3_act_2(:)'; 
            BSD_act_2_ROI4(ROI2_idx).ROI2 = BSD_ROI4_act_2(:)';  
            BSD_act_2_ROI5(ROI2_idx).ROI2 = BSD_ROI5_act_2(:)';   

            BSD_act_3_ROI2(ROI2_idx).ROI2 = BSD_ROI2_act_3(:)';   
            BSD_act_3_ROI3(ROI2_idx).ROI2 = BSD_ROI3_act_3(:)'; 
            BSD_act_3_ROI4(ROI2_idx).ROI2 = BSD_ROI4_act_3(:)';  
            BSD_act_3_ROI5(ROI2_idx).ROI2 = BSD_ROI5_act_3(:)';   

            BSD_act_4_ROI2(ROI2_idx).ROI2 = BSD_ROI2_act_4(:)';   
            BSD_act_4_ROI3(ROI2_idx).ROI2 = BSD_ROI3_act_4(:)'; 
            BSD_act_4_ROI4(ROI2_idx).ROI2 = BSD_ROI4_act_4(:)';  
            BSD_act_4_ROI5(ROI2_idx).ROI2 = BSD_ROI5_act_4(:)';   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store results in matrices for further analysis 
            if exist('aux_BSD_ROI2_2','var') == 1
                ROI2_mtx_idx = ROI2_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI2_mtx(ROI2_mtx_idx,k) = aux_BSD_ROI2_2(1,k);
                    
                    VBSD_avg_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_avg_2(1,k);
                    VBSD_act_1_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_1(1,k);
                    VBSD_act_2_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_2(1,k);
                    VBSD_act_3_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_3(1,k);
                    VBSD_act_4_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_4(1,k);

                    BSD_act_1_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_1(1,k);
                    BSD_act_2_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_2(1,k);
                    BSD_act_3_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_3(1,k);
                    BSD_act_4_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI2(ROI2_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI2{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI2_mtx(ROI2_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI2_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI3_2','var') == 1
                ROI3_mtx_idx = ROI3_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI3_mtx(ROI3_mtx_idx,k) = aux_BSD_ROI3_2(1,k);
                    
                    VBSD_avg_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_avg_2(1,k);
                    VBSD_act_1_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_1(1,k);
                    VBSD_act_2_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_2(1,k);
                    VBSD_act_3_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_3(1,k);
                    VBSD_act_4_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_4(1,k);

                    BSD_act_1_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_1(1,k);
                    BSD_act_2_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_2(1,k);
                    BSD_act_3_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_3(1,k);
                    BSD_act_4_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI3(ROI3_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI3{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI3_mtx(ROI3_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI3_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI4_2','var') == 1
                ROI4_mtx_idx = ROI4_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI4_mtx(ROI4_mtx_idx,k) = aux_BSD_ROI4_2(1,k);
                    
                    VBSD_avg_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_avg_2(1,k);
                    VBSD_act_1_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_1(1,k);
                    VBSD_act_2_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_2(1,k);
                    VBSD_act_3_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_3(1,k);
                    VBSD_act_4_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_4(1,k);

                    BSD_act_1_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_1(1,k);
                    BSD_act_2_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_2(1,k);
                    BSD_act_3_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_3(1,k);
                    BSD_act_4_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI4(ROI4_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI4{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI4_mtx(ROI4_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI4_2(i-NUM_in+1,k); 
                end
            end
            end   
            if exist('aux_BSD_ROI5_2','var') == 1
                ROI5_mtx_idx = ROI5_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI5_mtx(ROI5_mtx_idx,k) = aux_BSD_ROI5_2(1,k);
                    
                    VBSD_avg_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_avg_2(1,k);
                    VBSD_act_1_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_1(1,k);
                    VBSD_act_2_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_2(1,k);
                    VBSD_act_3_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_3(1,k);
                    VBSD_act_4_ROI5_mtx(ROI5_mtx_idx,k) = VBSD_ROI5_act_4(1,k);

                    BSD_act_1_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_1(1,k);
                    BSD_act_2_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_2(1,k);
                    BSD_act_3_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_3(1,k);
                    BSD_act_4_ROI5_mtx(ROI5_mtx_idx,k) = BSD_ROI5_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI5(ROI5_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI5{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI5_mtx(ROI5_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI5_2(i-NUM_in+1,k); 
                end
            end
            end
    % save results
            cd(dir_save)
            save('BSDX_median.mat','BSDX_median','-v7.3');
            save(['BSD_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI2','-v7.3'); 
            save(['BSD_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI3','-v7.3'); 
            save(['BSD_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4','-v7.3');  
            save(['BSD_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI5','-v7.3');   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save results
            save(['VBSD_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI2_2','-v7.3');  
            save(['VBSD_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI3_2','-v7.3'); 
            save(['VBSD_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI4_2','-v7.3'); 
            save(['VBSD_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI5_2','-v7.3');   

    %     save(['VBSD_avg_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_avg_2','-v7.3');     
    %     save(['VBSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_avg_2','-v7.3');   
    % 
    %     save(['VBSD_act_1_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_1','-v7.3');
    %     save(['VBSD_act_1_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_1','-v7.3');   
    %     
    %     save(['VBSD_act_2_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_2','-v7.3');
    %     save(['VBSD_act_2_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_2','-v7.3');   
    % 
    %     save(['VBSD_act_3_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_3','-v7.3');
    %     save(['VBSD_act_3_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_3','-v7.3');   
    %    
    %     save(['VBSD_act_4_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_4','-v7.3'); 
    %     save(['VBSD_act_4_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_4','-v7.3');
    %     save(['VBSD_act_4_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI5_act_4','-v7.3');   

            save(['alpha_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI2','-v7.3');  
            save(['alpha_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI3','-v7.3'); 
            save(['alpha_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI4','-v7.3');  
            save(['alpha_ins_ROI5_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI5','-v7.3');   

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(ROI{case_idx},'ROI1')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear bsd_sums aux_BSD_ROI1_2 aux_BSD_ROI2_2 aux_BSD_ROI3_2 aux_BSD_ROI4_2   
            clear centers_ROI4 radii_ROI4 bsd_sums_ROI4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = NUM_in:NUM
%% for methodology description, see details in oldROI6 part.
    % single valued functions per unit length for ROI4/4
                for k = 2:round(nbin_radii_hg)+1
                    count_ROI1 = 0;
                    count_ROI2 = 0;
                    count_ROI3 = 0;
                    count_ROI4 = 0; % no bubbles initially in each bin/each image
                    for j = 1:numel(centers_xy.(im_name{i-NUM_in+1})(:,2))
                        if centers_xy.(im_name{i-NUM_in+1})(j,2)  <= 1.*image_length_mm/4 
    %new:ROI1
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k)                             
                                count_ROI1 = count_ROI1 + 1;
                                BSD_ROI1{i-NUM_in+1}(1,k-1) =  count_ROI1; % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI1{i-NUM_in+1}(1,k-1) =  count_ROI1; % add 0 if it doesn't find radius for the kth bin
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 1.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 2.*image_length_mm/4 
    %new:ROI2
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) 
                                count_ROI2 = count_ROI2 + 1;
                                BSD_ROI2{i-NUM_in+1}(1,k-1) =  count_ROI2;  
                            else
                                BSD_ROI2{i-NUM_in+1}(1,k-1) =  count_ROI2;  
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 2.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 3.*image_length_mm/4
    %new:ROI3
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k)
                                count_ROI3 = count_ROI3 + 1;
                                BSD_ROI3{i-NUM_in+1}(1,k-1) =  count_ROI3;  
                            else
                                BSD_ROI3{i-NUM_in+1}(1,k-1) =  count_ROI3;  
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 3.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 4.*image_length_mm/4
    %new:ROI4
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI4 = count_ROI4 + 1;
                                BSD_ROI4{i-NUM_in+1}(1,k-1) =  count_ROI4;  
                            else
                                BSD_ROI4{i-NUM_in+1}(1,k-1) =  count_ROI4;  
                            end
                        end
                    end
                end
            end
            for i=NUM_in:NUM
                    if isempty(BSD_ROI1{i-NUM_in+1}) 
                        BSD_ROI1{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI2{i-NUM_in+1}) 
                        BSD_ROI2{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI3{i-NUM_in+1}) 
                        BSD_ROI3{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI4{i-NUM_in+1}) 
                        BSD_ROI4{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
            end
            clear radii_ROI4 centers_ROI4
            for i=NUM_in:NUM
    % VBSD per unit length
%% for methodology description, see details in oldROI5 part.
    % 1. sort radii for each new ROI/each image according to the BSD_ROI1-4 findings
                bsd_sums = [1,sum(BSD_ROI1{1,i-NUM_in+1}),sum(BSD_ROI1{1,i-NUM_in+1}+BSD_ROI2{1,i-NUM_in+1}),...
                sum(BSD_ROI1{1,i-NUM_in+1}+BSD_ROI2{1,i-NUM_in+1}+BSD_ROI3{1,i-NUM_in+1}),...
                sum(BSD_ROI1{1,i-NUM_in+1}+BSD_ROI2{1,i-NUM_in+1}+BSD_ROI3{1,i-NUM_in+1}+BSD_ROI4{1,i-NUM_in+1})]; % indices to calculate Volume for each ROI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1 find total bubble numbers per BSD of new ROI
    %             bsd_toal_sums.(im_name{i-NUM_in+1}) = sum(BSD_ROI4{1,i-NUM_in+1});
    % 1.2 find radii and centers corresponding per BSD in new ROI
                radii_ROI4.(im_name{i-NUM_in+1})  = radii_mm_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5));
                centers_ROI4.(im_name{i-NUM_in+1})(:,1)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5),1);
                centers_ROI4.(im_name{i-NUM_in+1})(:,2)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5),2);
                centers_ROI4.(im_name{i-NUM_in+1})       = round(centers_ROI4.(im_name{i-NUM_in+1}),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort 
                vol_b_xy_ROI1.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2)));
                vol_b_xy_ROI2.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3)));
                vol_b_xy_ROI3.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4)));
                vol_b_xy_ROI4.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5)));
    % 2. calculate VBSD for each new ROI/each bin/each image 
                count_vol_ROI1{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI1            
                count_vol_ROI2{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI2
                count_vol_ROI3{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI3
                count_vol_ROI4{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI4
                for k = 2:numel(BSDX_median)+1
    % calculate volume as sum of volumes of bubbles of each bin found in BSD_ROI
                    count_vol_ROI1{i-NUM_in+1}(1,k)  = count_vol_ROI1{i-NUM_in+1}(1,k-1) + BSD_ROI1{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI1{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI1.(im_name{i-NUM_in+1})...           % find number of bubbles per bin in specific ROI
                        (count_vol_ROI1{i-NUM_in+1}(1,k-1):count_vol_ROI1{i-NUM_in+1}(1,k)-1));  % calculate the sum of volumes of bubbles in same bin
                    count_vol_ROI2{i-NUM_in+1}(1,k)  = count_vol_ROI2{i-NUM_in+1}(1,k-1) + BSD_ROI2{i-NUM_in+1}(1,k-1);  
                    aux_VBSD_ROI2{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI2.(im_name{i-NUM_in+1})...
                        (count_vol_ROI2{i-NUM_in+1}(1,k-1):count_vol_ROI2{i-NUM_in+1}(1,k)-1)); 
                    count_vol_ROI3{i-NUM_in+1}(1,k)  = count_vol_ROI3{i-NUM_in+1}(1,k-1) + BSD_ROI3{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI3{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI3.(im_name{i-NUM_in+1})...
                        (count_vol_ROI3{i-NUM_in+1}(1,k-1):count_vol_ROI3{i-NUM_in+1}(1,k)-1)); 
                    count_vol_ROI4{i-NUM_in+1}(1,k)  = count_vol_ROI4{i-NUM_in+1}(1,k-1) + BSD_ROI4{i-NUM_in+1}(1,k-1);  
                    aux_VBSD_ROI4{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI4.(im_name{i-NUM_in+1})...
                        (count_vol_ROI4{i-NUM_in+1}(1,k-1):count_vol_ROI4{i-NUM_in+1}(1,k)-1)); 
                end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate instanteneous void fraction for ROI1-8
                aux_alpha_ins_ROI1{i-NUM_in+1} = sum(aux_VBSD_ROI1{i-NUM_in+1});
                aux_alpha_ins_ROI2{i-NUM_in+1} = sum(aux_VBSD_ROI2{i-NUM_in+1});
                aux_alpha_ins_ROI3{i-NUM_in+1} = sum(aux_VBSD_ROI3{i-NUM_in+1});
                aux_alpha_ins_ROI4{i-NUM_in+1} = sum(aux_VBSD_ROI4{i-NUM_in+1});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end

    % save txt files    
    % indices of images to use in order to aling videos in time (see bubble_time_analysis.m for more info) 
            if strcmp(phase{case_idx},'peak')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'trough')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'posp2')
                   time_idx_in = 1;
                   time_idx_end = NUM;

            elseif strcmp(phase{case_idx},'minp2')
                   time_idx_in = 1;
                   time_idx_end = NUM;
            end
            % align instant data in time 
            for i = time_idx_in: time_idx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                radii_ins_ROI4.(im_name{i+1-time_idx_in}) = radii_ROI4.(im_name{i});
                centers_ins_ROI4.(im_name{i+1-time_idx_in}) = centers_ROI4.(im_name{i});
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cd(['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\',wave,'\',date,'\server\entire_video\ROI1-8\text files']);
            for i = 1:(time_idx_end - time_idx_in)
            fid=fopen(['centers_ins_ROI4.',(im_name{i}),'.old',ROI{case_idx},'.',run{case_idx},'.txt'],'wt');
            header1 = 'y(mm)';
            header2 = 'x(mm)';
            fprintf(fid, [ ' ' header1 ' ' header2 '\r\n']);
            fprintf(fid, '%4.1f %4.1f\r\n',centers_ROI4.(im_name{i})');
            fclose(fid);
            end
            for i = 1:(time_idx_end - time_idx_in)
            fid=fopen(['radii_ins_ROI4.',(im_name{i}),'.old',ROI{case_idx},'.',run{case_idx},'.txt'],'wt');
            header1 = 'r(mm)';
            fprintf(fid, [ ' ' header1 '\r\n']);
            fprintf(fid, '%4.2f\r\n',radii_ROI4.(im_name{i}));
            fclose(fid);
            end
            save(['radii_ins_ROI4_',phase{case_idx},'.old',ROI{case_idx},'.',run{case_idx},'.mat'],'-struct','radii_ROI4');
            save(['centers_ins_ROI4_',phase{case_idx},'.old',ROI{case_idx},'.',run{case_idx},'.mat'],'-struct','centers_ROI4');
            cd(dir_save)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %     save(['BSD_avg_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI1_2','-v7.3');   
    %     save(['BSD_avg_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI2_2','-v7.3');   
    %     save(['BSD_avg_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI3_2','-v7.3');  
    %     save(['BSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4_2','-v7.3');  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged BSD 
            for i=NUM_in:NUM
                for k = 1:numel(BSDX_median)
                    aux_BSD_ROI1(i-NUM_in+1,k) = (BSD_ROI1{1,i-NUM_in+1}(k));  % place BSD values of each image/rep. in a 2D matrix 
                    aux_BSD_ROI2(i-NUM_in+1,k) = (BSD_ROI2{1,i-NUM_in+1}(k)); 
                    aux_BSD_ROI3(i-NUM_in+1,k) = (BSD_ROI3{1,i-NUM_in+1}(k));  
                    aux_BSD_ROI4(i-NUM_in+1,k) = (BSD_ROI4{1,i-NUM_in+1}(k));  
                end
            end
            for k = 1:numel(BSDX_median)
    % find mean per radius bin for each case and place them in cell matrix
                aux_BSD_ROI1_2(1,k) = mean(aux_BSD_ROI1(:,k));
                aux_BSD_ROI2_2(1,k) = mean(aux_BSD_ROI2(:,k));
                aux_BSD_ROI3_2(1,k) = mean(aux_BSD_ROI3(:,k));
                aux_BSD_ROI4_2(1,k) = mean(aux_BSD_ROI4(:,k));
            end
            ROI1_idx = ROI1_idx + 1; % new ROI1 index
    % store in BSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages 
            BSD_avg_ROI1(ROI1_idx).ROI1 = aux_BSD_ROI1_2(:)';   
            BSD_avg_ROI2(ROI1_idx).ROI1 = aux_BSD_ROI2_2(:)';   
            BSD_avg_ROI3(ROI1_idx).ROI1 = aux_BSD_ROI3_2(:)';  
            BSD_avg_ROI4(ROI1_idx).ROI1 = aux_BSD_ROI4_2(:)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged VBSD
            for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_VBSD_ROI1_2(i-NUM_in+1,k) = (aux_VBSD_ROI1{1,i-NUM_in+1}(k));  % place VBSD values of each image/rep. in a 2D matrix  
                aux_VBSD_ROI2_2(i-NUM_in+1,k) = (aux_VBSD_ROI2{1,i-NUM_in+1}(k)); 
                aux_VBSD_ROI3_2(i-NUM_in+1,k) = (aux_VBSD_ROI3{1,i-NUM_in+1}(k));  
                aux_VBSD_ROI4_2(i-NUM_in+1,k) = (aux_VBSD_ROI4{1,i-NUM_in+1}(k)); 
            end
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:numel(BSDX_median)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% aux_VBSD_2 is equivalent to aux_BSD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find mean per radius bin for each case and place them in cell matrix
                VBSD_ROI1_avg_2(1,k) = mean(aux_VBSD_ROI1_2(:,k));
                VBSD_ROI2_avg_2(1,k) = mean(aux_VBSD_ROI2_2(:,k));
                VBSD_ROI3_avg_2(1,k) = mean(aux_VBSD_ROI3_2(:,k));
                VBSD_ROI4_avg_2(1,k) = mean(aux_VBSD_ROI4_2(:,k));
    % find time averages for segments of videos
                VBSD_ROI1_act_1(1,k) = mean(aux_VBSD_ROI1_2(1:round(numel(aux_VBSD_ROI1_2(:,k))/4),k));
                VBSD_ROI2_act_1(1,k) = mean(aux_VBSD_ROI2_2(1:round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_1(1,k) = mean(aux_VBSD_ROI3_2(1:round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_1(1,k) = mean(aux_VBSD_ROI4_2(1:round(numel(aux_VBSD_ROI4_2(:,k))/4),k));

                VBSD_ROI1_act_2(1,k) = mean(aux_VBSD_ROI1_2(round(numel(aux_VBSD_ROI1_2(:,k))/4):2*round(numel(aux_VBSD_ROI1_2(:,k))/4),k));
                VBSD_ROI2_act_2(1,k) = mean(aux_VBSD_ROI2_2(round(numel(aux_VBSD_ROI2_2(:,k))/4):2*round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_2(1,k) = mean(aux_VBSD_ROI3_2(round(numel(aux_VBSD_ROI3_2(:,k))/4):2*round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_2(1,k) = mean(aux_VBSD_ROI4_2(round(numel(aux_VBSD_ROI4_2(:,k))/4):2*round(numel(aux_VBSD_ROI4_2(:,k))/4),k));

                VBSD_ROI1_act_3(1,k) = mean(aux_VBSD_ROI1_2(2*round(numel(aux_VBSD_ROI1_2(:,k))/4):3*round(numel(aux_VBSD_ROI1_2(:,k))/4),k));
                VBSD_ROI2_act_3(1,k) = mean(aux_VBSD_ROI2_2(2*round(numel(aux_VBSD_ROI2_2(:,k))/4):3*round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_3(1,k) = mean(aux_VBSD_ROI3_2(2*round(numel(aux_VBSD_ROI3_2(:,k))/4):3*round(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_3(1,k) = mean(aux_VBSD_ROI4_2(2*round(numel(aux_VBSD_ROI4_2(:,k))/4):3*round(numel(aux_VBSD_ROI4_2(:,k))/4),k));

                VBSD_ROI1_act_4(1,k) = mean(aux_VBSD_ROI1_2(3*round(numel(aux_VBSD_ROI1_2(:,k))/4):4*(numel(aux_VBSD_ROI1_2(:,k))/4),k));
                VBSD_ROI2_act_4(1,k) = mean(aux_VBSD_ROI2_2(3*round(numel(aux_VBSD_ROI2_2(:,k))/4):4*(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_4(1,k) = mean(aux_VBSD_ROI3_2(3*round(numel(aux_VBSD_ROI3_2(:,k))/4):4*(numel(aux_VBSD_ROI3_2(:,k))/4),k));
                VBSD_ROI4_act_4(1,k) = mean(aux_VBSD_ROI4_2(3*round(numel(aux_VBSD_ROI4_2(:,k))/4):4*(numel(aux_VBSD_ROI4_2(:,k))/4),k));
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------

                BSD_ROI1_act_1(1,k) = mean(aux_BSD_ROI1(1:round(numel(aux_BSD_ROI1(:,k))/4),k));
                BSD_ROI2_act_1(1,k) = mean(aux_BSD_ROI2(1:round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_1(1,k) = mean(aux_BSD_ROI3(1:round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_1(1,k) = mean(aux_BSD_ROI4(1:round(numel(aux_BSD_ROI4(:,k))/4),k));


                BSD_ROI1_act_2(1,k) = mean(aux_BSD_ROI1(round(numel(aux_BSD_ROI1(:,k))/4):2*round(numel(aux_BSD_ROI1(:,k))/4),k));
                BSD_ROI2_act_2(1,k) = mean(aux_BSD_ROI2(round(numel(aux_BSD_ROI2(:,k))/4):2*round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_2(1,k) = mean(aux_BSD_ROI3(round(numel(aux_BSD_ROI3(:,k))/4):2*round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_2(1,k) = mean(aux_BSD_ROI4(round(numel(aux_BSD_ROI4(:,k))/4):2*round(numel(aux_BSD_ROI4(:,k))/4),k));


                BSD_ROI1_act_3(1,k) = mean(aux_BSD_ROI1(2*round(numel(aux_BSD_ROI1(:,k))/4):3*round(numel(aux_BSD_ROI1(:,k))/4),k));
                BSD_ROI2_act_3(1,k) = mean(aux_BSD_ROI2(2*round(numel(aux_BSD_ROI2(:,k))/4):3*round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_3(1,k) = mean(aux_BSD_ROI3(2*round(numel(aux_BSD_ROI3(:,k))/4):3*round(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_3(1,k) = mean(aux_BSD_ROI4(2*round(numel(aux_BSD_ROI4(:,k))/4):3*round(numel(aux_BSD_ROI4(:,k))/4),k));


                BSD_ROI1_act_4(1,k) = mean(aux_BSD_ROI1(3*round(numel(aux_BSD_ROI1(:,k))/4):4*(numel(aux_BSD_ROI1(:,k))/4),k));
                BSD_ROI2_act_4(1,k) = mean(aux_BSD_ROI2(3*round(numel(aux_BSD_ROI2(:,k))/4):4*(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_4(1,k) = mean(aux_BSD_ROI3(3*round(numel(aux_BSD_ROI3(:,k))/4):4*(numel(aux_BSD_ROI3(:,k))/4),k));
                BSD_ROI4_act_4(1,k) = mean(aux_BSD_ROI4(3*round(numel(aux_BSD_ROI4(:,k))/4):4*(numel(aux_BSD_ROI4(:,k))/4),k));

            end

    % store in VBSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages  
            VBSD_avg_ROI1(ROI1_idx).ROI1 = VBSD_ROI1_avg_2(:)';   
            VBSD_avg_ROI2(ROI1_idx).ROI1 = VBSD_ROI2_avg_2(:)';   
            VBSD_avg_ROI3(ROI1_idx).ROI1 = VBSD_ROI3_avg_2(:)';  
            VBSD_avg_ROI4(ROI1_idx).ROI1 = VBSD_ROI4_avg_2(:)';  

            VBSD_act_1_ROI1(ROI1_idx).ROI1 = VBSD_ROI1_act_1(:)'; 
            VBSD_act_1_ROI2(ROI1_idx).ROI1 = VBSD_ROI2_act_1(:)';  
            VBSD_act_1_ROI3(ROI1_idx).ROI1 = VBSD_ROI3_act_1(:)';   
            VBSD_act_1_ROI4(ROI1_idx).ROI1 = VBSD_ROI4_act_1(:)';   

            VBSD_act_2_ROI1(ROI1_idx).ROI1 = VBSD_ROI1_act_2(:)';   
            VBSD_act_2_ROI2(ROI1_idx).ROI1 = VBSD_ROI2_act_2(:)';   
            VBSD_act_2_ROI3(ROI1_idx).ROI1 = VBSD_ROI3_act_2(:)'; 
            VBSD_act_2_ROI4(ROI1_idx).ROI1 = VBSD_ROI4_act_2(:)';  

            VBSD_act_3_ROI1(ROI1_idx).ROI1 = BSD_ROI1_act_3(:)';   
            VBSD_act_3_ROI2(ROI1_idx).ROI1 = BSD_ROI2_act_3(:)';   
            VBSD_act_3_ROI3(ROI1_idx).ROI1 = BSD_ROI3_act_3(:)'; 
            VBSD_act_3_ROI4(ROI1_idx).ROI1 = BSD_ROI4_act_3(:)';  

            VBSD_act_4_ROI1(ROI1_idx).ROI1 = BSD_ROI1_act_4(:)';   
            VBSD_act_4_ROI2(ROI1_idx).ROI1 = BSD_ROI2_act_4(:)';   
            VBSD_act_4_ROI3(ROI1_idx).ROI1 = BSD_ROI3_act_4(:)'; 
            VBSD_act_4_ROI4(ROI1_idx).ROI1 = BSD_ROI4_act_4(:)';  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            BSD_act_1_ROI1(ROI1_idx).ROI1 = BSD_ROI1_act_1(:)'; 
            BSD_act_1_ROI2(ROI1_idx).ROI1 = BSD_ROI2_act_1(:)';  
            BSD_act_1_ROI3(ROI1_idx).ROI1 = BSD_ROI3_act_1(:)';   
            BSD_act_1_ROI4(ROI1_idx).ROI1 = BSD_ROI4_act_1(:)';   

            BSD_act_2_ROI1(ROI1_idx).ROI1 = BSD_ROI1_act_2(:)';   
            BSD_act_2_ROI2(ROI1_idx).ROI1 = BSD_ROI2_act_2(:)';   
            BSD_act_2_ROI3(ROI1_idx).ROI1 = BSD_ROI3_act_2(:)'; 
            BSD_act_2_ROI4(ROI1_idx).ROI1 = BSD_ROI4_act_2(:)';  

            BSD_act_3_ROI1(ROI1_idx).ROI1 = BSD_ROI1_act_3(:)';   
            BSD_act_3_ROI2(ROI1_idx).ROI1 = BSD_ROI2_act_3(:)';   
            BSD_act_3_ROI3(ROI1_idx).ROI1 = BSD_ROI3_act_3(:)'; 
            BSD_act_3_ROI4(ROI1_idx).ROI1 = BSD_ROI4_act_3(:)';  

            BSD_act_4_ROI1(ROI1_idx).ROI1 = BSD_ROI1_act_4(:)';   
            BSD_act_4_ROI2(ROI1_idx).ROI1 = BSD_ROI2_act_4(:)';   
            BSD_act_4_ROI3(ROI1_idx).ROI1 = BSD_ROI3_act_4(:)'; 
            BSD_act_4_ROI4(ROI1_idx).ROI1 = BSD_ROI4_act_4(:)'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store results in matrices for further analysis 
            if exist('aux_BSD_ROI1_2','var') == 1
                ROI1_mtx_idx = ROI1_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI1_mtx(ROI1_mtx_idx,k) = aux_BSD_ROI1_2(1,k);
                    
                    VBSD_avg_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_avg_2(1,k);
                    
                    VBSD_act_1_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_act_1(1,k);
                    VBSD_act_2_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_act_2(1,k);
                    VBSD_act_3_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_act_3(1,k);
                    VBSD_act_4_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_act_4(1,k);

                    BSD_act_1_ROI1_mtx(ROI1_mtx_idx,k) = BSD_ROI1_act_1(1,k);
                    BSD_act_2_ROI1_mtx(ROI1_mtx_idx,k) = BSD_ROI1_act_2(1,k);
                    BSD_act_3_ROI1_mtx(ROI1_mtx_idx,k) = BSD_ROI1_act_3(1,k);
                    BSD_act_4_ROI1_mtx(ROI1_mtx_idx,k) = BSD_ROI1_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI1(ROI1_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI1{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI1_mtx(ROI1_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI1_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI2_2','var') == 1
                ROI2_mtx_idx = ROI2_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI2_mtx(ROI2_mtx_idx,k) = aux_BSD_ROI2_2(1,k);
                    VBSD_avg_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_avg_2(1,k);
                    VBSD_act_1_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_1(1,k);
                    VBSD_act_2_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_2(1,k);
                    VBSD_act_3_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_3(1,k);
                    VBSD_act_4_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_4(1,k);

                    BSD_act_1_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_1(1,k);
                    BSD_act_2_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_2(1,k);
                    BSD_act_3_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_3(1,k);
                    BSD_act_4_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI2(ROI2_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI2{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI2_mtx(ROI2_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI2_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI3_2','var') == 1
                ROI3_mtx_idx = ROI3_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI3_mtx(ROI3_mtx_idx,k) = aux_BSD_ROI3_2(1,k);
                    VBSD_avg_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_avg_2(1,k);
                    VBSD_act_1_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_1(1,k);
                    VBSD_act_2_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_2(1,k);
                    VBSD_act_3_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_3(1,k);
                    VBSD_act_4_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_4(1,k);

                    BSD_act_1_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_1(1,k);
                    BSD_act_2_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_2(1,k);
                    BSD_act_3_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_3(1,k);
                    BSD_act_4_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI3(ROI3_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI3{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI3_mtx(ROI3_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI3_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI4_2','var') == 1
                ROI4_mtx_idx = ROI4_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI4_mtx(ROI4_mtx_idx,k) = aux_BSD_ROI4_2(1,k);
                    VBSD_avg_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_avg_2(1,k);
                    VBSD_act_1_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_1(1,k);
                    VBSD_act_2_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_2(1,k);
                    VBSD_act_3_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_3(1,k);
                    VBSD_act_4_ROI4_mtx(ROI4_mtx_idx,k) = VBSD_ROI4_act_4(1,k);

                    BSD_act_1_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_1(1,k);
                    BSD_act_2_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_2(1,k);
                    BSD_act_3_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_3(1,k);
                    BSD_act_4_ROI4_mtx(ROI4_mtx_idx,k) = BSD_ROI4_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI4(ROI4_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI4{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI4_mtx(ROI4_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI4_2(i-NUM_in+1,k); 
                end
            end
            end
      % save results
            cd(dir_save)
            save('BSDX_median.mat','BSDX_median','-v7.3');
            save(['BSD_ins_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI1','-v7.3');   
            save(['BSD_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI2','-v7.3'); 
            save(['BSD_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI3','-v7.3'); 
            save(['BSD_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4','-v7.3');  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save results
            save(['VBSD_ins_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI1_2','-v7.3');   
            save(['VBSD_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI2_2','-v7.3');  
            save(['VBSD_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI3_2','-v7.3'); 
            save(['VBSD_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI4_2','-v7.3'); 

    %     save(['VBSD_avg_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_avg_2','-v7.3');     
    %     save(['VBSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_avg_2','-v7.3');   
    % 
    %     save(['VBSD_act_1_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_1','-v7.3');
    %     save(['VBSD_act_1_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_1','-v7.3');   
    % 
    %     save(['VBSD_act_2_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_2','-v7.3');
    %     save(['VBSD_act_2_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_2','-v7.3');   
    % 
    %     save(['VBSD_act_3_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_3','-v7.3');
    %     save(['VBSD_act_3_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_3','-v7.3');   
    %    
    %     save(['VBSD_act_4_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_4','-v7.3'); 
    %     save(['VBSD_act_4_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_4','-v7.3');
    %     save(['VBSD_act_4_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_4','-v7.3');   

            save(['alpha_ins_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI1','-v7.3');  
            save(['alpha_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI2','-v7.3');  
            save(['alpha_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI3','-v7.3'); 
            save(['alpha_ins_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI4','-v7.3');  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(ROI{case_idx},'ROI0')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear bsd_sums 
            clear aux_BSD_ROI1_0 bsd_sums aux_BSD_ROI1_2 aux_BSD_ROI2_2 aux_BSD_ROI3_2   
            clear centers_ROI4 radii_ROI4 bsd_sums_ROI4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = NUM_in:NUM
%% for methodology description, see details in oldROI6 part.
    % single valued functions per unit length for ROI4/4
                for k = 2:round(nbin_radii_hg)+1
                    count_ROI0 = 0;
                    count_ROI1 = 0;
                    count_ROI2 = 0;
                    count_ROI3 = 0; % no bubbles initially in each bin/each image
                    for j = 1:numel(centers_xy.(im_name{i-NUM_in+1})(:,2))
                        if centers_xy.(im_name{i-NUM_in+1})(j,2)  <=0
    %new:ROI0
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) % statement to check if bubble falls within the Kth bin size
                                count_ROI0 = count_ROI0 + 1;
                                BSD_ROI0{i-NUM_in+1}(1,k-1) =  count_ROI0;  
                            else
                                BSD_ROI0{i-NUM_in+1}(1,k-1) =  count_ROI0;  
                            end
                        
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  >0 && centers_xy.(im_name{i-NUM_in+1})(j,2)<= 1.*image_length_mm/4 
    %new:ROI1
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k)                             
                                count_ROI1 = count_ROI1 + 1;
                                BSD_ROI1{i-NUM_in+1}(1,k-1) =  count_ROI1; % add 1 for each detected radius in the kth bin
                            else
                                BSD_ROI1{i-NUM_in+1}(1,k-1) =  count_ROI1; % add 0 if it doesn't find radius for the kth bin
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 1.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 2.*image_length_mm/4 
    %new:ROI2
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k) 
                                count_ROI2 = count_ROI2 + 1;
                                BSD_ROI2{i-NUM_in+1}(1,k-1) =  count_ROI2;  
                            else
                                BSD_ROI2{i-NUM_in+1}(1,k-1) =  count_ROI2;  
                            end
                        elseif centers_xy.(im_name{i-NUM_in+1})(j,2)  > 2.*image_length_mm/4 && centers_xy.(im_name{i-NUM_in+1})(j,2) <= 3.*image_length_mm/4
    %new:ROI3
                            if radii_mm_xy.(im_name{i-NUM_in+1})(j)>=rad_bin(k-1) && radii_mm_xy.(im_name{i-NUM_in+1})(j)<rad_bin(k)
                                count_ROI3 = count_ROI3 + 1;
                                BSD_ROI3{i-NUM_in+1}(1,k-1) =  count_ROI3;  
                            else
                                BSD_ROI3{i-NUM_in+1}(1,k-1) =  count_ROI3;  
                            end
                        end
                    end
                end
            end
            for i=NUM_in:NUM
                    if isempty(BSD_ROI0{i-NUM_in+1}) 
                        BSD_ROI0{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI1{i-NUM_in+1}) 
                        BSD_ROI1{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI2{i-NUM_in+1}) 
                        BSD_ROI2{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
                    if isempty(BSD_ROI3{i-NUM_in+1}) 
                        BSD_ROI3{i-NUM_in+1}= zeros(1,numel(BSDX_median));
                    end
            end
            clear radii_ROI3 centers_ROI3
            for i=NUM_in:NUM
    % VBSD per unit length
%% for methodology description, see details in oldROI5 part.
    % 1. sort radii for each new ROI/each image according to the BSD_ROI1-4 findings
                bsd_sums = [1,sum(BSD_ROI0{1,i-NUM_in+1}),sum(BSD_ROI0{1,i-NUM_in+1}+BSD_ROI1{1,i-NUM_in+1}),...
                sum(BSD_ROI0{1,i-NUM_in+1}+BSD_ROI1{1,i-NUM_in+1}+BSD_ROI2{1,i-NUM_in+1}),...
                sum(BSD_ROI0{1,i-NUM_in+1}+BSD_ROI1{1,i-NUM_in+1}+BSD_ROI2{1,i-NUM_in+1}+BSD_ROI3{1,i-NUM_in+1})]; % indices to calculate Volume for each ROI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1 find total bubble numbers per BSD of new ROI
    %             bsd_toal_sums.(im_name{i-NUM_in+1}) = sum(BSD_ROI4{1,i-NUM_in+1});
    % 1.2 find radii and centers corresponding per BSD in new ROI
    %% 
                radii_ROI4.(im_name{i-NUM_in+1})  = radii_mm_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5));
                centers_ROI4.(im_name{i-NUM_in+1})(:,1)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5),1);
                centers_ROI4.(im_name{i-NUM_in+1})(:,2)  = centers_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5),2);
                centers_ROI4.(im_name{i-NUM_in+1})       = round(centers_ROI4.(im_name{i-NUM_in+1}),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort 
                vol_b_xy_ROI0.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(1):bsd_sums(2)));
                vol_b_xy_ROI1.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(2)+1:bsd_sums(3)));
                vol_b_xy_ROI2.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(3)+1:bsd_sums(4)));
                vol_b_xy_ROI3.(im_name{i-NUM_in+1}) = sort(vol_b_xy.(im_name{i-NUM_in+1})(bsd_sums(4)+1:bsd_sums(5)));
    % 2. calculate VBSD for each new ROI/each bin/each image 
                count_vol_ROI0{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI4
                count_vol_ROI1{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI1            
                count_vol_ROI2{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI2
                count_vol_ROI3{i-NUM_in+1}(1,1) = 1; % index of first bubble of first bin ROI3
                for k = 2:numel(BSDX_median)+1
    % calculate volume as sum of volumes of bubbles of each bin found in BSD_ROI
                    count_vol_ROI0{i-NUM_in+1}(1,k)  = count_vol_ROI0{i-NUM_in+1}(1,k-1) + BSD_ROI0{i-NUM_in+1}(1,k-1);  
                    aux_VBSD_ROI0{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI0.(im_name{i-NUM_in+1})...
                        (count_vol_ROI0{i-NUM_in+1}(1,k-1):count_vol_ROI0{i-NUM_in+1}(1,k)-1));     
    
                    count_vol_ROI1{i-NUM_in+1}(1,k)  = count_vol_ROI1{i-NUM_in+1}(1,k-1) + BSD_ROI1{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI1{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI1.(im_name{i-NUM_in+1})...           % find number of bubbles per bin in specific ROI
                        (count_vol_ROI1{i-NUM_in+1}(1,k-1):count_vol_ROI1{i-NUM_in+1}(1,k)-1));  % calculate the sum of volumes of bubbles in same bin
                    count_vol_ROI2{i-NUM_in+1}(1,k)  = count_vol_ROI2{i-NUM_in+1}(1,k-1) + BSD_ROI2{i-NUM_in+1}(1,k-1);  
                    aux_VBSD_ROI2{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI2.(im_name{i-NUM_in+1})...
                        (count_vol_ROI2{i-NUM_in+1}(1,k-1):count_vol_ROI2{i-NUM_in+1}(1,k)-1)); 
                    count_vol_ROI3{i-NUM_in+1}(1,k)  = count_vol_ROI3{i-NUM_in+1}(1,k-1) + BSD_ROI3{i-NUM_in+1}(1,k-1); 
                    aux_VBSD_ROI3{i-NUM_in+1}(1,k-1) = sum(vol_b_xy_ROI3.(im_name{i-NUM_in+1})...
                        (count_vol_ROI3{i-NUM_in+1}(1,k-1):count_vol_ROI3{i-NUM_in+1}(1,k)-1)); 

                end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate instanteneous void fraction for ROI1-8
                aux_alpha_ins_ROI0{i-NUM_in+1} = sum(aux_VBSD_ROI0{i-NUM_in+1});
                aux_alpha_ins_ROI1{i-NUM_in+1} = sum(aux_VBSD_ROI1{i-NUM_in+1});
                aux_alpha_ins_ROI2{i-NUM_in+1} = sum(aux_VBSD_ROI2{i-NUM_in+1});
                aux_alpha_ins_ROI3{i-NUM_in+1} = sum(aux_VBSD_ROI3{i-NUM_in+1});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %     save(['BSD_avg_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI1_2','-v7.3');   
    %     save(['BSD_avg_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI2_2','-v7.3');   
    %     save(['BSD_avg_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI3_2','-v7.3');  
    %     save(['BSD_avg_ROI4_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI4_2','-v7.3');  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged BSD 
            for i=NUM_in:NUM
                for k = 1:numel(BSDX_median)
                    aux_BSD_ROI0(i-NUM_in+1,k) = (BSD_ROI0{1,i-NUM_in+1}(k));  
                    aux_BSD_ROI1(i-NUM_in+1,k) = (BSD_ROI1{1,i-NUM_in+1}(k));  % place BSD values of each image/rep. in a 2D matrix 
                    aux_BSD_ROI2(i-NUM_in+1,k) = (BSD_ROI2{1,i-NUM_in+1}(k)); 
                    aux_BSD_ROI3(i-NUM_in+1,k) = (BSD_ROI3{1,i-NUM_in+1}(k));  
                end
            end
            for k = 1:numel(BSDX_median)
    % find mean per radius bin for each case and place them in cell matrix
                aux_BSD_ROI0_2(1,k) = mean(aux_BSD_ROI0(:,k));
                aux_BSD_ROI1_2(1,k) = mean(aux_BSD_ROI1(:,k));
                aux_BSD_ROI2_2(1,k) = mean(aux_BSD_ROI2(:,k));
                aux_BSD_ROI3_2(1,k) = mean(aux_BSD_ROI3(:,k));
            end
            ROI0_idx = ROI0_idx + 1; % new ROI0 index
    % store in BSD_ROI0-9.ROI0-6 type of struct
    % use later for ensenble averages 
            BSD_avg_ROI0(ROI0_idx).ROI1 = aux_BSD_ROI0_2(:)';
            BSD_avg_ROI1(ROI0_idx).ROI1 = aux_BSD_ROI1_2(:)';   
            BSD_avg_ROI2(ROI0_idx).ROI1 = aux_BSD_ROI2_2(:)';   
            BSD_avg_ROI3(ROI0_idx).ROI1 = aux_BSD_ROI3_2(:)';  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time averaged VBSD
            for i=NUM_in:NUM
            for k = 1:numel(BSDX_median)
                aux_VBSD_ROI0_2(i-NUM_in+1,k) = (aux_VBSD_ROI0{1,i-NUM_in+1}(k)); 
                aux_VBSD_ROI1_2(i-NUM_in+1,k) = (aux_VBSD_ROI1{1,i-NUM_in+1}(k));  % place VBSD values of each image/rep. in a 2D matrix  
                aux_VBSD_ROI2_2(i-NUM_in+1,k) = (aux_VBSD_ROI2{1,i-NUM_in+1}(k)); 
                aux_VBSD_ROI3_2(i-NUM_in+1,k) = (aux_VBSD_ROI3{1,i-NUM_in+1}(k));  
            end
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:numel(BSDX_median)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% aux_VBSD_2 is equivalent to aux_BSD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find mean per radius bin for each case and place them in cell matrix
                VBSD_ROI0_avg_2(1,k) = mean(aux_VBSD_ROI0_2(:,k));
                VBSD_ROI1_avg_2(1,k) = mean(aux_VBSD_ROI1_2(:,k));
                VBSD_ROI2_avg_2(1,k) = mean(aux_VBSD_ROI2_2(:,k));
                VBSD_ROI3_avg_2(1,k) = mean(aux_VBSD_ROI3_2(:,k));
    % find time averages for segments of videos
                VBSD_ROI0_act_1(1,k) = mean(aux_VBSD_ROI0_2(1:round(numel(aux_VBSD_ROI0_2(:,k))/4),k));
                VBSD_ROI1_act_1(1,k) = mean(aux_VBSD_ROI1_2(1:round(numel(aux_VBSD_ROI1_2(:,k))/4),k));
                VBSD_ROI2_act_1(1,k) = mean(aux_VBSD_ROI2_2(1:round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_1(1,k) = mean(aux_VBSD_ROI3_2(1:round(numel(aux_VBSD_ROI3_2(:,k))/4),k));

                VBSD_ROI0_act_2(1,k) = mean(aux_VBSD_ROI0_2(round(numel(aux_VBSD_ROI0_2(:,k))/4):2*round(numel(aux_VBSD_ROI0_2(:,k))/4),k));
                VBSD_ROI1_act_2(1,k) = mean(aux_VBSD_ROI1_2(round(numel(aux_VBSD_ROI1_2(:,k))/4):2*round(numel(aux_VBSD_ROI1_2(:,k))/4),k));
                VBSD_ROI2_act_2(1,k) = mean(aux_VBSD_ROI2_2(round(numel(aux_VBSD_ROI2_2(:,k))/4):2*round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_2(1,k) = mean(aux_VBSD_ROI3_2(round(numel(aux_VBSD_ROI3_2(:,k))/4):2*round(numel(aux_VBSD_ROI3_2(:,k))/4),k));

                VBSD_ROI0_act_3(1,k) = mean(aux_VBSD_ROI0_2(2*round(numel(aux_VBSD_ROI0_2(:,k))/4):3*round(numel(aux_VBSD_ROI0_2(:,k))/4),k));
                VBSD_ROI1_act_3(1,k) = mean(aux_VBSD_ROI1_2(2*round(numel(aux_VBSD_ROI1_2(:,k))/4):3*round(numel(aux_VBSD_ROI1_2(:,k))/4),k));
                VBSD_ROI2_act_3(1,k) = mean(aux_VBSD_ROI2_2(2*round(numel(aux_VBSD_ROI2_2(:,k))/4):3*round(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_3(1,k) = mean(aux_VBSD_ROI3_2(2*round(numel(aux_VBSD_ROI3_2(:,k))/4):3*round(numel(aux_VBSD_ROI3_2(:,k))/4),k));

                VBSD_ROI0_act_4(1,k) = mean(aux_VBSD_ROI0_2(3*round(numel(aux_VBSD_ROI0_2(:,k))/4):4*(numel(aux_VBSD_ROI0_2(:,k))/4),k));
                VBSD_ROI1_act_4(1,k) = mean(aux_VBSD_ROI1_2(3*round(numel(aux_VBSD_ROI1_2(:,k))/4):4*(numel(aux_VBSD_ROI1_2(:,k))/4),k));
                VBSD_ROI2_act_4(1,k) = mean(aux_VBSD_ROI2_2(3*round(numel(aux_VBSD_ROI2_2(:,k))/4):4*(numel(aux_VBSD_ROI2_2(:,k))/4),k));
                VBSD_ROI3_act_4(1,k) = mean(aux_VBSD_ROI3_2(3*round(numel(aux_VBSD_ROI3_2(:,k))/4):4*(numel(aux_VBSD_ROI3_2(:,k))/4),k));
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
                BSD_ROI0_act_1(1,k) = mean(aux_BSD_ROI0(1:round(numel(aux_BSD_ROI0(:,k))/4),k));
                BSD_ROI1_act_1(1,k) = mean(aux_BSD_ROI1(1:round(numel(aux_BSD_ROI1(:,k))/4),k));
                BSD_ROI2_act_1(1,k) = mean(aux_BSD_ROI2(1:round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_1(1,k) = mean(aux_BSD_ROI3(1:round(numel(aux_BSD_ROI3(:,k))/4),k));


                BSD_ROI0_act_2(1,k) = mean(aux_BSD_ROI0(round(numel(aux_BSD_ROI0(:,k))/4):2*round(numel(aux_BSD_ROI0(:,k))/4),k));
                BSD_ROI1_act_2(1,k) = mean(aux_BSD_ROI1(round(numel(aux_BSD_ROI1(:,k))/4):2*round(numel(aux_BSD_ROI1(:,k))/4),k));
                BSD_ROI2_act_2(1,k) = mean(aux_BSD_ROI2(round(numel(aux_BSD_ROI2(:,k))/4):2*round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_2(1,k) = mean(aux_BSD_ROI3(round(numel(aux_BSD_ROI3(:,k))/4):2*round(numel(aux_BSD_ROI3(:,k))/4),k));

                
                BSD_ROI0_act_3(1,k) = mean(aux_BSD_ROI0(2*round(numel(aux_BSD_ROI0(:,k))/4):3*round(numel(aux_BSD_ROI0(:,k))/4),k));
                BSD_ROI1_act_3(1,k) = mean(aux_BSD_ROI1(2*round(numel(aux_BSD_ROI1(:,k))/4):3*round(numel(aux_BSD_ROI1(:,k))/4),k));
                BSD_ROI2_act_3(1,k) = mean(aux_BSD_ROI2(2*round(numel(aux_BSD_ROI2(:,k))/4):3*round(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_3(1,k) = mean(aux_BSD_ROI3(2*round(numel(aux_BSD_ROI3(:,k))/4):3*round(numel(aux_BSD_ROI3(:,k))/4),k));


                BSD_ROI0_act_4(1,k) = mean(aux_BSD_ROI0(3*round(numel(aux_BSD_ROI0(:,k))/4):4*(numel(aux_BSD_ROI0(:,k))/4),k));
                BSD_ROI1_act_4(1,k) = mean(aux_BSD_ROI1(3*round(numel(aux_BSD_ROI1(:,k))/4):4*(numel(aux_BSD_ROI1(:,k))/4),k));
                BSD_ROI2_act_4(1,k) = mean(aux_BSD_ROI2(3*round(numel(aux_BSD_ROI2(:,k))/4):4*(numel(aux_BSD_ROI2(:,k))/4),k));
                BSD_ROI3_act_4(1,k) = mean(aux_BSD_ROI3(3*round(numel(aux_BSD_ROI3(:,k))/4):4*(numel(aux_BSD_ROI3(:,k))/4),k));

            end

    % store in VBSD_ROI1-8.ROI-5 type of struct
    % use later for ensenble averages  
            VBSD_avg_ROI0(ROI0_idx).ROI0 = VBSD_ROI0_avg_2(:)';  
            VBSD_avg_ROI1(ROI0_idx).ROI0 = VBSD_ROI1_avg_2(:)';   
            VBSD_avg_ROI2(ROI0_idx).ROI0 = VBSD_ROI2_avg_2(:)';   
            VBSD_avg_ROI3(ROI0_idx).ROI0 = VBSD_ROI3_avg_2(:)';  

            VBSD_act_1_ROI0(ROI0_idx).ROI0 = VBSD_ROI0_act_1(:)';   
            VBSD_act_1_ROI1(ROI0_idx).ROI0 = VBSD_ROI1_act_1(:)'; 
            VBSD_act_1_ROI2(ROI0_idx).ROI0 = VBSD_ROI2_act_1(:)';  
            VBSD_act_1_ROI3(ROI0_idx).ROI0 = VBSD_ROI3_act_1(:)';   

            VBSD_act_2_ROI0(ROI0_idx).ROI0 = VBSD_ROI0_act_2(:)';  
            VBSD_act_2_ROI1(ROI0_idx).ROI0 = VBSD_ROI1_act_2(:)';   
            VBSD_act_2_ROI2(ROI0_idx).ROI0 = VBSD_ROI2_act_2(:)';   
            VBSD_act_2_ROI3(ROI0_idx).ROI0 = VBSD_ROI3_act_2(:)'; 

            VBSD_act_3_ROI0(ROI0_idx).ROI0 = BSD_ROI0_act_3(:)';  
            VBSD_act_3_ROI1(ROI0_idx).ROI0 = BSD_ROI1_act_3(:)';   
            VBSD_act_3_ROI2(ROI0_idx).ROI0 = BSD_ROI2_act_3(:)';   
            VBSD_act_3_ROI3(ROI0_idx).ROI0 = BSD_ROI3_act_3(:)'; 
            VBSD_act_4_ROI0(ROI0_idx).ROI0 = BSD_ROI0_act_4(:)';  
            VBSD_act_4_ROI1(ROI0_idx).ROI0 = BSD_ROI1_act_4(:)';   
            VBSD_act_4_ROI2(ROI0_idx).ROI0 = BSD_ROI2_act_4(:)';   
            VBSD_act_4_ROI3(ROI0_idx).ROI0 = BSD_ROI3_act_4(:)'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            BSD_act_1_ROI0(ROI0_idx).ROI0 = BSD_ROI0_act_1(:)';   
            BSD_act_1_ROI1(ROI0_idx).ROI0 = BSD_ROI1_act_1(:)'; 
            BSD_act_1_ROI2(ROI0_idx).ROI0 = BSD_ROI2_act_1(:)';  
            BSD_act_1_ROI3(ROI0_idx).ROI0 = BSD_ROI3_act_1(:)';  
            
            BSD_act_2_ROI0(ROI0_idx).ROI0 = BSD_ROI0_act_2(:)';  
            BSD_act_2_ROI1(ROI0_idx).ROI0 = BSD_ROI1_act_2(:)';   
            BSD_act_2_ROI2(ROI0_idx).ROI0 = BSD_ROI2_act_2(:)';   
            BSD_act_2_ROI3(ROI0_idx).ROI0 = BSD_ROI3_act_2(:)'; 

            BSD_act_3_ROI0(ROI0_idx).ROI0 = BSD_ROI0_act_3(:)';  
            BSD_act_3_ROI1(ROI0_idx).ROI0 = BSD_ROI1_act_3(:)';   
            BSD_act_3_ROI2(ROI0_idx).ROI0 = BSD_ROI2_act_3(:)';   
            BSD_act_3_ROI3(ROI0_idx).ROI0 = BSD_ROI3_act_3(:)'; 
            
            BSD_act_4_ROI0(ROI0_idx).ROI0 = BSD_ROI0_act_4(:)'; 
            BSD_act_4_ROI1(ROI0_idx).ROI0 = BSD_ROI1_act_4(:)';   
            BSD_act_4_ROI2(ROI0_idx).ROI0 = BSD_ROI2_act_4(:)';   
            BSD_act_4_ROI3(ROI0_idx).ROI0 = BSD_ROI3_act_4(:)'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store results in matrices for further analysis 
            if exist('aux_BSD_ROI0_2','var') == 1
                ROI0_mtx_idx = ROI0_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI0_mtx(ROI0_mtx_idx,k) = aux_BSD_ROI0_2(1,k);
                    VBSD_avg_ROI0_mtx(ROI0_mtx_idx,k) = VBSD_ROI0_avg_2(1,k);
                    VBSD_act_1_ROI0_mtx(ROI0_mtx_idx,k) = VBSD_ROI0_act_1(1,k);
                    VBSD_act_2_ROI0_mtx(ROI0_mtx_idx,k) = VBSD_ROI0_act_2(1,k);
                    VBSD_act_3_ROI0_mtx(ROI0_mtx_idx,k) = VBSD_ROI0_act_3(1,k);
                    VBSD_act_4_ROI0_mtx(ROI0_mtx_idx,k) = VBSD_ROI0_act_4(1,k);

                    BSD_act_1_ROI0_mtx(ROI0_mtx_idx,k) = BSD_ROI0_act_1(1,k);
                    BSD_act_2_ROI0_mtx(ROI0_mtx_idx,k) = BSD_ROI0_act_2(1,k);
                    BSD_act_3_ROI0_mtx(ROI0_mtx_idx,k) = BSD_ROI0_act_3(1,k);
                    BSD_act_4_ROI0_mtx(ROI0_mtx_idx,k) = BSD_ROI0_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI0(ROI0_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI0{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI0_mtx(ROI0_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI0_2(i-NUM_in+1,k); 
                end
            end
            end
    
            if exist('aux_BSD_ROI1_2','var') == 1
                ROI1_mtx_idx = ROI1_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI1_mtx(ROI1_mtx_idx,k) = aux_BSD_ROI1_2(1,k);
                    
                    VBSD_avg_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_avg_2(1,k);
                    
                    VBSD_act_1_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_act_1(1,k);
                    VBSD_act_2_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_act_2(1,k);
                    VBSD_act_3_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_act_3(1,k);
                    VBSD_act_4_ROI1_mtx(ROI1_mtx_idx,k) = VBSD_ROI1_act_4(1,k);

                    BSD_act_1_ROI1_mtx(ROI1_mtx_idx,k) = BSD_ROI1_act_1(1,k);
                    BSD_act_2_ROI1_mtx(ROI1_mtx_idx,k) = BSD_ROI1_act_2(1,k);
                    BSD_act_3_ROI1_mtx(ROI1_mtx_idx,k) = BSD_ROI1_act_3(1,k);
                    BSD_act_4_ROI1_mtx(ROI1_mtx_idx,k) = BSD_ROI1_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI1(ROI1_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI1{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI1_mtx(ROI1_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI1_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI2_2','var') == 1
                ROI2_mtx_idx = ROI2_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI2_mtx(ROI2_mtx_idx,k) = aux_BSD_ROI2_2(1,k);
                    VBSD_avg_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_avg_2(1,k);
                    VBSD_act_1_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_1(1,k);
                    VBSD_act_2_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_2(1,k);
                    VBSD_act_3_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_3(1,k);
                    VBSD_act_4_ROI2_mtx(ROI2_mtx_idx ,k) = VBSD_ROI2_act_4(1,k);

                    BSD_act_1_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_1(1,k);
                    BSD_act_2_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_2(1,k);
                    BSD_act_3_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_3(1,k);
                    BSD_act_4_ROI2_mtx(ROI2_mtx_idx ,k) = BSD_ROI2_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI2(ROI2_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI2{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI2_mtx(ROI2_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI2_2(i-NUM_in+1,k); 
                end
            end
            end
            if exist('aux_BSD_ROI3_2','var') == 1
                ROI3_mtx_idx = ROI3_mtx_idx + 1;
                for k = 1:numel(BSDX_median)
        % place BSD in matrix to use for ensemble avrgs
                    BSD_avg_ROI3_mtx(ROI3_mtx_idx,k) = aux_BSD_ROI3_2(1,k);
                    VBSD_avg_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_avg_2(1,k);
                    VBSD_act_1_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_1(1,k);
                    VBSD_act_2_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_2(1,k);
                    VBSD_act_3_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_3(1,k);
                    VBSD_act_4_ROI3_mtx(ROI3_mtx_idx,k) = VBSD_ROI3_act_4(1,k);

                    BSD_act_1_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_1(1,k);
                    BSD_act_2_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_2(1,k);
                    BSD_act_3_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_3(1,k);
                    BSD_act_4_ROI3_mtx(ROI3_mtx_idx,k) = BSD_ROI3_act_4(1,k);
                end
                for i = NUM_in : NUM
                    alpha_ins_ROI3(ROI3_mtx_idx,i-NUM_in+1) = aux_alpha_ins_ROI3{i-NUM_in+1};
                end
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
            for k = 1 : numel(BSDX_median)
                for i = NUM_in :NUM
                   VBSD_inst_ROI3_mtx(ROI3_mtx_idx,i-NUM_in+1,k) = aux_VBSD_ROI3_2(i-NUM_in+1,k); 
                end
            end
            end

      % save results
            cd(dir_save)
            save('BSDX_median.mat','BSDX_median','-v7.3');
            save(['BSD_ins_ROI0_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI0','-v7.3');  
            save(['BSD_ins_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI1','-v7.3');   
            save(['BSD_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI2','-v7.3'); 
            save(['BSD_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_BSD_ROI3','-v7.3'); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save results
            save(['VBSD_ins_ROI0_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI0_2','-v7.3'); 
            save(['VBSD_ins_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI1_2','-v7.3');   
            save(['VBSD_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI2_2','-v7.3');  
            save(['VBSD_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'aux_VBSD_ROI3_2','-v7.3'); 

    %     save(['VBSD_avg_ROI0_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI0_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_avg_2','-v7.3');   
    %     save(['VBSD_avg_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_avg_2','-v7.3');     
    % 
    %     save(['VBSD_act_1_ROI0_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI0_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_1','-v7.3');   
    %     save(['VBSD_act_1_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_1','-v7.3');
    % 
    %     save(['VBSD_act_2_ROI0_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI0_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_2','-v7.3');   
    %     save(['VBSD_act_2_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_2','-v7.3');
    % 
    %     save(['VBSD_act_3_ROI0_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI0_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_3','-v7.3');   
    %     save(['VBSD_act_3_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_3','-v7.3');
    %    
    %     save(['VBSD_act_4_ROI0_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI4_act_0','-v7.3');   
    %     save(['VBSD_act_4_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI1_act_4','-v7.3');   
    %     save(['VBSD_act_4_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI2_act_4','-v7.3'); 
    %     save(['VBSD_act_4_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'VBSD_ROI3_act_4','-v7.3');

            save(['alpha_ins_ROI0_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI0','-v7.3');  
            save(['alpha_ins_ROI1_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI1','-v7.3');  
            save(['alpha_ins_ROI2_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI2','-v7.3');  
            save(['alpha_ins_ROI3_',wave,'_',phase{case_idx},'_',ROI{case_idx},'_',run{case_idx},'.mat'],'alpha_ins_ROI3','-v7.3'); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








        else
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------END OF TIME AVERAGED ANALYSIS----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------ENSEMBLE AVERAGE STATISTICS------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find average over NUM_case repetitions for k bin radii 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if avg == 1
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% FOR ANALYTICAL METHOD OF DERIVING BSD 
    for k = 1: numel(BSDX_median) % number of bins
% BSD analytical 
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI0_mtx','var') == 1
            BSD_mean_ROI0(k) = mean(BSD_avg_ROI0_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI0(k) = mean(BSD_act_1_ROI0_mtx(:,k)); 
            BSD_mean_act_2_ROI0(k) = mean(BSD_act_2_ROI0_mtx(:,k)); 
            BSD_mean_act_3_ROI0(k) = mean(BSD_act_3_ROI0_mtx(:,k)); 
            BSD_mean_act_4_ROI0(k) = mean(BSD_act_4_ROI0_mtx(:,k)); 
            
            VBSD_mean_ROI0(k) = mean(VBSD_avg_ROI0_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI0(k) = mean(VBSD_act_1_ROI0_mtx(:,k)); 
            VBSD_mean_act_2_ROI0(k) = mean(VBSD_act_2_ROI0_mtx(:,k)); 
            VBSD_mean_act_3_ROI0(k) = mean(VBSD_act_3_ROI0_mtx(:,k)); 
            VBSD_mean_act_4_ROI0(k) = mean(VBSD_act_4_ROI0_mtx(:,k)); 
            
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
%            VBSD_inst_ROI1_matrix = VBSD_inst_ROI1_mtx; 
            
        else 
            disp('only one repetition for this ROI');
        end

        if exist('BSD_avg_ROI1_mtx','var') == 1
            BSD_mean_ROI1(k) = mean(BSD_avg_ROI1_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI1(k) = mean(BSD_act_1_ROI1_mtx(:,k)); 
            BSD_mean_act_2_ROI1(k) = mean(BSD_act_2_ROI1_mtx(:,k)); 
            BSD_mean_act_3_ROI1(k) = mean(BSD_act_3_ROI1_mtx(:,k)); 
            BSD_mean_act_4_ROI1(k) = mean(BSD_act_4_ROI1_mtx(:,k)); 
            
            VBSD_mean_ROI1(k) = mean(VBSD_avg_ROI1_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI1(k) = mean(VBSD_act_1_ROI1_mtx(:,k)); 
            VBSD_mean_act_2_ROI1(k) = mean(VBSD_act_2_ROI1_mtx(:,k)); 
            VBSD_mean_act_3_ROI1(k) = mean(VBSD_act_3_ROI1_mtx(:,k)); 
            VBSD_mean_act_4_ROI1(k) = mean(VBSD_act_4_ROI1_mtx(:,k)); 
            
    %% place instant BSD/VBSD in matrix to use for instanteneous ensemble avrgs 
%            VBSD_inst_ROI1_matrix = VBSD_inst_ROI1_mtx; 
            
        else 
            disp('only one repetition for this ROI');
        end
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI2_mtx','var') == 1
            BSD_mean_ROI2(k) = mean(BSD_avg_ROI2_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI2(k) = mean(BSD_act_1_ROI2_mtx(:,k)); 
            BSD_mean_act_2_ROI2(k) = mean(BSD_act_2_ROI2_mtx(:,k)); 
            BSD_mean_act_3_ROI2(k) = mean(BSD_act_3_ROI2_mtx(:,k)); 
            BSD_mean_act_4_ROI2(k) = mean(BSD_act_4_ROI2_mtx(:,k)); 
            
            VBSD_mean_ROI2(k) = mean(VBSD_avg_ROI2_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI2(k) = mean(VBSD_act_1_ROI2_mtx(:,k)); 
            VBSD_mean_act_2_ROI2(k) = mean(VBSD_act_2_ROI2_mtx(:,k)); 
            VBSD_mean_act_3_ROI2(k) = mean(VBSD_act_3_ROI2_mtx(:,k)); 
            VBSD_mean_act_4_ROI2(k) = mean(VBSD_act_4_ROI2_mtx(:,k));   
        else 
            disp('only one repetition for this ROI');
        end
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI3_mtx','var') == 1
            BSD_mean_ROI3(k) = mean(BSD_avg_ROI3_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI3(k) = mean(BSD_act_1_ROI3_mtx(:,k)); 
            BSD_mean_act_2_ROI3(k) = mean(BSD_act_2_ROI3_mtx(:,k)); 
            BSD_mean_act_3_ROI3(k) = mean(BSD_act_3_ROI3_mtx(:,k)); 
            BSD_mean_act_4_ROI3(k) = mean(BSD_act_4_ROI3_mtx(:,k));
            
            VBSD_mean_ROI3(k) = mean(VBSD_avg_ROI3_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI3(k) = mean(VBSD_act_1_ROI3_mtx(:,k)); 
            VBSD_mean_act_2_ROI3(k) = mean(VBSD_act_2_ROI3_mtx(:,k)); 
            VBSD_mean_act_3_ROI3(k) = mean(VBSD_act_3_ROI3_mtx(:,k)); 
            VBSD_mean_act_4_ROI3(k) = mean(VBSD_act_4_ROI3_mtx(:,k));
        else 
            disp('only one repetition for this ROI');
        end
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI4_mtx','var') == 1
            BSD_mean_ROI4(k) = mean(BSD_avg_ROI4_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI4(k) = mean(BSD_act_1_ROI4_mtx(:,k)); 
            BSD_mean_act_2_ROI4(k) = mean(BSD_act_2_ROI4_mtx(:,k)); 
            BSD_mean_act_3_ROI4(k) = mean(BSD_act_3_ROI4_mtx(:,k)); 
            BSD_mean_act_4_ROI4(k) = mean(BSD_act_4_ROI4_mtx(:,k));
            
            VBSD_mean_ROI4(k) = mean(VBSD_avg_ROI4_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI4(k) = mean(VBSD_act_1_ROI4_mtx(:,k)); 
            VBSD_mean_act_2_ROI4(k) = mean(VBSD_act_2_ROI4_mtx(:,k)); 
            VBSD_mean_act_3_ROI4(k) = mean(VBSD_act_3_ROI4_mtx(:,k)); 
            VBSD_mean_act_4_ROI4(k) = mean(VBSD_act_4_ROI4_mtx(:,k));
        else 
            disp('only one repetition for this ROI');
        end
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI5_mtx','var') == 1
            BSD_mean_ROI5(k) = mean(BSD_avg_ROI5_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI5(k) = mean(BSD_act_1_ROI5_mtx(:,k)); 
            BSD_mean_act_2_ROI5(k) = mean(BSD_act_2_ROI5_mtx(:,k)); 
            BSD_mean_act_3_ROI5(k) = mean(BSD_act_3_ROI5_mtx(:,k)); 
            BSD_mean_act_4_ROI5(k) = mean(BSD_act_4_ROI5_mtx(:,k));
            
            VBSD_mean_ROI5(k) = mean(VBSD_avg_ROI5_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI5(k) = mean(VBSD_act_1_ROI5_mtx(:,k)); 
            VBSD_mean_act_2_ROI5(k) = mean(VBSD_act_2_ROI5_mtx(:,k)); 
            VBSD_mean_act_3_ROI5(k) = mean(VBSD_act_3_ROI5_mtx(:,k)); 
            VBSD_mean_act_4_ROI5(k) = mean(VBSD_act_4_ROI5_mtx(:,k));
        else 
            disp('only one repetition for this ROI');
        end
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI6_mtx','var') == 1
            BSD_mean_ROI6(k) = mean(BSD_avg_ROI6_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI6(k) = mean(BSD_act_1_ROI6_mtx(:,k)); 
            BSD_mean_act_2_ROI6(k) = mean(BSD_act_2_ROI6_mtx(:,k)); 
            BSD_mean_act_3_ROI6(k) = mean(BSD_act_3_ROI6_mtx(:,k)); 
            BSD_mean_act_4_ROI6(k) = mean(BSD_act_4_ROI6_mtx(:,k));
            
            VBSD_mean_ROI6(k) = mean(VBSD_avg_ROI6_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI6(k) = mean(VBSD_act_1_ROI6_mtx(:,k)); 
            VBSD_mean_act_2_ROI6(k) = mean(VBSD_act_2_ROI6_mtx(:,k)); 
            VBSD_mean_act_3_ROI6(k) = mean(VBSD_act_3_ROI6_mtx(:,k)); 
            VBSD_mean_act_4_ROI6(k) = mean(VBSD_act_4_ROI6_mtx(:,k)); 
        else 
            disp('only one repetition for this ROI');
        end
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI7_mtx','var') == 1
            BSD_mean_ROI7(k) = mean(BSD_avg_ROI7_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI7(k) = mean(BSD_act_1_ROI7_mtx(:,k)); 
            BSD_mean_act_2_ROI7(k) = mean(BSD_act_2_ROI7_mtx(:,k)); 
            BSD_mean_act_3_ROI7(k) = mean(BSD_act_3_ROI7_mtx(:,k)); 
            BSD_mean_act_4_ROI7(k) = mean(BSD_act_4_ROI7_mtx(:,k));
            
            VBSD_mean_ROI7(k) = mean(VBSD_avg_ROI7_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI7(k) = mean(VBSD_act_1_ROI7_mtx(:,k)); 
            VBSD_mean_act_2_ROI7(k) = mean(VBSD_act_2_ROI7_mtx(:,k)); 
            VBSD_mean_act_3_ROI7(k) = mean(VBSD_act_3_ROI7_mtx(:,k)); 
            VBSD_mean_act_4_ROI7(k) = mean(VBSD_act_4_ROI7_mtx(:,k));
        else 
            disp('only one repetition for this ROI');
        end
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI8_mtx','var') == 1
            BSD_mean_ROI8(k) = mean(BSD_avg_ROI8_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI8(k) = mean(BSD_act_1_ROI8_mtx(:,k)); 
            BSD_mean_act_2_ROI8(k) = mean(BSD_act_2_ROI8_mtx(:,k)); 
            BSD_mean_act_3_ROI8(k) = mean(BSD_act_3_ROI8_mtx(:,k)); 
            BSD_mean_act_4_ROI8(k) = mean(BSD_act_4_ROI8_mtx(:,k));
            
            VBSD_mean_ROI8(k) = mean(VBSD_avg_ROI8_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI8(k) = mean(VBSD_act_1_ROI8_mtx(:,k)); 
            VBSD_mean_act_2_ROI8(k) = mean(VBSD_act_2_ROI8_mtx(:,k)); 
            VBSD_mean_act_3_ROI8(k) = mean(VBSD_act_3_ROI8_mtx(:,k)); 
            VBSD_mean_act_4_ROI8(k) = mean(VBSD_act_4_ROI8_mtx(:,k));
        else 
            disp('only one repetition for this ROI');
        end
% check if there are repetitions over which to average
        if exist('BSD_avg_ROI9_mtx','var') == 1
            BSD_mean_ROI9(k) = mean(BSD_avg_ROI9_mtx(:,k)); % ensemble average of BSD 
            
            BSD_mean_act_1_ROI9(k) = mean(BSD_act_1_ROI9_mtx(:,k)); 
            BSD_mean_act_2_ROI9(k) = mean(BSD_act_2_ROI9_mtx(:,k)); 
            BSD_mean_act_3_ROI9(k) = mean(BSD_act_3_ROI9_mtx(:,k)); 
            BSD_mean_act_4_ROI9(k) = mean(BSD_act_4_ROI9_mtx(:,k));
            
            VBSD_mean_ROI9(k) = mean(VBSD_avg_ROI9_mtx(:,k)); 
            
            VBSD_mean_act_1_ROI9(k) = mean(VBSD_act_1_ROI9_mtx(:,k)); 
            VBSD_mean_act_2_ROI9(k) = mean(VBSD_act_2_ROI9_mtx(:,k)); 
            VBSD_mean_act_3_ROI9(k) = mean(VBSD_act_3_ROI9_mtx(:,k)); 
            VBSD_mean_act_4_ROI9(k) = mean(VBSD_act_4_ROI9_mtx(:,k));
        else 
            disp('only one repetition for this ROI');
        end
    end
    for i = NUM_in:NUM % time count
%% find ensemble averages of instanteneous void fraction a(x,t) and
%% instanteneous Volume Size Distribution
        if exist('BSD_avg_ROI0_mtx','var') == 1
            alpha_mean_ROI0(i-NUM_in+1) = mean(alpha_ins_ROI0(:,i-NUM_in+1));
            if numel(VBSD_inst_ROI0_mtx(1,:,:)) == 1
                for k = 1 : numel(BSDX_median)
                    VBSD_inst_avg_ROI0(i-NUM_in+1,k) = mean(VBSD_inst_ROI0_mtx(:,i-NUM_in+1,k));
                end
            else
                for k = 1 : numel(BSDX_median)
                    VBSD_inst_avg_ROI0(i-NUM_in+1,k) = mean(VBSD_inst_ROI0_mtx(:,i-NUM_in+1,k));
                end
            end
        end

        if exist('BSD_avg_ROI1_mtx','var') == 1
            alpha_mean_ROI1(i-NUM_in+1) = mean(alpha_ins_ROI1(:,i-NUM_in+1));
            if numel(VBSD_inst_ROI1_mtx(1,:,:)) == 1
                for k = 1 : numel(BSDX_median)
                    VBSD_inst_avg_ROI1(i-NUM_in+1,k) = mean(VBSD_inst_ROI1_mtx(:,i-NUM_in+1,k));
                end
            else
                for k = 1 : numel(BSDX_median)
                    VBSD_inst_avg_ROI1(i-NUM_in+1,k) = mean(VBSD_inst_ROI1_mtx(:,i-NUM_in+1,k));
                end
            end
        end
        if exist('BSD_avg_ROI2_mtx','var') == 1
            alpha_mean_ROI2(i-NUM_in+1) = mean(alpha_ins_ROI2(:,i-NUM_in+1));
            for k = 1 : numel(BSDX_median)
               VBSD_inst_avg_ROI2(i-NUM_in+1,k) = mean(VBSD_inst_ROI2_mtx(:,i-NUM_in+1,k));
            end
        end
        if exist('BSD_avg_ROI3_mtx','var') == 1
            alpha_mean_ROI3(i-NUM_in+1) = mean(alpha_ins_ROI3(:,i-NUM_in+1));
            for k = 1 : numel(BSDX_median)
               VBSD_inst_avg_ROI3(i-NUM_in+1,k) = mean(VBSD_inst_ROI3_mtx(:,i-NUM_in+1,k));
            end
        end
        if exist('BSD_avg_ROI4_mtx','var') == 1
            alpha_mean_ROI4(i-NUM_in+1) = mean(alpha_ins_ROI4(:,i-NUM_in+1));
            for k = 1 : numel(BSDX_median)
               VBSD_inst_avg_ROI4(i-NUM_in+1,k) = mean(VBSD_inst_ROI4_mtx(:,i-NUM_in+1,k));
            end
        end
        if exist('BSD_avg_ROI5_mtx','var') == 1
            alpha_mean_ROI5(i-NUM_in+1) = mean(alpha_ins_ROI5(:,i-NUM_in+1));
            for k = 1 : numel(BSDX_median)
               VBSD_inst_avg_ROI5(i-NUM_in+1,k) = mean(VBSD_inst_ROI5_mtx(:,i-NUM_in+1,k));
            end
        end
        if exist('BSD_avg_ROI6_mtx','var') == 1
            alpha_mean_ROI6(i-NUM_in+1) = mean(alpha_ins_ROI6(:,i-NUM_in+1));
            for k = 1 : numel(BSDX_median)
               VBSD_inst_avg_ROI6(i-NUM_in+1,k) = mean(VBSD_inst_ROI6_mtx(:,i-NUM_in+1,k));
            end
        end
        if exist('BSD_avg_ROI7_mtx','var') == 1
            alpha_mean_ROI7(i-NUM_in+1) = mean(alpha_ins_ROI7(:,i-NUM_in+1));
            for k = 1 : numel(BSDX_median)
               VBSD_inst_avg_ROI7(i-NUM_in+1,k) = mean(VBSD_inst_ROI7_mtx(:,i-NUM_in+1,k));
            end
        end
        if exist('BSD_avg_ROI8_mtx','var') == 1
            alpha_mean_ROI8(i-NUM_in+1) = mean(alpha_ins_ROI8(:,i-NUM_in+1));
            if numel(VBSD_inst_ROI8_mtx(1,:,:)) == 1
                for k = 1 : numel(BSDX_median)
                   VBSD_inst_avg_ROI8(i-NUM_in+1,k) = (VBSD_inst_ROI8_mtx(:,i-NUM_in+1,k));
                end
            else
                for k = 1 : numel(BSDX_median)
                   VBSD_inst_avg_ROI8(i-NUM_in+1,k) = mean(VBSD_inst_ROI8_mtx(:,i-NUM_in+1,k));
                end
            end
        end
        if exist('BSD_avg_ROI9_mtx','var') == 1
            alpha_mean_ROI9(i-NUM_in+1) = mean(alpha_ins_ROI9(:,i-NUM_in+1));
            if numel(VBSD_inst_ROI9_mtx(1,:,:)) == 1
                for k = 1 : numel(BSDX_median)
                   VBSD_inst_avg_ROI9(i-NUM_in+1,k) = (VBSD_inst_ROI9_mtx(:,i-NUM_in+1,k));
                end
            else
                for k = 1 : numel(BSDX_median)
                   VBSD_inst_avg_ROI9(i-NUM_in+1,k) = mean(VBSD_inst_ROI9_mtx(:,i-NUM_in+1,k));
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(dir_save_rep);
    if exist('BSD_mean_ROI0','var') == 1
        save(['BSD_mean_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI0','-v7.3');  
        
        save(['BSD_mean_act_1_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI0','-v7.3');   
        save(['BSD_mean_act_2_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI0','-v7.3');   
        save(['BSD_mean_act_3_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI0','-v7.3');   
        save(['BSD_mean_act_4_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI0','-v7.3'); 
        
        save(['VBSD_mean_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI0','-v7.3');   
        save(['VBSD_mean_act_1_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI0','-v7.3');   
        save(['VBSD_mean_act_2_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI0','-v7.3');   
        save(['VBSD_mean_act_3_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI0','-v7.3');   
        save(['VBSD_mean_act_4_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI0','-v7.3'); 
        save(['alpha_mean_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI0','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI0_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI0','-v7.3'); 
% instanteneous volume all runs:        
        save(['VBSD_inst_ROI0_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI0_mtx','-v7.3'); 
    end
    
    if exist('BSD_mean_ROI1','var') == 1
        save(['BSD_mean_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI1','-v7.3');  
        
        save(['BSD_mean_act_1_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI1','-v7.3');   
        save(['BSD_mean_act_2_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI1','-v7.3');   
        save(['BSD_mean_act_3_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI1','-v7.3');   
        save(['BSD_mean_act_4_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI1','-v7.3'); 
        
        save(['VBSD_mean_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI1','-v7.3');   
        save(['VBSD_mean_act_1_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI1','-v7.3');   
        save(['VBSD_mean_act_2_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI1','-v7.3');   
        save(['VBSD_mean_act_3_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI1','-v7.3');   
        save(['VBSD_mean_act_4_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI1','-v7.3'); 
        save(['alpha_mean_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI1','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI1_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI1','-v7.3'); 
% instanteneous volume all runs:        
        save(['VBSD_inst_ROI1_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI1_mtx','-v7.3'); 
    end
    if exist('BSD_mean_ROI2','var') == 1
        save(['BSD_mean_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI2','-v7.3');  
        save(['BSD_mean_act_1_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI2','-v7.3');   
        save(['BSD_mean_act_2_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI2','-v7.3');   
        save(['BSD_mean_act_3_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI2','-v7.3');   
        save(['BSD_mean_act_4_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI2','-v7.3'); 
        
        save(['VBSD_mean_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI2','-v7.3');   
        save(['VBSD_mean_act_1_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI2','-v7.3');   
        save(['VBSD_mean_act_2_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI2','-v7.3');   
        save(['VBSD_mean_act_3_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI2','-v7.3');   
        save(['VBSD_mean_act_4_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI2','-v7.3');  
        save(['alpha_mean_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI2','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI2_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI2','-v7.3'); 
% instanteneous volume all runs:        
        save(['VBSD_inst_ROI2_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI2_mtx','-v7.3'); 
    end
    if exist('BSD_mean_ROI3','var') == 1
        save(['BSD_mean_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI3','-v7.3');   
        save(['BSD_mean_act_1_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI3','-v7.3');   
        save(['BSD_mean_act_2_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI3','-v7.3');   
        save(['BSD_mean_act_3_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI3','-v7.3');   
        save(['BSD_mean_act_4_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI3','-v7.3'); 
        
        save(['VBSD_mean_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI3','-v7.3');   
        save(['VBSD_mean_act_1_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI3','-v7.3');   
        save(['VBSD_mean_act_2_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI3','-v7.3');   
        save(['VBSD_mean_act_3_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI3','-v7.3');   
        save(['VBSD_mean_act_4_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI3','-v7.3');   
        save(['alpha_mean_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI3','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI3_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI3','-v7.3'); 
% instanteneous volume all runs:        
        save(['VBSD_inst_ROI3_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI3_mtx','-v7.3'); 
    end
    if exist('BSD_mean_ROI4','var') == 1
        save(['BSD_mean_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI4','-v7.3'); 
        save(['BSD_mean_act_1_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI4','-v7.3');   
        save(['BSD_mean_act_2_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI4','-v7.3');   
        save(['BSD_mean_act_3_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI4','-v7.3');   
        save(['BSD_mean_act_4_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI4','-v7.3');
        
        save(['VBSD_mean_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI4','-v7.3');   
        save(['VBSD_mean_act_1_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI4','-v7.3');   
        save(['VBSD_mean_act_2_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI4','-v7.3');   
        save(['VBSD_mean_act_3_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI4','-v7.3');   
        save(['VBSD_mean_act_4_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI4','-v7.3');  
        save(['alpha_mean_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI4','-v7.3'); 
% instanteneous ensemble average volume per radius:
       save(['VBSD_inst_avg_ROI4_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI4','-v7.3'); 
% instanteneous volume all runs:        
        save(['VBSD_inst_ROI4_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI4_mtx','-v7.3'); 
    end
    if exist('BSD_mean_ROI5','var') == 1
        save(['BSD_mean_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI5','-v7.3');  
        save(['BSD_mean_act_1_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI5','-v7.3');   
        save(['BSD_mean_act_2_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI5','-v7.3');   
        save(['BSD_mean_act_3_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI5','-v7.3');   
        save(['BSD_mean_act_4_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI5','-v7.3'); 
        
        save(['VBSD_mean_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI5','-v7.3');   
        save(['VBSD_mean_act_1_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI5','-v7.3');   
        save(['VBSD_mean_act_2_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI5','-v7.3');   
        save(['VBSD_mean_act_3_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI5','-v7.3');   
        save(['VBSD_mean_act_4_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI5','-v7.3'); 
        save(['alpha_mean_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI5','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI5_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI5','-v7.3'); 
% instanteneous volume all runs:        
        save(['VBSD_inst_ROI5_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI5_mtx','-v7.3'); 
    end
    if exist('BSD_mean_ROI6','var') == 1
        save(['BSD_mean_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI6','-v7.3'); 
        save(['BSD_mean_act_1_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI6','-v7.3');   
        save(['BSD_mean_act_2_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI6','-v7.3');   
        save(['BSD_mean_act_3_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI6','-v7.3');   
        save(['BSD_mean_act_4_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI6','-v7.3'); 
        
        save(['VBSD_mean_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI6','-v7.3');   
        save(['VBSD_mean_act_1_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI6','-v7.3');   
        save(['VBSD_mean_act_2_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI6','-v7.3');   
        save(['VBSD_mean_act_3_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI6','-v7.3');   
        save(['VBSD_mean_act_4_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI6','-v7.3');   
        save(['alpha_mean_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI6','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI6_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI6','-v7.3'); 
% instanteneous volume all runs:        
        save(['VBSD_inst_ROI6_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI6_mtx','-v7.3'); 
    end
    if exist('BSD_mean_ROI7','var') == 1
        save(['BSD_mean_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI7','-v7.3');  
        save(['BSD_mean_act_1_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI7','-v7.3');   
        save(['BSD_mean_act_2_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI7','-v7.3');   
        save(['BSD_mean_act_3_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI7','-v7.3');   
        save(['BSD_mean_act_4_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI7','-v7.3'); 
        
        save(['VBSD_mean_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI7','-v7.3');   
        save(['VBSD_mean_act_1_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI7','-v7.3');   
        save(['VBSD_mean_act_2_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI7','-v7.3');   
        save(['VBSD_mean_act_3_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI7','-v7.3');   
        save(['VBSD_mean_act_4_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI7','-v7.3');  
        save(['alpha_mean_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI7','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI7_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI7','-v7.3'); 
% instanteneous volume all runs:        
        save(['VBSD_inst_ROI7_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI7_mtx','-v7.3'); 
    end
    if exist('BSD_mean_ROI8','var') == 1
        save(['BSD_mean_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI8','-v7.3'); 
        save(['BSD_mean_act_1_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI8','-v7.3');   
        save(['BSD_mean_act_2_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI8','-v7.3');   
        save(['BSD_mean_act_3_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI8','-v7.3');   
        save(['BSD_mean_act_4_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI8','-v7.3'); 
        
        save(['VBSD_mean_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI8','-v7.3');   
        save(['VBSD_mean_act_1_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI8','-v7.3');   
        save(['VBSD_mean_act_2_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI8','-v7.3');   
        save(['VBSD_mean_act_3_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI8','-v7.3');   
        save(['VBSD_mean_act_4_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI8','-v7.3');  
        save(['alpha_mean_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI8','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI8_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI8','-v7.3'); 
% instanteneous volume  all runs:        
        save(['VBSD_inst_ROI8_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI8_mtx','-v7.3'); 
    end
    if exist('BSD_mean_ROI9','var') == 1
        save(['BSD_mean_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_ROI9','-v7.3'); 
        save(['BSD_mean_act_1_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_1_ROI9','-v7.3');   
        save(['BSD_mean_act_2_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_2_ROI9','-v7.3');   
        save(['BSD_mean_act_3_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_3_ROI9','-v7.3');   
        save(['BSD_mean_act_4_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'BSD_mean_act_4_ROI9','-v7.3'); 
        
        save(['VBSD_mean_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_ROI9','-v7.3');   
        save(['VBSD_mean_act_1_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_1_ROI9','-v7.3');   
        save(['VBSD_mean_act_2_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_2_ROI9','-v7.3');   
        save(['VBSD_mean_act_3_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_3_ROI9','-v7.3');   
        save(['VBSD_mean_act_4_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_mean_act_4_ROI9','-v7.3');  
        save(['alpha_mean_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'alpha_mean_ROI9','-v7.3'); 
% instanteneous ensemble average volume per radius:
        save(['VBSD_inst_avg_ROI9_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_avg_ROI9','-v7.3'); 
% instanteneous volume  all runs:        
        save(['VBSD_inst_ROI9_mtx_',wave,'_',phase{case_idx},'_',comments,'_.mat'],'VBSD_inst_ROI9_mtx','-v7.3'); 
    end
    
cd(dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif avg == 0

    disp('goodbye');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










