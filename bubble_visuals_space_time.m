%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------------BUBBLE VISUALS --------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visualisation obtained from bubble_space/time_ANALYSIS_freshwater.m
%%% load data and ask for ploting options.
%%% some calculations are made within the script: some max,mins,statistical
%%% analysis and averages.
%%% plot steady BSD for one wave, intanteneous or averaged for many waves
%%% see bubble_analysis.m for more details on derivations of bubble statistics 
%%% 
clear;
clc;
close;
beep off;
set(0,'DefaultFigureVisible','on'); % graphics/plots on
% parameter alocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 0;   % (0) windows or (1) linux
disp('P?');
%% ----------------------------PLOT OPTIONS--------------------------------
% P = 0 ;
% P = 0; % plots for paper 0
% P = 2; % plots for seminar
P = 3; % plots for paper 1
%% P = 1: volume of bubbles and void fraction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1: volume of bubbles time series and void fraction time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NUM_in = 1;   % *first image  index (check file names in folder)*
NUM = 7500;    %   last image index
% disp('how many repetitions or cases?');
% NUM_case = input('');
NUM_case = 1;  % choose number of cases (repetitions, phase shifts, amps, spectra)
NUM_wave = 16; % total nunber of wave cases
markers = {'d','x','+','o'};
%% ------------------------------Parameters--------------------------------
data_type = 'freshwater';
wave_dir_load = ["GW_175" "GW_18" "JS_15" "JS_155"];
wave = {'GW_175','GW_18','JS_15','JS_155'};
im_name = cell(1,NUM);
run = {'run1','run2','run3'};
run_max = numel(run);
phase = {'peak','trough','minp2','posp2'};
ROI = {'ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI7','ROI8'};
% ROI = {'ROI0','ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI7','ROI8','ROI9'};
ROI_max = numel(ROI);
ROI_old = {'ROI1','ROI2','ROI3','ROI4','ROI5','ROI5'};
% ROI_old = {'ROI0','ROI1','ROI2','ROI3','ROI4','ROI5','ROI6'};
dir_general = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\',data_type];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
per_bin = (2.4*10.^-1);   % norm factor radii bin in mm
% per_bin = 1;  % norm factor not used
vol_ROI = 350.*50.*100;  % volume of each rect ROI in mm^3;
vol_fraction = 350.*50.*20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vol_ROI_mm = vol_ROI.*10^9; % volume of each ROI in mm^3
vol_cube = 10^9;    % 1m^3 in mm^3
vol_factor = vol_cube./vol_ROI; % N is divided by vol_roi/vol_cube so--> N(r)*vol_factor
norm_bin_size = 1./per_bin;   % norm factor volume: N is multiplied with this as well
% .*vol_factor.*norm_bin_size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_p = 1 ; % period of central frequency in sec
lambda= 1450; % wavelength of central frequency in mm
U_c = 1.46;   % phase speed calculated by newton algorithm in dispersion.m
frame_rate = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%% wave parameters 
bulge_area = [334,60,260,170;540,233,480,270,;100,23,50,15;176,22,145,65]; % max bulge area for each case (calculated from wave crests visuals/averages)
% bulge_area = [318,50,240,160;540,255,616,181,;103,10,33,15;144,24,146,72]; % some newer calculations (very similar tho)

% bulge_area = [600,60,260,170;540,233,480,170,;100,22,50,15;176,22,145,65]; % *** not used **8max bulge area for each case (calculated from wave crests and mainprogram.m for spillers)

for wave_idx = 1:4
    for phase_idx = 1:4
        vol_bulge.(wave{wave_idx}).(phase{phase_idx})= 350*bulge_area(wave_idx,phase_idx); % average calculated volume of each wave bulge
    end
end
% these are calculated in wave_crest_visuals.m
%--------------------------------------------------------------------------
gamma_error = [0.13 0.05 0.15 0.13;0.1 0.26 0.06 0.2 ; 0.03 0.05 0.02 0.1;0.08 0.16 0.11 0.13]; 
gamma_error_mean = mean(gamma_error(:));
%--------------------------------------------------------------------------
vol_bulge_error = 350.*[55 10 50 45; 17 50 42 55; .5 6 3.5 3.5; 22 9 9.9 50];
bulge_vol = 350.*bulge_area; % same as vol_bulge but in a matrix
%--------------------------------------------------------------------------
vol_bulge_error_mean = mean(vol_bulge_error(:));
%--------------------------------------------------------------------------
break_loc = {189./lambda,132./lambda,158./lambda,239./lambda;161./lambda,170./lambda,128./lambda,231./lambda;200./lambda,190./lambda,211./lambda,225./lambda;204./lambda,187./lambda,183./lambda,230./lambda}; % non-dimensional L (breaking distance from x=0) 
% break_loc = break_loc';
% break_dist = [224,144,200,256;196,75,146,269;234,263,301,416;278,237,256,320]; % in mm
%--------------------------------------------------------------------------
bulge_length = [26,12.5,23,20;34,22.6,31.5,19.5;15,9.8,13,8;22,6.6,20.4,14.5]; % in mm 
%--------------------------------------------------------------------------
bulge_height = [19,8.2,23.7,14;26.7,19.5,27.1,16.1;10.3,3.43,4.6,3;16,4.5,14.6,7.5]; % in mm 
%--------------------------------------------------------------------------
bulge_length_error = [2.6,.7,5.8,3.1;1.2,2.3,4.5,1.1;1.2,2.3,0.9,1.6;0.9,3,0.8,6];
bulge_length_error_mean = mean(bulge_length_error(:));
%--------------------------------------------------------------------------
bulge_height_error = [4.6,0.6,10,2;2.4,7,5.5,3.8;0.3,2.4,0.5,0.6;1.8,0.4,1.5,6];
%--------------------------------------------------------------------------
timeorig = [63.738,62.892,63.272,64.256;63.674,62.866,63.212,64.154;63.826,63.294,63.536,64.18;63.808,63.284,63.51,64.114]; % in seconds
%--------------------------------------------------------------------------
deltat = [0.086,0.066,0.1,0.078;0.92,0.042,0.092,0.09;0.042,0.022,0.016,0.044;0.08,0.028,0.074,0.094];
%--------------------------------------------------------------------------
gamma = bulge_height./bulge_length; % dz/dx at breaking 
%--------------------------------------------------------------------------
beta = 0.4.*(gamma - 0.08).^(5/2); % breaking strength 
%--------------------------------------------------------------------------
beta_error =  0.4.*(gamma_error).^(5/2);
beta_error_mean = mean(beta_error(:));
%--------------------------------------------------------------------------
hsq = (pi.*(bulge_height.*bulge_length)./4); % pi*h^2/4 elliptical cross section area of entrained air 
%--------------------------------------------------------------------------
hsq_error = pi.*(bulge_length_error.* bulge_height_error)./4 ;
%--------------------------------------------------------------------------
hsq_vol = hsq.*350; % *width of flume
%--------------------------------------------------------------------------
hsq_vol_error = hsq_error.*350;
%--------------------------------------------------------------------------
crest_v = [1320,1220,1320,1130;1300,1270,1520,1330;1130,1010,1180,1200;1270,980,1100,1180]; % crest velocity in mm/s , at breaking calculated from crestsmainprogram_velocities.m
%--------------------------------------------------------------------------
rho = 0.001; % in g/mm^3 
g = 981; % in mm/s^2
epsilon_l = (rho/g).*beta.*(crest_v.^5); % energy diss rate in kg/mm^3* mm/s^2
%--------------------------------------------------------------------------
z_max_u = [93,85.3,94.7,87.8;92.3,88,94.1,86;75.1,66.4,74.5,73.2;74.8,61.9,74.8,73]; %max surf elevation upstream breaking
z_max_d = [81.9,72,80.5,76.8;83.7,75.4,83.8,85;74.3,58.8,67.7,72.8;74.5,56,70.4,71.6]; %max surf elevation downstream breaking
%--------------------------------------------------------------------------
We = 4.7; sigma = 73;
a_h = 2^(-8/5).*(sigma.*We/rho)^(3/5).*epsilon_l.^(-2/5); % Hinze threshold;
%--------------------------------------------------------------------------
stp_linear = [0.38,0.260,0.39,0.09;0.37,0.31,0.38,0.13;0.40,0.11,0.32,0.07;0.38,0.12,0.31,0.07];  % steepness calculated from dispersion.m
stp_nonlinear = [1.02,0.62,0.80,0.78;0.99,0.63,0.91,0.81;0.70,0.51,0.63,0.69;0.69,0.66,0.53,0.66];
%--------------------------------------------------------------------------
beta_stpln = 0.4.*(stp_linear - 0.08).^(5/2); % breaking strength 
epsilon_bar_ln = (rho/g).*beta_stpln.*(crest_v.^5); % energy diss rate in kg/mm^3* mm/s^2
%--------------------------------------------------------------------------

%% -------------------------------time-------------------------------------
%% GW
for wave_idx = 1:2 
    for phase_idx = 1
        for ROI_idx = 1:8 
            time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = (127301:134800)./frame_rate; 
        end
    end
end
for wave_idx = 1:2 
    for phase_idx = 2
        for ROI_idx = 1:8 
            time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:)  = (125801:133300)./frame_rate;  
        end
    end
end
for wave_idx = 1:2 
    for phase_idx = 3
        for ROI_idx = 1:8 
            time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:)  = (126501:134000)./frame_rate;
        end
    end
end
for wave_idx = 1:2 
    for phase_idx = 4
        for ROI_idx = 1:8 
            time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:)  = (128501:136000)./frame_rate;
        end
    end
end
%% JS
for wave_idx = 3 : 4 
    for phase_idx = 1
        for ROI_idx = 1:8 
            time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = (127301:134800)./frame_rate; 
        end
    end
end  
for wave_idx = 3 : 4 
    for phase_idx = 2
        for ROI_idx = 1:8 
            time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:)  = (126001:133500)./frame_rate; 
        end
    end
end 
for wave_idx = 3 : 4 
    for phase_idx = 3
        for ROI_idx = 1:8             
            time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:)  = (126501:134000)./frame_rate;
        end
    end
end
for wave_idx = 3 : 4 
    for phase_idx = 4
        for ROI_idx = 1:8         
            time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:)  = (128001:135500)./frame_rate;
        end
    end
end

for wave_idx = 1:4 
    for phase_idx = 1:4
        for ROI_idx = 1:8 
            time_stamp.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = time_stamp_aux.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) - (timeorig(wave_idx,phase_idx));
%             ./deltat(wave_idx,phase_idx)
        end
    end
end
%% -----------------------------space--------------------------------------
% x/lambda or bulge length -coordinate array 
for wave_idx = 1:4
for  phase_idx = 1:4
%     ./bulge_length(wave_idx,phase_idx);
% -break_dist(wave_idx,phase_idx)
    x_dir_mtx.(wave{wave_idx}).(phase{phase_idx})(1,1)= 0 ;
    for ROI_idx = 2:8
        x_dir_mtx.(wave{wave_idx}).(phase{phase_idx})(1,ROI_idx) = x_dir_mtx.(wave{wave_idx}).(phase{phase_idx})(1,ROI_idx-1) - 50;
        x_dir_mtx.(wave{wave_idx}).(phase{phase_idx})(1,ROI_idx) = round(x_dir_mtx.(wave{wave_idx}).(phase{phase_idx})(1,ROI_idx),3);
    end
end
end
for wave_idx = 1:4
for phase_idx = 1:4
    x_dir_mtx_left.(wave{wave_idx}).(phase{phase_idx}) = - x_dir_mtx.(wave{wave_idx}).(phase{phase_idx});
%     x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx}) = x_dir_mtx_left.(wave{wave_idx}).(phase{phase_idx})./lambda;
    x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx}) = round([0.035;0.070;0.105;0.140;0.175;0.210;0.245;0.280;0.315]+break_loc{wave_idx,phase_idx},3);
%     x_dir_mtx_leftnd_median.(wave{wave_idx}).(phase{phase_idx}) = [0.018;

end
end
for wave_idx = 1:4
for phase_idx = 1:4
    s = 1;
    for j = 1 : numel(x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx}))-1 
    x_dir_mtx_leftnd_median.(wave{wave_idx}).(phase{phase_idx})(s) = round(median(x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(j:j+1)),3);
    s = s + 1;
    end
end
end

%% ---------------------------Load Data------------------------------------
%% ---------------------------Plot Options---------------------------------
% load and save depending on the option
%% VOLUME AND VOID FRACTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if P == 1
%% --------------- load wave data:-----------------------------------------
    dir_wave = 'C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\freshwater';
    cd(dir_wave)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear dir;
        comments = ['comparison'];
        for wave_idx = 1:4
%% directory for wave repetitions
          dir_reps_aux = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\',...
              data_type,'\',wave_dir_load(wave_idx),'\ROI1-8\wavereps']; 
          dir_reps(wave_idx) = strcat(dir_reps_aux(1),dir_reps_aux(2),dir_reps_aux(3),dir_reps_aux(4),dir_reps_aux(5));
%% directory for instant data
          dir_ins_aux = ['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\data\',...
              data_type,'\',wave_dir_load(wave_idx),'\ROI1-8\']; 
          dir_ins(wave_idx) = strcat(dir_ins_aux(1),dir_ins_aux(2),dir_ins_aux(3),dir_ins_aux(4),dir_ins_aux(5));
        end
%% directory to save figs.
        dir_save =['C:\Users\Konstaninos\OneDrive - University College London\Desktop\Work\third year\laboratory\bubbles\plots\',data_type,['\P',num2str(P)]];
%% load bubble volume (alpha)
for wave_idx = 1:4
    for  phase_idx = 1:4
        alpha_mean_ROI = {['alpha_mean_ROI1'],['alpha_mean_ROI2'],['alpha_mean_ROI3'],['alpha_mean_ROI4'],['alpha_mean_ROI5'],['alpha_mean_ROI6'],['alpha_mean_ROI7'],['alpha_mean_ROI8']};
        cd(dir_reps(wave_idx));
            for ROI_idx = 1:8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                aux_alpha.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = load(['alpha_mean_',ROI{ROI_idx},'_',wave{wave_idx},'_',phase{phase_idx},'__.mat']);
                alpha_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) =  aux_alpha.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).(alpha_mean_ROI{ROI_idx})(1:NUM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
    end
end
% find void fraction (normalised volume)
for wave_idx = 1:4
    for  phase_idx = 1:4   
        for ROI_idx = 1:8
%            void_fraction_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = alpha_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})./vol_fraction;
           void_fraction_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = alpha_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})./vol_bulge.(wave{wave_idx}).(phase{phase_idx});
           void_fraction_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = void_fraction_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx});
        
           alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = alpha_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx});
           % subtract image bias (spots on camera chip recognised as
           % bubbles) ---> systematic error
           alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1:400) = 0;
        end
    end
end

    clear aux_alpha
%% load instant bubble size ditributions (BSD) 
%first load radius bins
    for  phase_idx = 1:4
        cd(dir_ins(1));
        auxBSDX_median = load('BSDX_median.mat');
        BSDX_median.(phase{phase_idx}) = auxBSDX_median.BSDX_median;
    end
    for i = 1:NUM
        BSDX_median_mtx(i,:) =  BSDX_median.(phase{1})(:);
    end
% time for 3D plots
for wave_idx = 1 : 4 
    for phase_idx = 1:4
        for k = 1:numel(BSDX_median.(phase{1}))
            time_mtx_.(wave{wave_idx}).(phase{phase_idx})(:,k) = time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,:);
        end
    end
end
% ins BSD and ins VBSD are loaded in different ways because they are saved differently in spacetime analysis. This should be fixed in spacetime analysis.m
% this is too complicated here (but works) because the BSD files were not saved in spacetime.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this part is way complicated for no other reason that this hasn't been done in spacetime.m
    auxBSD_ROI = {'aux_BSD_ROI1','aux_BSD_ROI2','aux_BSD_ROI3','aux_BSD_ROI4','aux_BSD_ROI5','aux_BSD_ROI6','aux_BSD_ROI7','aux_BSD_ROI8'};
    for wave_idx = 1:4
        cd(dir_ins(wave_idx));
        for run_idx = 1:run_max
            for  phase_idx = 1:4
                for ROI_idx = 1:4
                    auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{1}) = load(['BSD_ins_',ROI{ROI_idx},'_',wave{wave_idx},'_',phase{phase_idx},'_',ROI_old{1},'_',run{run_idx},'.mat']);
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{1}).(run{run_idx}) =  auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{1}).(auxBSD_ROI{ROI_idx});
                end
                for ROI_idx = 2:5
                    auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{2}) = load(['BSD_ins_',ROI{ROI_idx},'_',wave{wave_idx},'_',phase{phase_idx},'_',ROI_old{2},'_',run{run_idx},'.mat']);
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{2}).(run{run_idx}) =  auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{2}).(auxBSD_ROI{ROI_idx});     
                end
                for ROI_idx = 3:6
                    auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{3}) = load(['BSD_ins_',ROI{ROI_idx},'_',wave{wave_idx},'_',phase{phase_idx},'_',ROI_old{3},'_',run{run_idx},'.mat']);
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{3}).(run{run_idx}) =  auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{3}).(auxBSD_ROI{ROI_idx});
                end
                for ROI_idx = 4:7
                    auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{4}) = load(['BSD_ins_',ROI{ROI_idx},'_',wave{wave_idx},'_',phase{phase_idx},'_',ROI_old{4},'_',run{run_idx},'.mat']);
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{4}).(run{run_idx})=  auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{4}).(auxBSD_ROI{ROI_idx});
                end
                for ROI_idx = 5:8
                    auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{5}) = load(['BSD_ins_',ROI{ROI_idx},'_',wave{wave_idx},'_',phase{phase_idx},'_',ROI_old{5},'_',run{run_idx},'.mat']);
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{5}).(run{run_idx}) =  auxBSD.(phase{phase_idx}).(ROI{ROI_idx}).(ROI_old{5}).(auxBSD_ROI{ROI_idx});
                end
            end
        end
   end
   clear aux_BSD
%% load instant volume size distributions (VBSD)
    for wave_idx = 1:4
        for  phase_idx = 1:4
            VBSD_inst_avg_ROI = {['VBSD_inst_avg_ROI1'],['VBSD_inst_avg_ROI2'],['VBSD_inst_avg_ROI3'],['VBSD_inst_avg_ROI4'],...
    ['VBSD_inst_avg_ROI5'],['VBSD_inst_avg_ROI6'],['VBSD_inst_avg_ROI7'],['VBSD_inst_avg_ROI8']};   

            VBSD_inst_mtx =  {['VBSD_inst_ROI1_mtx'],['VBSD_inst_ROI2_mtx'],['VBSD_inst_ROI3_mtx'],['VBSD_inst_ROI4_mtx'],...
            ['VBSD_inst_ROI5_mtx'],['VBSD_inst_ROI6_mtx'],['VBSD_inst_ROI7_mtx'],['VBSD_inst_ROI8_mtx']};
                cd(dir_reps(wave_idx));
            for ROI_idx = 1:8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                auxVBSD_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = load(['VBSD_inst_avg_',ROI{ROI_idx},'_',wave{wave_idx},'_',phase{phase_idx},'__.mat']);
                VBSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) =  auxVBSD_ins_avg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).(VBSD_inst_avg_ROI{ROI_idx})(1:NUM,:);

                auxVBSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = load(['VBSD_inst_',ROI{ROI_idx},'_mtx_',wave{wave_idx},'_',phase{phase_idx},'__.mat']);
                VBSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) =  auxVBSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).(VBSD_inst_mtx{ROI_idx});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end
    clear auxVBSD_ins_avg auxVBSD_ins_avg_mtx
%% find ensemble averages of instant BSD
% 1. put all oldROI.run into same matrices  
% 2. find average for each time step and each radiu bin
    for wave_idx = 1:4
        for  phase_idx = 1:4
%% ROI1
                runs = 1;
                for ROI_old_idx = 1:1
                    for run_idx = 1:run_max
                    BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{1})(runs,:,:) =  ...
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{1}).(ROI_old{ROI_old_idx}).(run{run_idx}); 
                    runs = runs + 1;
                    end
                end
             for ROI_old_idx = 1:1
% find ensemble average 
                for k = 1:numel(BSDX_median.(phase{1}))
                          BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{1})(:,k) = mean(BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{1})(:,:,k));                  
                end
            end
% ROI2
            runs = 1;
            for ROI_old_idx = 1:2
                for run_idx = 1:run_max
                BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{2})(runs,:,:) =  ...
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{2}).(ROI_old{ROI_old_idx}).(run{run_idx});  
                runs = runs + 1;
                end
            end
            for ROI_old_idx = 1:2
% find ensemble average 
                for k = 1:numel(BSDX_median.(phase{1}))
                          BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{2})(:,k) = mean(BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{2})(:,:,k));                  
                end
            end

% ROI3
            runs = 1;
            for ROI_old_idx = 1:3
                for run_idx = 1:run_max
                BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{3})(runs,:,:) =  ...
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{3}).(ROI_old{ROI_old_idx}).(run{run_idx});   
                runs = runs + 1;
                end
            end
            for ROI_old_idx = 1:3
% find ensemble average 
                for k = 1:numel(BSDX_median.(phase{1}))
                          BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{3})(:,k) = mean(BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{3})(:,:,k));                  
                end
            end
% ROI4
            runs = 1;
            for ROI_old_idx = 1:4
                for run_idx = 1:run_max
                BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{4})(runs,:,:) =  ...
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{4}).(ROI_old{ROI_old_idx}).(run{run_idx});  
                runs = runs + 1;
                end
            end
            for ROI_old_idx = 1:4
% find ensemble average 
                for k = 1:numel(BSDX_median.(phase{1}))
                          BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{4})(:,k) = mean(BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{4})(:,:,k));                  
                end
            end
% ROI5
            runs = 1;
            for ROI_old_idx = 2:5
                for run_idx = 1:run_max
                BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{5})(runs,:,:) =  ...
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{5}).(ROI_old{ROI_old_idx}).(run{run_idx}); 
                runs = runs + 1;
                end
            end
            for ROI_old_idx = 2:5
% find ensemble average 
                for k = 1:numel(BSDX_median.(phase{1}))
                          BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{5})(:,k) = mean(BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{5})(:,:,k));                  
                end
            end
% ROI6
            runs = 1;
            for ROI_old_idx = 3:5
                for run_idx = 1:run_max
                BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{6})(runs,:,:) =  ...
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{6}).(ROI_old{ROI_old_idx}).(run{run_idx});   
                runs = runs + 1;
                end
            end
            for ROI_old_idx = 3:5
% find ensemble average 
                for k = 1:numel(BSDX_median.(phase{1}))
                          BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{6})(:,k) = mean(BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{6})(:,:,k));                  
                end
            end            
% ROI7
            runs = 1;
            for ROI_old_idx = 4:5
                for run_idx = 1:run_max
                BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{7})(runs,:,:) =  ...
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{7}).(ROI_old{ROI_old_idx}).(run{run_idx});   
                runs = runs + 1;
                end
            end
            for ROI_old_idx = 4:5            
% find ensemble average 
                for k = 1:numel(BSDX_median.(phase{1}))
                          BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{7})(:,k) = mean(BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{7})(:,:,k));                  
                end
            end             
% ROI8
            runs  = 1;
            for ROI_old_idx = 5:5
                for run_idx = 1:run_max
                BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{8})(runs,:,:) =  ...
                    BSD_ins_noavg.(wave{wave_idx}).(phase{phase_idx}).(ROI{8}).(ROI_old{ROI_old_idx}).(run{run_idx}); 
                runs = runs + 1;
                end
            end
            for ROI_old_idx = 5:5
% find ensemble average 
                for k = 1:numel(BSDX_median.(phase{1}))
                          BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{8})(:,k) = mean(BSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{8})(:,:,k));                  
                end
            end
        end
    end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% subtract initial image bias from BSDs
for wave_idx = 1:4
    for  phase_idx = 1:4
        for ROI_idx = 1:8
%             BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,:) = BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,:) - BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,:);
%             BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,:) = BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,:) - BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(2,:);

            BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1:400,:) = 0;
%             BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(7300:end,:) = 0;
        end
    end
end
for wave_idx = 1:4
    for  phase_idx = 1:4
        for ROI_idx = 1:8
            for i = NUM_in:NUM
                for k = 1 : numel(BSDX_median.(phase{1}))
                    if BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k)< 0 
                        BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) = 0;
                    end
                end
            end
         end
    end
end

%% the same for VBSD
for wave_idx = 1:4
    for  phase_idx = 1:4
 % find ensemble average 
        for ROI_idx = 1:8
            for k = 1:numel(BSDX_median.(phase{1}))
                      VBSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(VBSD_ins_noavg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,:,k));                  
            end
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% FILTERED DATA
 %% bubble size distribution and bubble volume smoothing:

% smooth in time
% smooth
for wave_idx = 1:4
    for phase_idx = 1:4
        for ROI_idx = 1:8
            alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = smoothdata(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:),'movmedian',11);
        end
%             for i = 1 : NUM
%                 alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(:,i) = smoothdata(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(:,i),'movmedian',3);
%             end
    end
end
% noise 
for wave_idx = 1:4
    for phase_idx = 1:4
        for ROI_idx = 1:8
            alpha_ins_noise.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:)  - alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:);
        end
    end
end
% std of noise and threshold
for wave_idx = 1:4
    for phase_idx = 1:4
        for ROI_idx = 1:8
            alpha_ins_thres.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = 3.*std(alpha_ins_noise.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:));
        end
    end
end


% % filter in space
% for wave_idx = 1:4
%     for phase_idx = 1:4
%         for i = 1:NUM
%                  alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(:,i) = ...
%         smooth(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(:,i),'moving',windwoROI)';
% %                  alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(:,i) = round(alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(:,i),2);
%         end
%     end
% end
for wave_idx = 1:4
    for phase_idx  = 1:4
        for ROI_idx = 1:8
                alpha_ins_avg_max_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = max(alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:));
                alpha_tot_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1) = sum(alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:));
                alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1) = mean(alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:));

        end
    end
end
windowSize = 2; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
% smooth BSD
for wave_idx = 1:4
    for phase_idx = 1:4
        for ROI_idx = 1:8
            for i = 1:NUM
%             for k = 1:numel(BSDX_median.(phase{1}))
%                     BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:) = filter(b,a,BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:));
                    BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:) = smoothdata(BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:),'loess',1);
                    VBSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:) = smoothdata(VBSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:),'loess',3);

            end
        end
    end
end
% BSD noise
% for wave_idx = 1:4
%     for phase_idx = 1:4
%         for ROI_idx = 1:8
%             for k = 1:numel(BSDX_median.(phase{1}))
%                 BSD_noise.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) - BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k);
%                 BSD_noise100.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = BSD_noise.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k)./max(BSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k)).*100;
%             end
%         end
%     end
% end
% volume max, mean etc

for wave_idx = 1:4
    for phase_idx  = 1:4
        for i = NUM_in:NUM
            alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,i) = sum(alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(:,i));
        end
    end
end
for wave_idx = 1:4
    for phase_idx  = 1:4
        alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx}) = max(alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,:));
    end
end
% for wave_idx = 1:4
%     for phase_idx  = 1:4
%             alpha_tot_max_smooth_ROI_idx.(wave{wave_idx}).(phase{phase_idx})= find(alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx}) == alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,:));
%     end
% end
% % mean of all ROI
% for wave_idx = 1:4
%     for phase_idx  = 1:4
%             alpha_tot_smooth_ROImean.(wave{wave_idx}).(phase{phase_idx})(1,1) = (alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,alpha_ins_max_smooth_idx.(wave{wave_idx}).(phase{phase_idx})(1,:))+...
%                 alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,alpha_ins_max_smooth_idx.(wave{wave_idx}).(phase{phase_idx})(2,:))+...
%                 alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,alpha_ins_max_smooth_idx.(wave{wave_idx}).(phase{phase_idx})(3,:))+...
%                 alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,alpha_ins_max_smooth_idx.(wave{wave_idx}).(phase{phase_idx})(4,:))+...
%                 alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,alpha_ins_max_smooth_idx.(wave{wave_idx}).(phase{phase_idx})(5,:))+...
%                 alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,alpha_ins_max_smooth_idx.(wave{wave_idx}).(phase{phase_idx})(6,:))+...
%                 alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,alpha_ins_max_smooth_idx.(wave{wave_idx}).(phase{phase_idx})(7,:))+...
%                 alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,alpha_ins_max_smooth_idx.(wave{wave_idx}).(phase{phase_idx})(8,:)))./8;
%     
%     end
% end

% %void fraction max, mean etc
% for wave_idx = 1:4
%     for phase_idx  = 1:4
%         for ROI_idx = 1:8
%                 void_fraction_ins_max_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = max(void_fraction_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:));
%                 void_fraction_tot_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1) = sum(void_fraction_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:));
%                 
%         end
%         void_fraction_max_smooth.(wave{wave_idx}).(phase{phase_idx})(1,1) = mean(void_fraction_ins_max_smooth.(wave{wave_idx}).(phase{phase_idx})(:,:));
% 
%     end
% end
%% split BSD into bubble size ranges
for wave_idx = 1:4
    for phase_idx  = 1:4
        for ROI_idx = 1:8
% r = min - 1.5 mm (1-6)
            for i = NUM_in:NUM
                BSD_min.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = sum(BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,1:6));
                BSD_min.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = smoothdata(BSD_min.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i),'movmedian',21);
            end
        end
    end
end
% r = 1.7-3.9 (7-16)
for wave_idx = 1:4
    for phase_idx  = 1:4
        for ROI_idx = 1:8
            for i = NUM_in:NUM
                BSD_med.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = sum(BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,7:16));
                BSD_med.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = smoothdata(BSD_med.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i),'movmedian',21);
            end
        end
    end
end
% r = 4.1-7.28 (17-30)
for wave_idx = 1:4
    for phase_idx  = 1:4
        for ROI_idx = 1:8
            for i = NUM_in:NUM
                BSD_lrg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = sum(BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,17:end));
                BSD_lrg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = smoothdata(BSD_lrg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i),'movmedian',31);
            end
        end
    end
end
% add all ROI
for wave_idx = 1:4
    for phase_idx  = 1:4
            for i = NUM_in:NUM
                BSD_minR.(wave{wave_idx}).(phase{phase_idx})(i) = 0;
                BSD_medR.(wave{wave_idx}).(phase{phase_idx})(i) = 0;
                BSD_lrgR.(wave{wave_idx}).(phase{phase_idx})(i) = 0;
                for ROI_idx = 1:8
                    BSD_minR.(wave{wave_idx}).(phase{phase_idx})(i) = BSD_minR.(wave{wave_idx}).(phase{phase_idx})(i) + BSD_min.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i);
                    BSD_medR.(wave{wave_idx}).(phase{phase_idx})(i) = BSD_medR.(wave{wave_idx}).(phase{phase_idx})(i) + BSD_med.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i);
                    BSD_lrgR.(wave{wave_idx}).(phase{phase_idx})(i) = BSD_lrgR.(wave{wave_idx}).(phase{phase_idx})(i) + BSD_lrg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i);
                end
            end
    end
end
%--------------------------------------------------------------------------
%% segragate instant BSDs for 0.07 sec - every 1 sec (= wave period) starting from the first time volume that is larger than 10^3 
% some time/space arrays
for wave_idx = 1:4
    for phase_idx  = 1:4
        for ROI_idx = 1:8
            if wave_idx ==3 && phase_idx == 2 
                if ROI_idx == 5 || ROI_idx == 6 || ROI_idx ==7
                    idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = 580;
                else
            	idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = find(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) >= 5*10^2 ...
                & alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) <= 10^3);
                end
            
            elseif wave_idx == 4 && phase_idx == 2 
                if ROI_idx == 6 || ROI_idx == 7 || ROI_idx ==8
                    idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = 800;
                else
                    idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = find(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) >= 5*10^2 ...
                & alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) <= 10^3);
                end 
            elseif wave_idx ==3 && phase_idx == 3
                    idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = find(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) >= 5*10^2 ...
                & alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) <= 10^3);
            elseif wave_idx ==4 && phase_idx == 3
                    idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = find(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) >= 5*10^2 ...
                & alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) <= 10^3);            
            elseif wave_idx ==1 && phase_idx == 2
                    idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = find(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) >= 7.5*10^2 ...
                & alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) <= 2.5*10^3); 
            else
                
            idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = find(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) >= 5*10^2 ...
                & alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) <= 5*10^3);
            end

        end
    end
end
% select BSDs within those time periods found above
        clear  BSD_time_sel_1  BSD_time_sel_2  BSD_time_sel_3  BSD_time_sel_4
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for k = 1:numel(BSDX_median.(phase{phase_idx}))
        BSD_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1):idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+200,k));                    

        BSD_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2200,k));

        BSD_time_sel_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+4000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+4200,k));

        BSD_time_sel_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+6000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+6200,k));
        end
    end
end
end
% do the same selection for void fraction
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        alpha_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx, ...
            idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1):idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+200);                    

        alpha_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx, ...
            idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2200);

        alpha_time_sel_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx, ...
            idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+4000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+4200);

        alpha_time_sel_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx, ...
            idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+6000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+6200);
    end
end
end

% and for VBSD
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for k = 1:numel(BSDX_median.(phase{phase_idx}))
        VBSD_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (VBSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1):idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+200,k));                    

        VBSD_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (VBSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2200,k));

        VBSD_time_sel_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (VBSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+4000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+4200,k));

        VBSD_time_sel_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (VBSD_ins_avg_mtx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+6000:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+6200,k));
        end
    end
end
end

% find mean of those selections
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
                for k = 1:numel(BSDX_median.(phase{phase_idx}))
                    BSD_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(BSD_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
                    BSD_time_sel_mean_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(BSD_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
                    BSD_time_sel_mean_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(BSD_time_sel_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
                    BSD_time_sel_mean_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(BSD_time_sel_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));

                    VBSD_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(VBSD_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
                    VBSD_time_sel_mean_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(VBSD_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
                    VBSD_time_sel_mean_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(VBSD_time_sel_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
                    VBSD_time_sel_mean_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = mean(VBSD_time_sel_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
                
                end
                alpha_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = mean(alpha_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:));
                alpha_time_sel_mean_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = mean(alpha_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:));
                alpha_time_sel_mean_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = mean(alpha_time_sel_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:));
                alpha_time_sel_mean_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = mean(alpha_time_sel_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:));
    end
end
end
%% find bsd etc putside those selected time intervals
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for k = 1:numel(BSDX_median.(phase{phase_idx}))
        BSD_time_sel_1_bt.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+200:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2000,k));                    

        BSD_time_sel_2_bt.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = (BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
            (idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2200:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+4000,k));

        end
    end
end
end
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        alpha_time_sel_1_bt.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx, ...
            idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+200:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2000);                    

        alpha_time_sel_2_bt.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx, ...
            idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+2200:idx_time.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1)+4000);

    end
end
end
%% burst rates = BSD1_sel/BDS2_sel
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for k = 1:numel(BSDX_median.(phase{phase_idx}))
        BSD_burst_rate_1_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k) = BSD_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k)./BSD_time_sel_mean_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k);
        BSD_burst_rate_2_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k) = BSD_time_sel_mean_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k)./BSD_time_sel_mean_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k);
        BSD_burst_rate_3_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k) = BSD_time_sel_mean_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k)./BSD_time_sel_mean_4.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k);

        end
    end
end
end
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for k = 1:numel(BSDX_median.(phase{phase_idx}))
            BSD_burst_rate_1_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,isnan(BSD_burst_rate_1_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k))) = 0;
            BSD_burst_rate_2_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,isnan(BSD_burst_rate_2_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k))) = 0;
            BSD_burst_rate_3_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,isnan(BSD_burst_rate_3_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k))) = 0;

%             if isinf(BSD_burst_rate_1_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k))
%             BSD_burst_rate_1_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = BSD_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k);
%             BSD_burst_rate_2_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = BSD_time_sel_mean_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k);
%             BSD_burst_rate_3_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k) = BSD_time_sel_mean_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k);
            
%             end
        end
    end
end
end
%% peak radius in the volume distribution
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for k = 1:numel(BSDX_median.(phase{phase_idx}))
             VBSD_smooth_mean.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k) = mean(VBSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
        end
    end
end
end
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for i = 1:NUM
         [VBSD_smooth_max.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i),VBSD_smooth_peakr.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)]...
             = max(VBSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:));
        end
    end
end
end
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for i = 1:NUM
            Radiuspeak.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = BSDX_median.(phase{phase_idx})(1,VBSD_smooth_peakr.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i));
        end
    end
end
end
for wave_idx = 1:4
    for phase_idx  = 1:4
        for ROI_idx = 1:8
                Radiuspeak_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:) = smoothdata(Radiuspeak.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:),'movmean',31);
        end
    end
end
for wave_idx = 1:4
    for phase_idx  = 1:4
        ROI_plus = 0;
        for ROI_idx = 1:8
            Radiuspeak_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1:150+ROI_plus) = 0;
            ROI_plus = ROI_plus + 50;
        end
    end
end
for wave_idx = 1:4
    for phase_idx  = 1:4
        for ROI_idx = 1:8
            Radiuspeak_smooth_mxt.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:) = Radiuspeak_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:);
        end
    end
end

%% power law scaling for avg BSDs per ROI
% power law fit and slopes for total mean BSD
        Hz_down= 4; Hz_up = 30;
% first total time averaged BSD
for wave_idx = 1:4
    for phase_idx  = 1:4
             for ROI_idx = 1:8
                 for k = 1 : numel(BSDX_median.(phase{phase_idx}))
                    BSD_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k) = mean(BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:,k));
                 end
                    alpha_avg_total_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1) = mean(alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:));
             end
    end
end
for wave_idx = 1:4
    for phase_idx  = 1:4
        for ROI_idx = 1:8
    [f_LH.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}),gof.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}),...
        output_end.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})] = ...
        fit(BSDX_median.(phase{phase_idx})(Hz_down:Hz_up)',BSD_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(Hz_down:Hz_up)','power1');

% BSD_avg_total_smooth
    Q.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = round(f_LH.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).a(1),1);

    epsilon = 1.^(-1/3);  % this needs to be estimated + investigated: Q(t)(air flow rate)*?(t)(energy dissipation rate).^(-1/3)

    beta_active_theory = - 10./3;

    beta_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = f_LH.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).b;

    N_theory.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = Q.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).* epsilon.*BSDX_median.(phase{phase_idx})(Hz_down:Hz_up).^beta_active_theory;

    N_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = Q.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})* ...
        epsilon.*BSDX_median.(phase{phase_idx})(Hz_down:Hz_up).^beta_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx});
        end
    end
end
% JS waves have no slope break and no big bubbles:
% power law fit and slopes for total mean BSD

% e-folding BSDs (increase-decrease)
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        N_experimental_e.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = exp(1).*N_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx});   
        N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = N_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})./exp(1);

        N_theory_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = N_theory.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})./exp(1);
        N_theory_e.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = exp(1).*N_theory.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx});


    end
end
end
%% find percentage of Bsds within e-folds - ** a bit complique **
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for i = 1:NUM
        for k = 1:numel(BSDX_median.(phase{phase_idx}))- Hz_down +1
            % calculate difference between e-folds (higher-lower = positive number)
            N_experimental_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k) = ...
            (N_experimental_e.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k) - N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k))./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1);
            % calculate difference between data and lower e-fold (data-lower = neg,pos,zero?)
            BSD_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) = BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k+Hz_down-1)./alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,i)- N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1);
            % index indicating if points are inside (1) or outside (0) efolds
            numofbsdin.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) = 0;
       
            if  BSD_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k)<= N_experimental_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k)...
                    && BSD_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) >=0
            % points inside efolds
                numofbsdin.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) = 1;               
            elseif BSD_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) < 0
            % point outside efolds
            % do nothing
            elseif BSD_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) > N_experimental_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k)
            % point outside efolds
            % do nothing
            end
        end
        % check how many points (%) per BSD are inside efolds
        % find sum of ones per BSD
        numofbsdintot.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = sum(numofbsdin.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:)); 
        % find percentage of points inside per BSD - total :
        % Hinze_down:Hinze up = 27
        numofbsdintotper.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = numofbsdintot.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)./(Hz_up - Hz_down+1);
        end
        % index for how many BSD are inside with 50% of 70% points inside
        howmanyBSD.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = 0;
        for i = 1:NUM
            if numofbsdintotper.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) > 0.50 % number of bsd that have 50% or more of their points inside efolds
            % how many BSD out of NUM) have 50% of their points inside efolds
             howmanyBSD.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = howmanyBSD.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) + 1;
            end
        end
        % percentage of BSD with 50% (*or whatever %) points inside
        percentBSDin.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = howmanyBSD.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})./NUM;
   end
end
end
%% do the same for time-selected BSDs 
% first, add the fields of struct selected BSDs
for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        BSD_time_sell_alltogether.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = [BSD_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}); BSD_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})];
        alpha_time_sel_alltogether.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})= [alpha_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:);alpha_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:)];
    end
end
end

for wave_idx = 1:4
for phase_idx  = 1:4
    for ROI_idx = 1:8
        for i = 1:402
        for k = 1:numel(BSDX_median.(phase{phase_idx}))- Hz_down + 1
            % calculate difference between selected data and lower e-fold
            BSD_edifsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) = BSD_time_sell_alltogether.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k+Hz_down-1)./alpha_time_sel_alltogether.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) -...
                N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1);
            % index indicating if points are inside (1) or outside (0) efolds (default)
            numofbsdselin.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) = 0;
       
            if  BSD_edifsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k)<= N_experimental_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k)...
                    && BSD_edifsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) >=0
            % points inside efolds
                numofbsdselin.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) = 1;               
            elseif BSD_edifsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) < 0
            % point outside efolds
            % do nothing
            elseif BSD_edifsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k) > N_experimental_edif.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,k)
            % point outside efolds
            % do nothing
            end
        end
        % check how many points (%) per BSD are inside efolds
        % find sum 
        numofbsdselintot.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = sum(numofbsdselin.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:)); 
        % find percentage of points inside per BSD
        numofbsdselintotper.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) = numofbsdselintot.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)./(Hz_up - Hz_down + 1);
        end
        % index for how many BSD are inside with 70% points inside
        howmanyBSDsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = 0;
        for i = 1:402
            if numofbsdselintotper.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i) > 0.50 % number of bsd that have 50% or more of their points inside efolds
            % how many BSD have 70% of their points inside efolds
             howmanyBSDsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = howmanyBSDsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) + 1;
            end
        end
        % percentage of BSD with 50% (*or whatever %) points inside
        percentBSDselin.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}) = howmanyBSDsel.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})./402;
   end
end
end
%% BSD SUM ALL ROI analysis
%% find insant BSD sum of all ROI
for wave_idx = 1:4
for phase_idx  = 1:4
    for i = NUM_in:NUM
        for k = 1 :numel(BSDX_median.(phase{phase_idx}))
            BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) = 0;
            VBSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) = 0;
            
            for ROI_idx = 1:8
                BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) =   BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) + BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k);
                VBSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) =   VBSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) + VBSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,k);
            
            end
        end
    end
    
end
end
% power law scaling for summed over all ROI BSDs
% power law fit and slopes for total mean BSD
        Hz_down= 4; Hz_up = 30;
% first total time averaged BSD
for wave_idx = 1:4
    for phase_idx  = 1:4
         for k = 1 : numel(BSDX_median.(phase{phase_idx}))
            BSD_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,k) = mean(BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(:,k));
            VBSD_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,k) = mean(VBSDROIsum.(wave{wave_idx}).(phase{phase_idx})(:,k));

         end
    end
end
for wave_idx = 1:4
    for phase_idx  = 1:4
         alpha_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx}) = mean(alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,:));
    end
end
%% fit and slope beta for sums
for wave_idx = 1:4
    for phase_idx  = 1:4
    [f_LH_ROIsum.(wave{wave_idx}).(phase{phase_idx}),gof_ROIsum.(wave{wave_idx}).(phase{phase_idx}),output_end_ROIsum.(wave{wave_idx}).(phase{phase_idx})] = ...
        fit(BSDX_median.(phase{phase_idx})(Hz_down:Hz_up)',BSD_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx})(Hz_down:Hz_up)','power1');

% BSD_avg_total_smooth
    Q_ROIsum.(wave{wave_idx}).(phase{phase_idx}) = round(f_LH_ROIsum.(wave{wave_idx}).(phase{phase_idx}).a(1),1);
    beta_experimental_ROIsum.(wave{wave_idx}).(phase{phase_idx}) = f_LH_ROIsum.(wave{wave_idx}).(phase{phase_idx}).b;
    
    N_experimental_ROIsum.(wave{wave_idx}).(phase{phase_idx}) = Q_ROIsum.(wave{wave_idx}).(phase{phase_idx})* ...
        epsilon.*BSDX_median.(phase{phase_idx})(Hz_down:Hz_up).^beta_experimental_ROIsum.(wave{wave_idx}).(phase{phase_idx});
    end
end
% e-folding BSDs (increase-decrease)
for wave_idx = 1:4
for phase_idx  = 1:4
        N_experimental_e_ROIsum.(wave{wave_idx}).(phase{phase_idx})= exp(1).*N_experimental_ROIsum.(wave{wave_idx}).(phase{phase_idx});   
        N_experimental_e_half_ROIsum.(wave{wave_idx}).(phase{phase_idx}) = N_experimental_ROIsum.(wave{wave_idx}).(phase{phase_idx})./exp(1);

end
end
%% radius of peak vol - BSDROIsum
for wave_idx = 1:4
for phase_idx  = 1:4
        for i = 1:NUM
         [VBSD_smooth_maxROIsum.(wave{wave_idx}).(phase{phase_idx})(i),VBSD_smooth_peakrROIsum.(wave{wave_idx}).(phase{phase_idx})(i)]...
             = max(VBSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,:));
        end
end
end
for wave_idx = 1:4
for phase_idx  = 1:4
        for i = 1:NUM
            RadiuspeakROIsum.(wave{wave_idx}).(phase{phase_idx})(i) = BSDX_median.(phase{phase_idx})(1,VBSD_smooth_peakrROIsum.(wave{wave_idx}).(phase{phase_idx})(i));
        end
end
end
for wave_idx = 1:4
    for phase_idx  = 1:4
        Radiuspeak_smoothROIsum.(wave{wave_idx}).(phase{phase_idx})(:) = smoothdata(RadiuspeakROIsum.(wave{wave_idx}).(phase{phase_idx})(:),'movmean',31);
    end
end
for wave_idx = 1:4
    for phase_idx  = 1:4
        ROI_plus = 0;
            Radiuspeak_smoothROIsum.(wave{wave_idx}).(phase{phase_idx})(1:150+ROI_plus) = 0;
            ROI_plus = ROI_plus + 50;
    end
end
%% do the same for ROI SUM BSDs (percentage of bsd in mean fitted line)
for wave_idx = 1:4
    for phase_idx  = 1:4
        for i = 1:NUM
        for k = 1:Hz_up - Hz_down + 1
            % calculate difference between e-folds (higher-lower = positive number)
            N_experimental_edifROIsum.(wave{wave_idx}).(phase{phase_idx})(1,k) = ...
            (N_experimental_e_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,k) - N_experimental_e_half_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,k))./alpha_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx});
            
            % calculate difference between selected data and lower e-fold
            BSD_edifROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) = BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k+Hz_down-1)./alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(i) -...
                N_experimental_e_half_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,k)./alpha_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx});
            
            % index indicating if points are inside (1) or outside (0) efolds (default)
            numofbsdinROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) = 0;
       
            if  BSD_edifROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k)<= N_experimental_edifROIsum.(wave{wave_idx}).(phase{phase_idx})(1,k)...
                    && BSD_edifROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) >=0
            % points inside efolds
                numofbsdinROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) = 1;               
            elseif BSD_edifROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) < 0
            % point outside efolds
            % do nothing
            elseif BSD_edifROIsum.(wave{wave_idx}).(phase{phase_idx})(i,k) > N_experimental_edifROIsum.(wave{wave_idx}).(phase{phase_idx})(1,k)
            % point outside efolds
            % do nothing
            end
        end
        % check how many points (%) per BSD are inside efolds
        % find sum 
        numofbsdintotROIsum.(wave{wave_idx}).(phase{phase_idx})(i) = sum(numofbsdinROIsum.(wave{wave_idx}).(phase{phase_idx})(i,:)); 
        % find percentage of points inside per BSD
        numofbsdintotperROIsum.(wave{wave_idx}).(phase{phase_idx})(i) = numofbsdintotROIsum.(wave{wave_idx}).(phase{phase_idx})(i)./(Hz_up - Hz_down + 1);
        end
        
        % index for how many BSD are inside with 50% points inside
        howmanyBSDROIsum.(wave{wave_idx}).(phase{phase_idx}) = 0;
        
        for i = 1:NUM
            if numofbsdintotperROIsum.(wave{wave_idx}).(phase{phase_idx})(i) > 0.50 % number of bsd that have 50% or more of their points inside efolds
            % how many BSD have 50% of their points inside efolds
             howmanyBSDROIsum.(wave{wave_idx}).(phase{phase_idx}) = howmanyBSDROIsum.(wave{wave_idx}).(phase{phase_idx}) + 1;
            end
        end
        % percentage of BSD with 50% (*or whatever %) points inside
        percentBSDinROIsum.(wave{wave_idx}).(phase{phase_idx}) = howmanyBSDROIsum.(wave{wave_idx}).(phase{phase_idx})./NUM;
   end
end


%% INTEGRAL PROPERTIES
%% curve fitting to integral bubble-breaker properties
% bubble vol vs bulge vol
for wave_idx = 1:4
for phase_idx  = 1:4
     vol_bulge_mtx(wave_idx,phase_idx) = vol_bulge.(wave{wave_idx}).(phase{phase_idx});    
     alpha_tot_max_smooth_ROI_mtx(wave_idx,phase_idx) = alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx});
end
end
% reshape and fit straight line to data
vol_bulge_re= reshape(vol_bulge_mtx',[1,16]) ; % reshape into 1D array
alpha_tot_max_smooth_ROI_re= reshape(alpha_tot_max_smooth_ROI_mtx',[1,16]) ;
[pVV,SVV] = polyfit(vol_bulge_re,alpha_tot_max_smooth_ROI_re,1); % fit line
% [y_fitVV,deltaVV] = polyval(pVV,vol_bulge_re,SVV); % line eq. 
x_fitmarks_VV = linspace(0,max(vol_bulge_re)+3000, numel(vol_bulge_re));
[y_fitVV,deltaVV] = polyval(pVV,x_fitmarks_VV,SVV); % line eq. 

% bulge vol vs air vol
hsq_vol_re= reshape(hsq_vol',[1,16]) ; % reshape into 1D array
[pVA,SVA] = polyfit(hsq_vol_re,vol_bulge_re,1); % fit line
% [y_fitVV,deltaVV] = polyval(pVV,vol_bulge_re,SVV); % line eq. 
x_fitmarks_VA = linspace(0,max(hsq_vol_re)+15000, numel(hsq_vol_re));
[y_fitVA,deltaVA] = polyval(pVA,x_fitmarks_VA,SVA); % line eq. 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------PLOTS------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if P ==3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bubble vol. 2
%% BUBBLE VOLUME contour
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear c colorbar
%     sizefont1 = 18;
%     left = 0.20;
%     bottom = 0.23; %0.17
%     width = 0.70;
%     height = 0.65; %0.75
%     comments = '';
%     title_plot_P01 = 'void fraction';
%     P01 = figure('Name',[title_plot_P01,'_',comments],'NumberTitle','off');
%     set(P01,'units','normalized','outerposition',[0.1 0.1 0.5 0.8]); %0.65 0.7
%     set(P01,'PaperSize',[70 70]);
%     set(P01,'PaperUnits', 'centimeters');
%     set(P01,'Color',[1 1 1]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P011 = axes('Position',[left bottom width height]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     yyaxis right;
%     for wave_idx = 1
%         for phase_idx = 2
%                 title_save_void_fraction_xt= ['bubvol_',wave{wave_idx},'_',phase{phase_idx}, comments];
% %                 pcolor(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1:end),x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1:end),...
% %                     alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(1:end,1:end));
% %                 contourf(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1:end),x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1:end-1),...
% %                     alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(1:end,1:end),5);
% 
%                     hold on;
%                 contour(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1:end),x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1:end-1),...
%                     alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(1:end,1:end),35,'Linewidth',3);
%                % ./vol_bulge.(wave{wave_idx}).(phase{phase_idx}),10)
% %                 ./vol_fraction 
%         end
%     end
%     ylim([x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1) x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(end-1)+0.0001]);
%     colormap(parula);
%     set(gca,'colorscale','log')
% %     caxis([0 max(max(alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(:,:)))]);
%     caxis([10 10^3]);
%     ax = gca;
%     xt = ax.YTick;
%     yt = ax.YTick;
%     ax.YAxis(1).Color = 'black';
%     ax.YAxis(2).Color = 'black';
%     box on;
% % we need to create two y-axis (why? ask matlab) - this is the first one , placed on the right of the plot
%     set(ax,'ytick',[x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(2),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(3),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(4),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(5),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(6),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(7),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(8)])
%     names = {'ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI7','ROI8'};
%     hold on;
%     for wave_idx = 1
%         for phase_idx = 2
% 
% %             plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1:800),1.46.*time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1:800),'black','linewidth',2);
% %             contour(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1:end),x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1:end-1),...
% %                     alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(1:end,1:end),35);
%         end
%     end 
%     set(gca,'yticklabel',names,'FontSize',sizefont1);
%     set(gca,'TickDir','out');    
%     lighting phong;
%     lightangle(-90,30);
%     camlight('left');
%     shading interp;
%     EdgeColor = 'none';
%     yyaxis left
%     % so all this sh&t is needed to get both yaxes with the labels we want
%     % this is for left axis ; tank position info
%     ylim([x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1) x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(end-1)+0.0001]);
%     set(ax,'ytick',[x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(2),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(3),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(4),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(5),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(6),...
%         x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(7),...
%        x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(8)])
%     xlim([0 time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,end)+0.001]);
%     set(gca, 'xtick',[0 1 2 3 4]);
%     xlabel('$\boldmath{t''*f_c}$','FontSize',sizefont1,'interpreter', 'latex');
%     set(gca,'Fontsize',sizefont1);     
%     ylabel('tank position'); 
%     colormap(parula); 
%     c=colorbar; box on;
%     set(gca,'colorscale','log');
% %     caxis([0 max(max(alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(:,:)))]);
%     caxis([10 10^3]);
% %     c.YTicks = {10,100,1000};
% %     c.YTickLabel = {'10', '10^2', '10^3'};
% 
%     h=ylabel(c,'$\boldmath{V_{bls}(x'',t'') (mm^3)}$','Fontsize',sizefont1,'interpreter','latex');
% % % --------------------------------------------------------------------------
% % % save result
%     export_fig([dir_save,'\',title_save_void_fraction_xt],'-tif','-nocrop');
% % -------------------------------------------------------------------------- 
% %% BSD instant contour plot or 3D plot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cd(dir_save);
%     window = 1;
%     sizefont = 20;
%     sizefont1 = 28;
%     left = 0.22;
%     bottom = 0.17;
%     width = 0.7;
%     height = 0.75;
%     comments = '';
%     title_plot_P08 = 'BSD';
%     P08 = figure('Name',[title_plot_P08,'_',comments],'NumberTitle','off');
%     set(P08,'units','normalized','outerposition',[0.1 0.1 0.65 0.9]);
%     set(P08,'PaperSize',[70 70]);
%     set(P08,'PaperUnits', 'centimeters');
%     set(P08,'Color',[1 1 1]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P018 = axes('Position',[left bottom width height]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     cc = hsv(7500);
%     for wave_idx = 1
%         for phase_idx = 4
% 
%             for ROI_idx = 8
%                 title_save = ['bsd_instant_',wave{wave_idx},'_',phase{phase_idx},'_',ROI{ROI_idx} comments];
%                 for i = 1:NUM
% %                     if void_fraction_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,i)>0
%                         plot3(time_mtx_.(wave{wave_idx}).(phase{phase_idx})(i,:)',BSDX_median_mtx(i,:)',...
%                             (BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:).*vol_factor.*norm_bin_size),...
%                            '-','color',cc(i,:),'Linewidth',0.2);
% 
% %                         plot(BSDX_median_mtx(i,:),(BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:).*vol_factor.*norm_bin_size)./alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,i),...
% %                            '-','color',cc(i,:),'Linewidth',0.1);
% %                        ./alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,i)
% % ./alpha_ins_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,i)
%         %                set(gca,'ColorOrder',cc);
%                        hold on;
% %                     end
%                 end
% %                plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
% %                    (N_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
% %                    ,'color','black','linewidth',3)
% %                hold on
% %                plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
% %                    (N_experimental_e.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
% %                    ,'--','color','black','linewidth',2)
% %                plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
% %                    (N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
% %                    ,'--','color','black','linewidth',2)
%             end
% %         xlim([time_mtx_.(wave{wave_idx}).(phase{phase_idx})(1,1) time_mtx_.(wave{wave_idx}).(phase{phase_idx})(end,1)]);
% 
%         end
%     end
%     view([95,55,45]); %25 55 45   --->>95,55,45<<-----
%     set(gca,'Fontsize',sizefont);
%     ax = gca;
%     xax = ax.YAxis;  
%     yay = ax.YAxis;
%     box on; grid on;
%     set(xax,'TickDirection','out') 
%     set(yay,'TickDirection','out') 
%     set(gca,'ZScale','log','Fontsize',sizefont);
% %     set(gca,'XScale','log','Fontsize',sizefont);
% %     set(gca,'YScale','log','Fontsize',sizefont);
% %     ylabel('$\boldmath{d\bar{N}(r)}/dr$ ($\mathrm{mm^{-1}}$ $\mathrm{m{^{-3}}}$)','Fontsize', sizefont1,'interpreter', 'latex');
% %     ylabel('$\boldmath{(d\bar{N}(r)}/dr)/V_{bls}$','Fontsize', sizefont1,'interpreter', 'latex');
% %     xlabel('$\boldmath{r}$ (mm)','Fontsize', sizefont1,'interpreter', 'latex');
% 
%     zlabel('$\boldmath{d\bar{N}(r)}/dr$ ($\mathrm{mm^{-1}}$ $\mathrm{m{^{-3}}}$)','Fontsize', sizefont1,'interpreter', 'latex');
%     ylabel('$\boldmath{r}$ (mm)','Fontsize', sizefont1,'interpreter', 'latex');
%     xlabel('$\boldmath{{t''*f_c}}$','Fontsize', sizefont1,'interpreter', 'latex');
% %     set(gca,'Ydir','reverse');
%     set(gca,'Xdir','reverse');
%     yticks([0.5 1 2 5 10]);
%     xticks([0 1 2 3 4]);
%     zlim([10^0 10^6]);
%     zticks([10^0 10^1 10^2 10^3 10^4 10^5 10^6]);
% 
% %------------------------------------------------
% %     ylim([10^0 10^6]);
% %     yticks([10^0 10^1 10^2 10^3 10^4 10^5 10^6]);
% %------------------------------------------------
% %     ylim([10^-3 10^4]);
% %     yticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4]);
% % xticks([0.2 0.5 1 5 10]);
% %------------------------------------------------
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % save result
%     export_fig([dir_save,'\',title_save],'-png','-nocrop');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% BSD time evolution
%     title_plot_P04 = [];
%     P04 = figure('Name',[title_plot_P04,comments],'NumberTitle','off');
%     sizefont = 26;
%     sizefont1 = 26;
%     left = 0.25;
%     bottom = 0.17;
%     width = 0.60;
%     height = 0.75;
%     set(P04,'units','normalized','outerposition',[0.1 0. 0.45 0.95]);
%     set(P04,'PaperSize',[70 70]);
%     set(P04,'PaperUnits', 'centimeters');
%     set(P04,'Color',[1 1 1]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P014 = axes('Position',[left bottom width height]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear cc
%     cc = hsv(7500);
%     for wave_idx = 1
%         for phase_idx = 4
%             for ROI_idx = 8
%                 for i = 1: NUM
%                     title_save = ['BSD_normalised_sum_',wave{wave_idx},'_',phase{phase_idx},'_',ROI{ROI_idx}]; 
%                     hold on;
% %                     plot(BSDX_median.(phase{phase_idx}),((BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(i,:)).*...
% %                     vol_factor.*norm_bin_size)./alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,i)...
% %                     ,'color',cc(i,:,:),'linewidth',0.3);
% %                 ./alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,i)
% 
%                     plot(BSDX_median.(phase{phase_idx}),((BSD_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:)).*...
%                     vol_factor.*norm_bin_size)./alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,i)...
%                     ,'color',cc(i,:,:),'linewidth',0.3);
% 
%                 end
% %                 ./alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,i)
% 
%            plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
%            (N_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
%            ,'color','black','linewidth',3)
%             hold on
% 
%             plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
%            (N_experimental_e.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
%            ,'--','color','black','linewidth',2)
% 
%             plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
%            (N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
%            ,'--','color','black','linewidth',2)
% 
% %--------------------------------------------------------------------------
% %            plot(BSDX_median.(phase{phase_idx})(Hz_down:30),...
% %            (N_experimental_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,1:27).*vol_factor.*norm_bin_size)./...
% %            alpha_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx})...
% %            ,'color','black','linewidth',3)
% %             hold on;
% % %             N_experimental_ROIsu?,alpha_avg_ROIsum,N_experimental_e_half_ROIsum
% %             plot(BSDX_median.(phase{phase_idx})(Hz_down:30),...
% %            (N_experimental_e_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,1:27).*vol_factor.*norm_bin_size)./...
% %             alpha_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx})...
% %            ,'--','color','black','linewidth',2)
% % %             hold on;
% %             plot(BSDX_median.(phase{phase_idx})(Hz_down:30),...
% %            (N_experimental_e_half_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,1:27).*vol_factor.*norm_bin_size)./...
% %            alpha_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx})...
% %            ,'--','color','black','linewidth',2)
% 
% %            plot(BSDX_median.(phase{phase_idx}),...
% %            (BSD_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx}).*vol_factor.*norm_bin_size)./alpha_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,i)...
% %            ,'color','black','linewidth',3)
% %             hold on
% %             ./alpha_avg_ROIsum.(wave{wave_idx}).(phase{phase_idx})(1,i)
% %--------------------------------------------------------------------------
%              end
%         end
%     end   
%     set(gca,'Fontsize',sizefont);
%     grid on;
%     box on;
% %     ylabel('$\boldmath{d\bar{N}(r)}/dr$ ($\mathrm{mm^{-1}}$ $\mathrm{m{^{-3}}}$)','Fontsize', sizefont1,'interpreter', 'latex');
%     ylabel('$\boldmath{(d\bar{N}(r)}/dr)/V_{bls}$','Fontsize', sizefont1,'interpreter', 'latex');
%     xlabel('$\boldmath{r}$ (mm)','Fontsize', sizefont1,'interpreter', 'latex');
% 
%     xlim([0.2 10]);
%     ylim([10^-2 10^4]);
%     yticks([10^0 10^1 10^2 10^3 10^4]);
%     xticks([0.2 0.5 1 5 10]);
% %     xlabel('t(sec)','Fontsize',sizefont1); 
%     set(gca,'YScale','log','Fontsize',sizefont);
%     set(gca,'XScale','log','Fontsize',sizefont);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % save result
%     export_fig([dir_save,'\',title_save],'-png','-nocrop');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% BSD instant time selection (dT = wave period intervals)
%     title_plot_P04 = [];
%     P04 = figure('Name',[title_plot_P04,comments],'NumberTitle','off');
%     sizefont = 26;
%     sizefont1 = 26;
%     left = 0.25;
%     bottom = 0.17;
%     width = 0.60;
%     height = 0.75;
%     set(P04,'units','normalized','outerposition',[0.1 0. 0.45 0.95]);
%     set(P04,'PaperSize',[70 70]);
%     set(P04,'PaperUnits', 'centimeters');
%     set(P04,'Color',[1 1 1]);
%     diff_idx = [2];
%     diff_idx_p = [1,3,4];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P014 = axes('Position',[left bottom width height]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear cc
% %     aa = [winter(4),summer(4),gray(4),autumn(4)];
% %     bb = summer(4);
% %     cc = gray(4);
% %     dd = autumn(4);
%        colors.(wave{1})= winter(4);  colors.(wave{2}) = summer(4);  colors.(wave{3}) = gray(4);  colors.(wave{4}) = autumn(4);
%        colors.(wave{wave_idx})(phase_idx,:,:)
%       for wave_idx = 1:4
%         cc= jet(8);
%         jj = parula(8);
%         cc255 = 255.*cc;
%       end
% %     title_save = ['vbsd_ins_mean']; 
%     c_idx = 1;
%     for wave_idx = 1
%               if wave_idx == 2
%                 c_idx = 5;
%               elseif wave_idx == 3
%                 c_idx = 9;
%               elseif wave_idx ==4
%                 c_idx = 13;
%               end
%         for phase_idx = 1
%             for ROI_idx = 1:8
%                 for i = 1: 201
%                 title_save = ['avg_time_sel_',wave{wave_idx},'_',phase{phase_idx}]; 
% %                 plot(BSDX_median.(phase{phase_idx}),(BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(alpha_tot_max_smooth_ROI_idx.(wave{wave_idx}).(phase{phase_idx})-200:alpha_tot_max_smooth_ROI_idx.(wave{wave_idx}).(phase{phase_idx})+200,:).*...
% %                     vol_factor.*norm_bin_size)...
% %                    ,'^','color','blue','linewidth',0.5);
% %                                hold on;
% % 
% %                 plot(BSDX_median.(phase{phase_idx}),(BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(alpha_tot_max_smooth_ROI_idx.(wave{wave_idx}).(phase{phase_idx}),:).*...
% %                     vol_factor.*norm_bin_size)...
% %                    ,'-','color','black','linewidth',3);
% % %                cc(c_idx,:,:
% %                 c_idx = c_idx + 1;
% %                  plot(BSDX_median.(phase{phase_idx})(Hz_down:Hz_up),(130*N_theory.(wave{1}).(phase{1}).(ROI{4}).*vol_factor.*norm_bin_size)...
% %                      ...
% %                     ,'--','color','black','linewidth',3);  
% % 
% %                 plot(BSDX_median.(phase{phase_idx}),(BSDROIsum.(wave{3}).(phase{2})(alpha_tot_max_smooth_ROI_idx.(wave{3}).(phase{2}),:).*...
% %                     vol_factor.*norm_bin_size)...
% %                    ,'-','color','black','linewidth',3);
% %                 plot(BSDX_median.(phase{phase_idx}),(BSDROIsum.(wave{3}).(phase{2})(alpha_tot_max_smooth_ROI_idx.(wave{3}).(phase{2})-100:alpha_tot_max_smooth_ROI_idx.(wave{3}).(phase{2})+100,:).*...
% %                     vol_factor.*norm_bin_size)...
% %                    ,'^','color','green','linewidth',0.3);
% 
% 
% %                ./hsq_vol(wave_idx,phase_idx)                      
% %                ./gamma(wave_idx,phase_idx)
% %                ./alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})
% %                ./vol_bulge.(wave{wave_idx}).(phase{phase_idx})
% 
% %                 plot(BSDX_median.(phase{phase_idx}),(BSD_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:)...
% %                     .*vol_factor.*norm_bin_size)./alpha_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)...
% %                    ,'color',[253 210 51]/255,'linewidth',.5);
% % %                ./alpha_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)
% %                 hold on;
% %                 plot(BSDX_median.(phase{phase_idx}),(BSD_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:)...
% %                     .*vol_factor.*norm_bin_size)/alpha_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)...
% %                    ,'color','cyan','linewidth',0.5);
% 
% %                .*vol_factor.*norm_bin_size)
% %                ./alpha_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)
% %                 plot(BSDX_median.(phase{phase_idx}),(BSD_time_sel_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i,:).*...
% %                     vol_factor.*norm_bin_size)...
% %                    ,'color','red','linewidth',.5);
% %                ./alpha_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)
% 
% %                plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
% %                    (N_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
% %                    ,'color','black','linewidth',3)
% %                hold on
% %                plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
% %                    (N_experimental_e.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
% %                    ,'--','color','black','linewidth',2)
% %                plot(BSDX_median.(phase{phase_idx})(Hz_down:end),...
% %                    (N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)...
% %                    ,'--','color','black','linewidth',2)
% 
%                 plot(BSDX_median.(phase{phase_idx}),(BSD_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,:).*...
%                     vol_factor.*norm_bin_size)./alpha_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
%                    ,'color',cc(ROI_idx,:,:),'linewidth',3);
% % %                ./alpha_time_sel_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)
%                 hold on;
%                 plot(BSDX_median.(phase{phase_idx}),(BSD_time_sel_mean_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,:).*...
%                     vol_factor.*norm_bin_size)./alpha_time_sel_mean_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})...
%                    ,'color',cc(ROI_idx,:,:),'linewidth',3);
% 
% %                 plot(BSDX_median.(phase{phase_idx}),(BSD_time_sel_mean_3.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,:).*...
% %                     vol_factor.*norm_bin_size)...
% %                    ,'color','black','linewidth',1);
% %                ./alpha_time_sel_2.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i)
% 
% 
% % 
% %                  plot(BSDX_median.(phase{2})(Hz_down:Hz_up),(N_experimental.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)...
% %                      ,'color','black','linewidth',4); 
% %                  plot(BSDX_median.(phase{2})(Hz_down:Hz_up),(N_experimental_e.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)...
% %                      ,'--','color','black','linewidth',2);  
% %                  plot(BSDX_median.(phase{2})(Hz_down:Hz_up),(N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)...
% %                     ,'--','color','black','linewidth',2);  
% %                 ./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)
% % ./alpha_avg_total_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1)
%                 end
%             end
%         end
%     end
%     set(gca,'Fontsize',sizefont);
%     grid on;
%     box on;
% %     ylabel('$\mathrm{\bar{N}(r)}$ ($\mathrm{mm^{-1}}$ $\mathrm{m{^{-3}}}$)','Fontsize', sizefont1,'interpreter', 'latex');
%     ylabel('$\boldmath{(d\bar{N}(r)}/dr)/V_{bls}$','Fontsize', sizefont1,'interpreter', 'latex');
% %     ylabel('$\mathrm{\bar{N}_{max}(r)}/\mathrm{V_{bls}}$ ','Fontsize', sizefont1,'interpreter', 'latex');
% %     ylabel('$\boldmath{\mathrm{\bar{N}(r)}/\mathrm{V_{bls}}} (\mathrm{mm^{-1}}$ $\mathrm{m{^{-6}}})$ ','Fontsize', sizefont1,'interpreter', 'latex');
% 
% %     ylabel('$\boldmath{\mathrm{\bar{N}_{max}(r)}/\epsilon_l}$ ','Fontsize', sizefont1,'interpreter', 'latex');
% 
%     xlim([0.2 10]);
%     ylim([10^-2 10^3]);
%     yticks([10^0 10^1 10^2 10^3]);
% 
% %     ylim([10^-3 10^3]);
% %     ylim([0 10^6]);
% %     yticks([10^-1 10^1 10^3 10^6]);
%     xticks([0.2 0.5 1 5 10]);
%     xlabel('$r$(mm)','Fontsize',sizefont1,'interpreter', 'latex'); 
%     set(gca,'YScale','log','Fontsize',sizefont);
%     set(gca,'XScale','log','Fontsize',sizefont);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % save result
%     export_fig([dir_save,'\',title_save],'-png','-nocrop');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% VBSD 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     title_plot_P04 = [];
%     P04 = figure('Name',[title_plot_P04,comments],'NumberTitle','off');
%     sizefont = 26;
%     sizefont1 = 26;
%     left = 0.25;
%     bottom = 0.17;
%     width = 0.60;
%     height = 0.75;
%     set(P04,'units','normalized','outerposition',[0.1 0. 0.45 0.95]);
%     set(P04,'PaperSize',[70 70]);
%     set(P04,'PaperUnits', 'centimeters');
%     set(P04,'Color',[1 1 1]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P014 = axes('Position',[left bottom width height]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear cc
%     cc = jet(8);
%     jj = {'blue','red'};
%     for wave_idx = 3
%         for phase_idx = 2
%             for ROI_idx = 1:8
%                     title_save = ['VBSD_mean',wave{wave_idx},'_',phase{phase_idx},'_']; 
%                     hold on;
%                     plot(BSDX_median.(phase{phase_idx}),VBSD_smooth_mean.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,:)...
%                         ./alpha_ins_avg_total_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx)...
%                     ,'-','color',cc(ROI_idx,:,:),'linewidth',4);
% 
% %                     plot(BSDX_median.(phase{phase_idx}),VBSD_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1,:)...
% %                         ./alpha_time_sel_mean_1.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:)...
% %                     ,'-','color',cc(ROI_idx,:,:),'linewidth',5);
% 
%             end
%         end
%     end   
%     set(gca,'Fontsize',sizefont);
%     grid on;
%     box on;
% %     ylabel('$\boldmath{d\bar{N}(r)}/dr$ ($\mathrm{mm^{-1}}$ $\mathrm{m{^{-3}}}$)','Fontsize', sizefont1,'interpreter', 'latex');
%     ylabel('$\boldmath{V_{bls}(r)/V_{tot}}$','Fontsize', sizefont1,'interpreter', 'latex');
%     xlabel('$\boldmath{r}$ (mm)','Fontsize', sizefont1,'interpreter', 'latex');
% 
% %     ylim([0 0.15]);
%     xticks([0.2 1 3 5 7 10]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % save result
%     export_fig([dir_save,'\',title_save],'-png','-nocrop');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% peak volume radius 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    title_plot_P04 = [];
    P04 = figure('Name',[title_plot_P04,comments],'NumberTitle','off');
    sizefont = 22;
    sizefont1 = 22;
    left = 0.15;
    bottom = 0.32;
    width = 0.750;
    height = 0.6;
    set(P04,'units','normalized','outerposition',[0.1 0. 0.8 0.55]);
    set(P04,'PaperSize',[70 70]);
    set(P04,'PaperUnits', 'centimeters');
    set(P04,'Color',[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P014 = axes('Position',[left bottom width height]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear cc gca ax
%     cc = hsv(NUM);
%     jj = jet(8);
%     for wave_idx = 3
%         for phase_idx = 2
%             for ROI_idx = 1:8
%             title_save = ['radiuspeak_volume_',wave{wave_idx},'_',phase{phase_idx}]; 
% %             for i = 1:NUM
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,i),...
% %                    Radiuspeak_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(i),'o','color',cc(i,:,:));
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,:),...
% %                    Radiuspeak_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:),'-','color',jj(ROI_idx,:,:),'LineWidth',1);
% 
%             hold on;
%             xlim([time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1) time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,end)+ 0.05]);
% 
% %             end
%             end
%         end
%     end
%     set(gca,'Fontsize',sizefont);
%     ylim([0 7.5]);
%     yticks([0 1 2 3 4 5 6 7]);
%     xticks([0 1 2 3 4]);
%     xlabel('$\boldmath{t''*f_c}$','FontSize',sizefont1,'interpreter','latex');
%     ylabel({'radius of peak'; 'volume (mm)'},'Fontsize',sizefont1);
%     box on ; grid on;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for wave_idx = 2
        for phase_idx = 1
                title_save = ['radpeakvol_',wave{wave_idx},'_',phase{phase_idx}, comments];
                
                
%                 pcolor(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,:),x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1:end-1),...
%                     Radiuspeak_smooth_mxt.(wave{wave_idx}).(phase{phase_idx})(:,1:end));
                    hold on;
                contourf(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,:),x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1:end-1),...
                    Radiuspeak_smooth_mxt.(wave{wave_idx}).(phase{phase_idx})(:,1:end),5);
        end
    end
%     ylim([x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1) x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(end-1)+0.0001]);
    colormap(jet);
    caxis([0 7.5]);
    lighting phong;
    lightangle(-90,30); % -90,30
    camlight('left');
    shading interp;
    EdgeColor = 'none';
%     yyaxis left;

    ax = gca;
%     ax.YAxis(1).Color = 'black';
%     ax.YAxis(2).Color = 'black';
    set(ax,'Fontsize',sizefont1);
    set(ax,'TickDir','out'); 
    set(ax,'ytick',[x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(1),...
        x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(2),...
        x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(3),...
        x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(4),...
        x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(5),...
        x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(6),...
        x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(7),...
       x_dir_mtx_leftnd.(wave{wave_idx}).(phase{phase_idx})(8)])
%     names = {'ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI7','ROI8'};
%     set(ax,'yticklabel',names,'FontSize',sizefont1);
   
%     c=colorbar;
%     caxis([0 7.5]);
%     dx = 0.1; dy= 0.005;
%     a =  c.Position; %gets the positon and size of the color bar
%     set(c,'Position',[a(1)+dx a(2)+dy 0.015 0.6]);% To change size
%     h=ylabel(c,{'radius of peak'; 'volume (mm)'},'Fontsize',sizefont1);

%     ylim([0 7.5]);
%     yticks([0 1 2 3 4 5 6 7]);
    xticks([0 1 2 3 4]);
    xlabel('$\boldmath{t''*f_c}$','FontSize',sizefont1,'interpreter','latex');
    ylabel({'tank position'},'Fontsize',sizefont1);
%     box on ; grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % save result
    export_fig([dir_save,'\',title_save],'-png','-nocrop');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
% %% total instant volume timeseries 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     title_plot_P05 = [];
%     P05 = figure('Name',[title_plot_P05,comments],'NumberTitle','off');
%     sizefont = 24;
%     sizefont1 = 24;
%     left = 0.17;
%     bottom = 0.32;
%     width = 0.750;
%     height = 0.6;
%     set(P05,'units','normalized','outerposition',[0.1 0. 0.9 0.6]);
%     set(P05,'PaperSize',[70 70]);
%     set(P05,'PaperUnits', 'centimeters');
%     set(P05,'Color',[1 1 1]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P015 = axes('Position',[left bottom width height]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear cc
%     cc = hsv(NUM);
%     jj = jet(8);
%     for wave_idx = 3
%         for phase_idx = 2
%             for ROI_idx = 1:8
%             title_save = ['volume_tseries_',wave{wave_idx},'_',phase{phase_idx},'_']; 
% %             for i = 1:NUM
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,i),...
% %                     alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,i)./alpha_ins_avg_max_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx),'o','color',cc(i,:,:));
%                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,:),...
%                     alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,:),'-','color',jj(ROI_idx,:,:),'Linewidth',3);
%             hold on;
%             xlim([time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1) time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,end)+ 0.05]);
% %             end
%             end
%         end
%     end
%     set(gca,'Fontsize',sizefont);
% %     ylim([0 1.1]);
%     yticks([1000 2000 3000 4000 5000]);
% %     yticks([0 500 1000 1500]);
% 
%     xticks([0 1 2 3 4]);
%     xlabel('$\boldmath{t''*f_c}$','FontSize',sizefont1,'interpreter','latex');
% %     ylabel('$\boldmath{V_{bls}(t'')}/\boldmath{V_{max}}$','Fontsize',sizefont1,'interpreter','latex');
% 
%     ylabel('$\boldmath{V_{bls}(t'')} \mathrm{(mm^3)}$','Fontsize',sizefont1,'interpreter','latex');
% 
%     box on ; grid on;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % save result
%     export_fig([dir_save,'\',title_save],'-png','-nocrop');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% scatter plot peak radius - volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    title_plot_P05 = [];
    P05 = figure('Name',[title_plot_P05,comments],'NumberTitle','off');
    sizefont = 20;
    sizefont1 = 18;
    left = 0.17;
    bottom = 0.32;
    width = 0.750;
    height = 0.6;
    set(P05,'units','normalized','outerposition',[0.1 0. 0.4 0.8]);
    set(P05,'PaperSize',[70 70]);
    set(P05,'PaperUnits', 'centimeters');
    set(P05,'Color',[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P015 = axes('Position',[left bottom width height]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for wave_idx = 1:4
for phase_idx = 1:4
    for ROI_idx = 1:8
        Radiuspeak_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(Radiuspeak_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(:)==0) = nan;
    
    end
end
end
    clear cc
    cc = hsv(NUM);
    jj = jet(8);
    for wave_idx = 2
        for phase_idx = 1
            for ROI_idx = 1:8
            title_save = ['scatter_radi_vol_',wave{wave_idx},'_',phase{phase_idx},'_']; 
%             for i = 1:NUM
                scatter(alpha_ins_smooth.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1:10:end),...
                    Radiuspeak_smooth.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(1:10:end),'o','MarkerEdgeColor','white','MarkerFaceColor',jj(ROI_idx,:,:));
            hold on;
%             end
            end
        end
    end
    % for wave_idx = 2
    %     for phase_idx = 1
    %         title_save = ['scatter_radi_vol_sum_',wave{wave_idx},'_',phase{phase_idx},'_']; 
    %             scatter(alpha_tot_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,1:1:end),...
    %                 Radiuspeak_smoothROIsum.(wave{wave_idx}).(phase{phase_idx})(1:1:end),'o','MarkerFaceColor','black');
    %         hold on;
    %     end
    % end
    box on;grid on;
    set(gca,'Fontsize',sizefont);
    ylim([0 8]);
    xlim([10^1 10^4]);
    yticks([0 1 2 3 4 5 6 7 8]);
    xticks([10^1 10^2 10^3 10^4 10^5]);
    set(gca,'XScale','log','Fontsize',sizefont);
%      set(gca,'YScale','log','Fontsize',sizefont);
   
    ylabel({'radius of peak'; 'volume (mm)'},'Fontsize',sizefont1);
    xlabel('$\boldmath{V_{bls}} \mathrm{(mm^3)}$','Fontsize',sizefont1,'interpreter','latex');
%     xlabel('$\boldmath{V_{bls}(t'')}/\boldmath{V_{max}}$','Fontsize',sizefont1,'interpreter','latex');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % save result
    export_fig([dir_save,'\',title_save],'-png','-nocrop');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 






%% colorbar with time scale
    clear P07
    title_plot_P07 = [];
    P07 = figure('Name',[title_plot_P07,comments],'NumberTitle','off');
    sizefont = 24;
    sizefont1 = 34;
    left = 0.22;
    bottom = 0.17;
    width = 0.60;
    height = 0.80;
    limmax_B = 10^5;
    limmax_text =9.*10^4;
    set(P07,'units','normalized','outerposition',[0.4 0. 0.45 0.98]);
    set(P07,'PaperSize',[70 70]);
    set(P07,'PaperUnits', 'centimeters');
    set(P07,'Color',[1 1 1]);
    diff_idx = [2];
    diff_idx_p = [1,3,4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P017 = axes('Position',[left bottom width height]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for wave_idx = 1
    for phase_idx = 4
    title_save_colorbar_time = ['colorbar_time_',wave{wave_idx},'_',phase{phase_idx},'_']; 
    ax = axes;
    cc = colormap(hsv(7500));
    ax.Visible = 'off';
    xx = axes('CLim', [time_mtx_.(wave{wave_idx}).(phase{phase_idx})(1,1), time_mtx_.(wave{wave_idx}).(phase{phase_idx})(end,1)]);
    c = colorbar(xx);
    c.YTick = [time_mtx_.(wave{wave_idx}).(phase{phase_idx})(1,1), time_mtx_.(wave{wave_idx}).(phase{phase_idx})(2000,1),...
        time_mtx_.(wave{wave_idx}).(phase{phase_idx})(4000,1),time_mtx_.(wave{wave_idx}).(phase{phase_idx})(6000,1),time_mtx_.(wave{wave_idx}).(phase{phase_idx})(end,1)];
    set(c,'Location','west');
    set(gca,'FontSize',sizefont);
%     c.YTickLabel = {'0.0005', '1', '2', '3','3.75'};
%     h=ylabel(c,'$\mathrm{t*f_c}$','Fontsize',24,'interpreter','latex');  
    xx.Visible = 'off';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%save result
    export_fig([dir_save,'\',title_save_colorbar_time],'-png','-nocrop');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%% extras:
%% bubble bursting 
%     title_plot_P04 = [];
%     P04 = figure('Name',[title_plot_P04,comments],'NumberTitle','off');
%     sizefont = 26;
%     sizefont1 = 26;
%     left = 0.25;
%     bottom = 0.17;
%     width = 0.60;
%     height = 0.75;
%     limmax_B = 10^5;
%     limmax_text =9.*10^4;
%     set(P04,'units','normalized','outerposition',[0.1 0. 0.45 0.95]);
%     set(P04,'PaperSize',[70 70]);
%     set(P04,'PaperUnits', 'centimeters');
%     set(P04,'Color',[1 1 1]);
%     diff_idx = [2];
%     diff_idx_p = [1,3,4];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P014 = axes('Position',[left bottom width height]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 for ROI_idx = 4
%                 title_save = ['bubbleburst_trough']; 
% %                 for k = 1:numel(BSDX_median.(phase{1}))
% %                     plot(BSDX_median.(phase{1}),BSD_burst_rate_1_smooth.(wave{1}).(phase{1}).(ROI{ROI_idx})(1,:),'o','color','blue','MarkerSize',13,'linewidth',2);
%                     hold on;
%                     plot(BSDX_median.(phase{2}),BSD_burst_rate_1_smooth.(wave{1}).(phase{2}).(ROI{ROI_idx})(1,:),'o','color','red','MarkerSize',13,'linewidth',2);   
%                     set(gca,'Fontsize',sizefont);
%                     grid on;
%                     box on;
%                     xlim([0.2 10]);
%                     xticks([0.2 0.5 1 5 10]);
%                     xlabel('r(mm)','Fontsize',sizefont1); 
%                     ylabel('bursting factor','Fontsize', sizefont1);
%                     set(gca,'YScale','log','Fontsize',sizefont);
%                     set(gca,'XScale','log','Fontsize',sizefont);
%                     ylim([10^0 10^3]);
% 
% %                 end
%                 end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % save result
%     export_fig([dir_save,'\',title_save],'-png','-nocrop');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% exp fit and decay time
for wave_idx = 1:4
    for phase_idx = 1:4
% find time index of maximum number of bubbles: start of decay 
        BSD_minRd_idx.(wave{wave_idx}).(phase{phase_idx}) = find(BSD_minR.(wave{wave_idx}).(phase{phase_idx}) == max(BSD_minR.(wave{wave_idx}).(phase{phase_idx})));
        BSD_medRd_idx.(wave{wave_idx}).(phase{phase_idx}) = find(BSD_medR.(wave{wave_idx}).(phase{phase_idx}) == max(BSD_medR.(wave{wave_idx}).(phase{phase_idx})));
        BSD_lrgRd_idx.(wave{wave_idx}).(phase{phase_idx}) = find(BSD_lrgR.(wave{wave_idx}).(phase{phase_idx}) == max(BSD_lrgR.(wave{wave_idx}).(phase{phase_idx})));
    end
end
for wave_idx = 1:4
% shorten BSD to decay times
    for phase_idx = 1:4
        BSD_minRdec.(wave{wave_idx}).(phase{phase_idx}) = BSD_minR.(wave{wave_idx}).(phase{phase_idx})(BSD_minRd_idx.(wave{wave_idx}).(phase{phase_idx}):end);
        BSD_medRdec.(wave{wave_idx}).(phase{phase_idx}) = BSD_medR.(wave{wave_idx}).(phase{phase_idx})(BSD_medRd_idx.(wave{wave_idx}).(phase{phase_idx}):end);
        BSD_lrgRdec.(wave{wave_idx}).(phase{phase_idx}) = BSD_lrgR.(wave{wave_idx}).(phase{phase_idx})(BSD_lrgRd_idx.(wave{wave_idx}).(phase{phase_idx}):end);
    end
end

for wave_idx = 1:4
% shorten time arrays so they start from decay time s
    for phase_idx = 1:4
        time_stampdecmin.(wave{wave_idx}).(phase{phase_idx}) = time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_minRd_idx.(wave{wave_idx}).(phase{phase_idx}):end);
        time_stampdecmed.(wave{wave_idx}).(phase{phase_idx}) = time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_medRd_idx.(wave{wave_idx}).(phase{phase_idx}):end);
        time_stampdeclrg.(wave{wave_idx}).(phase{phase_idx}) = time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_lrgRd_idx.(wave{wave_idx}).(phase{phase_idx}):end);
    end
end
for wave_idx = 1:4
% smooth with interpolation: refine x grid (time)
    for phase_idx = 1:4
            tintarraymin.(wave{wave_idx}).(phase{phase_idx}) = linspace(time_stampdecmin.(wave{wave_idx}).(phase{phase_idx})(1),time_stampdecmin.(wave{wave_idx}).(phase{phase_idx})(end),100);
            tintarraymed.(wave{wave_idx}).(phase{phase_idx}) = linspace(time_stampdecmed.(wave{wave_idx}).(phase{phase_idx})(1),time_stampdecmed.(wave{wave_idx}).(phase{phase_idx})(end),70);
            tintarraylrg.(wave{wave_idx}).(phase{phase_idx}) = linspace(time_stampdeclrg.(wave{wave_idx}).(phase{phase_idx})(1),time_stampdeclrg.(wave{wave_idx}).(phase{phase_idx})(end),75);

    end
end
for wave_idx = 1:4
% smooth again with interpolation
    for phase_idx = 1:4
        BSD_minRintrp.(wave{wave_idx}).(phase{phase_idx}) = interp1(time_stampdecmin.(wave{wave_idx}).(phase{phase_idx}),...
            BSD_minRdec.(wave{wave_idx}).(phase{phase_idx}), tintarraymin.(wave{wave_idx}).(phase{phase_idx}),'spline');
        
        BSD_medRintrp.(wave{wave_idx}).(phase{phase_idx}) = interp1(time_stampdecmed.(wave{wave_idx}).(phase{phase_idx}),...
            BSD_medRdec.(wave{wave_idx}).(phase{phase_idx}), tintarraymed.(wave{wave_idx}).(phase{phase_idx}),'pchip');
        
        BSD_lrgRintrp.(wave{wave_idx}).(phase{phase_idx}) = interp1(time_stampdeclrg.(wave{wave_idx}).(phase{phase_idx}),...
            BSD_lrgRdec.(wave{wave_idx}).(phase{phase_idx}), tintarraylrg.(wave{wave_idx}).(phase{phase_idx}),'pchip');
    end
end
for wave_idx = 1:4
    for phase_idx = 1:4
% find local minimum in timeseries: end of decay time
        BSD_minRdm_lgc.(wave{wave_idx}).(phase{phase_idx}) = islocalmin(BSD_minRintrp.(wave{wave_idx}).(phase{phase_idx}),'MinProminence',0.2);
%         BSD_medRdm_min.(wave{wave_idx}).(phase{phase_idx}) = -findpeaks(-BSD_medRintrp.(wave{wave_idx}).(phase{phase_idx}),'MinPeakProminence',0.5);

        BSD_medRdm_lgc.(wave{wave_idx}).(phase{phase_idx}) = islocalmin(BSD_medRintrp.(wave{wave_idx}).(phase{phase_idx}));
        BSD_lrgRdm_lgc.(wave{wave_idx}).(phase{phase_idx}) = islocalmin(BSD_lrgRintrp.(wave{wave_idx}).(phase{phase_idx}));
    end
end
for wave_idx = 1:4
    for phase_idx = 1:4
        BSD_minRdm_m.(wave{wave_idx}).(phase{phase_idx})  = BSD_minRintrp.(wave{wave_idx}).(phase{phase_idx})(BSD_minRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
        BSD_medRdm_m.(wave{wave_idx}).(phase{phase_idx}) = BSD_medRintrp.(wave{wave_idx}).(phase{phase_idx})(BSD_medRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
        BSD_lrgRdm_m.(wave{wave_idx}).(phase{phase_idx})  = BSD_lrgRintrp.(wave{wave_idx}).(phase{phase_idx})(BSD_lrgRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
        
        tintarrayminm.(wave{wave_idx}).(phase{phase_idx}) = tintarraymin.(wave{wave_idx}).(phase{phase_idx})(BSD_minRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
        tintarraymedm.(wave{wave_idx}).(phase{phase_idx}) = tintarraymed.(wave{wave_idx}).(phase{phase_idx})(BSD_medRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
        tintarraylrgm.(wave{wave_idx}).(phase{phase_idx}) = tintarraylrg.(wave{wave_idx}).(phase{phase_idx})(BSD_lrgRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
    end
end
for wave_idx = 1:4
    for phase_idx = 1:4
        BSD_minRdm_idx.(wave{wave_idx}).(phase{phase_idx}) = find(BSD_minRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
        BSD_medRdm_idx.(wave{wave_idx}).(phase{phase_idx}) = find(BSD_medRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
        BSD_lrgRdm_idx.(wave{wave_idx}).(phase{phase_idx}) = find(BSD_lrgRdm_lgc.(wave{wave_idx}).(phase{phase_idx}));
    end
end
% exp fitting to decay
for wave_idx = 1:4
    for phase_idx = 1:4
        g = fittype('b*exp(-c*x)');
        f0minRd.(wave{wave_idx}).(phase{phase_idx}) = fit(tintarraymin.(wave{wave_idx}).(phase{phase_idx})(1:BSD_minRdm_idx.(wave{wave_idx}).(phase{phase_idx}))',...
            BSD_minRintrp.(wave{wave_idx}).(phase{phase_idx})(1:BSD_minRdm_idx.(wave{wave_idx}).(phase{phase_idx}))','exp1');
         
        BSD_minRd.(wave{wave_idx}).(phase{phase_idx}) = f0minRd.(wave{wave_idx}).(phase{phase_idx}).a.*exp(1).^(f0minRd.(wave{wave_idx}).(phase{phase_idx}).b.*...
             tintarraymin.(wave{wave_idx}).(phase{phase_idx}));
         

        f0medRd.(wave{wave_idx}).(phase{phase_idx}) = fit(tintarraymed.(wave{wave_idx}).(phase{phase_idx})(1:BSD_medRdm_idx.(wave{wave_idx}).(phase{phase_idx}))',...
            BSD_medRintrp.(wave{wave_idx}).(phase{phase_idx})(1:BSD_medRdm_idx.(wave{wave_idx}).(phase{phase_idx}))','exp1');
         
        BSD_medRd.(wave{wave_idx}).(phase{phase_idx}) = f0medRd.(wave{wave_idx}).(phase{phase_idx}).a.*exp(1).^(f0medRd.(wave{wave_idx}).(phase{phase_idx}).b.*...
             tintarraymed.(wave{wave_idx}).(phase{phase_idx}));
         
         
        f0lrgRd.(wave{wave_idx}).(phase{phase_idx}) = fit(tintarraylrg.(wave{wave_idx}).(phase{phase_idx})(1:BSD_lrgRdm_idx.(wave{wave_idx}).(phase{phase_idx}))',...
            BSD_lrgRintrp.(wave{wave_idx}).(phase{phase_idx})(1:BSD_lrgRdm_idx.(wave{wave_idx}).(phase{phase_idx}))','exp1');
         
        BSD_lrgRd.(wave{wave_idx}).(phase{phase_idx}) = f0lrgRd.(wave{wave_idx}).(phase{phase_idx}).a.*exp(1).^(f0lrgRd.(wave{wave_idx}).(phase{phase_idx}).b.*...
             tintarraylrg.(wave{wave_idx}).(phase{phase_idx}));
    end
end
%--------------------------------------------------------------------------
% growth, exp fitting to growth
%--------------------------------------------------------------------------
% for wave_idx = 1:4
%     for phase_idx = 1:4
%         g = fittype('b*exp(c*x)');
%         f0minRg.(wave{wave_idx}).(phase{phase_idx}) = fit(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1:BSD_minRd_idx.(wave{wave_idx}).(phase{phase_idx}))',...
%             BSD_smooth.(wave{wave_idx}).(phase{phase_idx})(1:BSD_minRd_idx.(wave{wave_idx}).(phase{phase_idx}))','exp1');
%          
%         BSD_minRd.(wave{wave_idx}).(phase{phase_idx}) = f0minRd.(wave{wave_idx}).(phase{phase_idx}).a.*exp(1).^(f0minRd.(wave{wave_idx}).(phase{phase_idx}).b.*...
%              tintarraymin.(wave{wave_idx}).(phase{phase_idx}));
%          
% 
% 
%          
%    
%         
%         
%     end
% end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 2D Plot of BSD(r,T):
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     cd(dir_save);
% %     comments_P12 = [''];
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     title_save= ['lifetime'];
% %     title_plot_P12 = [];
% %     clear P11
% %     P12 = figure('Name',[title_plot_P12,comments],'NumberTitle','off');
% %     sizefont = 22;
% %     sizefont1 = 22;
% %     left1 = 0.24;
% %     bottom1 = 0.20;
% %     width = 0.70;
% %     height = 0.70;
% %     cc = {'blue','blue'};
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     axes('Position',[left1 bottom1 width height]);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     set(P12,'units','normalized','outerposition',[0 0.08 0.40 0.9]);
% %     set(P12,'PaperSize',[6 7]);
% %     set(P12,'Color',[1 1 1]); 
% %     for wave_idx = 1
% %         for phase_idx = 1:2
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,:),BSD_minR.(wave{wave_idx}).(phase{phase_idx}),'-','linewidth',2,'color','blue')
% %                 hold on;
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_minRd_idx.(wave{wave_idx}).(phase{phase_idx}):end),...
% %                     BSD_minR.(wave{wave_idx}).(phase{phase_idx})(BSD_minRd_idx.(wave{wave_idx}).(phase{phase_idx}):end),'-','linewidth',5,'color','blue');
% %                 hold on;
% %                 scatter(tintarrayminm.(wave{wave_idx}).(phase{phase_idx})(1), ...
% %                     BSD_minRdm_m.(wave{wave_idx}).(phase{phase_idx})(1),350,'+','Linewidth',4,'MarkerEdgeColor','black');
% %                 hold on;
% %                 plot(tintarraymin.(wave{wave_idx}).(phase{phase_idx}),BSD_minRintrp.(wave{wave_idx}).(phase{phase_idx}),'linewidth',1.5,'color','red');
% %                 hold on;
% 
% %                 plot(tintarraymin.(wave{wave_idx}).(phase{phase_idx}),BSD_minRd.(wave{wave_idx}).(phase{phase_idx})...
% %                     ,'Linewidth',2,'color','black');
% 
%  
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,1:BSD_medRd_idx.(wave{wave_idx}).(phase{phase_idx})),...
% %                     BSD_medR.(wave{wave_idx}).(phase{phase_idx})(1:BSD_medRd_idx.(wave{wave_idx}).(phase{phase_idx})),'-','linewidth',2,'color','blue');
% %                 hold on;
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_medRd_idx.(wave{wave_idx}).(phase{phase_idx}):end),...
% %                     BSD_medR.(wave{wave_idx}).(phase{phase_idx})(BSD_medRd_idx.(wave{wave_idx}).(phase{phase_idx}):end),'-','linewidth',5,'color','blue');
% %                 hold on;
% %                 scatter(tintarraymedm.(wave{wave_idx}).(phase{phase_idx})(1), ...
% %                     BSD_medRdm_m.(wave{wave_idx}).(phase{phase_idx})(1),350,'+','Linewidth',4,'MarkerEdgeColor','black');
% %                 hold on;
% %                 plot(tintarraymed.(wave{wave_idx}).(phase{phase_idx}),BSD_medRintrp.(wave{wave_idx}).(phase{phase_idx}),'linewidth',1.5,'color','red');
% %                 hold on;
% 
% %                 plot(tintarraymed.(wave{wave_idx}).(phase{phase_idx}),BSD_medRd.(wave{wave_idx}).(phase{phase_idx})...
% %                     ,'Linewidth',1,'color','black');
% 
% 
% % % 
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,:),BSD_lrgR.(wave{wave_idx}).(phase{phase_idx}),'-','linewidth',2,'color','blue')
% %                 hold on;
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_lrgRd_idx.(wave{wave_idx}).(phase{phase_idx}):end),...
% %                     BSD_lrgR.(wave{wave_idx}).(phase{phase_idx})(BSD_lrgRd_idx.(wave{wave_idx}).(phase{phase_idx}):end),'-','linewidth',5,'color',cc{phase_idx});
% %                 hold on;
% %                 scatter(tintarraylrgm.(wave{wave_idx}).(phase{phase_idx})(1), ...
% %                     BSD_lrgRdm_m.(wave{wave_idx}).(phase{phase_idx})(1),350,'+','Linewidth',4,'MarkerEdgeColor','black');
% %                 hold on;
% %                 plot(tintarraylrg.(wave{wave_idx}).(phase{phase_idx}),BSD_lrgRintrp.(wave{wave_idx}).(phase{phase_idx}),'linewidth',1,'color','red');
% %                 hold on;
% % 
% %                 plot(tintarraylrg.(wave{wave_idx}).(phase{phase_idx}),BSD_lrgRd.(wave{wave_idx}).(phase{phase_idx})...
% %                     ,'Linewidth',1,'color','black');
% 
% 
% % all radii
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,:),BSD_all.(wave{wave_idx}).(phase{phase_idx}),'-','linewidth',2,'color','black');
% %                 hold on;
% 
% %             for ROI_idx = 5
%                 
% % % BSD timeseries for each ROI:                
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_mind_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),...
% %                     BSD_min.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(BSD_mind_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),'o','linewidth',1,'color','green');
% %                 hold on;
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_mind_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),...
% %                     BSD_minf.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}),'linewidth',2,'color','black');
% %                 
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_medd_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),...
% %                     BSD_med.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(BSD_medd_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),'o','linewidth',1,'color','green');
% %                 hold on;
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_medd_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),...
% %                     BSD_medf.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}),'linewidth',2,'color','black');
% %                 
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_lrgd_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),...
% %                     BSD_lrg.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx})(BSD_lrgd_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),'o','linewidth',1,'color','green');
% %                 hold on;
% %                 plot(time_stamp.(wave{wave_idx}).(phase{phase_idx})(1,BSD_lrgd_idx.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}):end),...
% %                     BSD_lrgf.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}),'linewidth',2,'color','black');
% %             end
%         end
%     end
%     set(gca,'Fontsize',sizefont);  
%     ylim([0 8]);
% %     yticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4]);
%     grid on; box on;
% %     set(gca,'YScale','log','Fontsize',sizefont);
% %     set(gca,'xScale','log','Fontsize',sizefont);
%     xticks([0 1 2 3 4]);
%     xlabel('t*f_c','FontSize',sizefont1);
%     ylabel('$\mathrm{N_{r,large}(t)}$','Fontsize', sizefont1,'interpreter', 'latex');

%     ylabel('$\mathrm{N_{med}}$','Fontsize', sizefont1,'interpreter', 'latex');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       export_fig([dir_save,'\',title_save],'-png','-nocrop');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

%% integral properties bubbles and bulges
        sizefont = 26;
        sizefont1 = 26;
        bottom = 0.24;
        left = 0.26;
        width = 0.7;
        height = 0.65;
        fig_title = ['whatev'];  
        a1 = figure('Name',[fig_title],'NumberTitle','off'); 
        axes('Position',[left bottom width height]); 
        set(a1,'units','normalized','outerposition',[0.1 0.1 0.5 0.85]);
        set(a1,'PaperSize',[8 10]);
        set(a1,'Color',[1 1 1]);
        fig_title = ['integrals_vol_vol']; 
%         colors = {[255,255,204]./255;[161,218,180]./255;[65,182,196]./255;[34,94,168]./255}; 
        colors = lines(4);
        for wave_idx = 1:4
            for phase_idx = 1:4
%                 scatter(stp_linear(wave_idx,phase_idx),alpha_tot_mean_smooth.(wave{wave_idx}).(phase{phase_idx}),62,markers{wave_idx},'MarkerEdgeColor','black','LineWidth',1)
                
%                 scatter(stp_linear(wave_idx,phase_idx),beta(wave_idx,phase_idx),62,markers{wave_idx},'MarkerEdgeColor','black','LineWidth',1)
%                 errorbar(stp_linear(wave_idx,phase_idx),beta(wave_idx,phase_idx),beta_error(wave_idx,phase_idx),'vertical','black');

%                 errorbar(hsq(wave_idx,phase_idx),alpha_tot_mean_smooth.(wave{wave_idx}).(phase{phase_idx}),hsq_error(wave_idx,phase_idx),'horizontal','black');

%                 scatter(epsilon_bar(wave_idx,phase_idx),alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx}),62,markers{wave_idx},'MarkerEdgeColor','black','LineWidth',1)


%                 scatter(hsq_vol(wave_idx,phase_idx),vol_bulge.(wave{wave_idx}).(phase{phase_idx}),72,markers{wave_idx},'MarkerEdgeColor','black','LineWidth',1)
%                 hold on;
%                 errorbar(hsq_vol(wave_idx,phase_idx),vol_bulge.(wave{wave_idx}).(phase{phase_idx}),vol_bulge_error(wave_idx,phase_idx),'vertical','black');
%                 errorbar(hsq_vol(wave_idx,phase_idx),vol_bulge.(wave{wave_idx}).(phase{phase_idx}),hsq_vol_error(wave_idx,phase_idx),'horizontal','black');
%                 plot(x_fitmarks_VA,y_fitVA,'black','linewidth',1);

%                 scatter(bulge_length(wave_idx,phase_idx),alpha_max_smooth.(wave{wave_idx}).(phase{phase_idx}),markers{wave_idx},'MarkerEdgeColor','black','LineWidth',1)

%                 scatter(bulge_length(wave_idx,phase_idx),alpha_max_smooth.(wave{wave_idx}).(phase{phase_idx}),markers{wave_idx},'MarkerEdgeColor','black','LineWidth',1)

%                 scatter(gamma(wave_idx,phase_idx),alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx}),62,markers{wave_idx},'MarkerEdgeColor','black','LineWidth',1)
%                 errorbar(1,1000,gamma_error_mean,'horizontal','black','linewidth',1);

                scatter(vol_bulge.(wave{wave_idx}).(phase{phase_idx}),alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})...
                    ,122,markers{phase_idx},'MarkerEdgeColor',colors(wave_idx,:,:),'LineWidth',2);
                hold on;
                errorbar(vol_bulge.(wave{wave_idx}).(phase{phase_idx}),alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx}),vol_bulge_error(wave_idx,phase_idx),'horizontal','black');
                plot(x_fitmarks_VV,y_fitVV,'black','linewidth',1);
                plot(x_fitmarks_VV,y_fitVV+2*deltaVV,'--','color','black');
                plot(x_fitmarks_VV,y_fitVV-2*deltaVV,'--','color','black');
                
%                 scatter(vol_bulge.(wave{2}).(phase{4}),alpha_tot_max_smooth_ROI.(wave{2}).(phase{4})...
%                     ,72,markers{2},'MarkerEdgeColor','red','LineWidth',1);
%                 
%                 scatter(vol_bulge.(wave{3}).(phase{2}),alpha_tot_max_smooth_ROI.(wave{3}).(phase{2})...
%                     ,72,markers{3},'MarkerEdgeColor','blue','LineWidth',1);

%                 scatter(beta(wave_idx,phase_idx),alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx}),82,markers{wave_idx},'MarkerEdgeColor','black','LineWidth',1)
%                 errorbar(beta(wave_idx,phase_idx),alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,1),beta_error(wave_idx,phase_idx),'horizontal','black','linewidth',1);

%                 errorbar(0.4,5000,beta_error_mean,'horizontal','black','linewidth',1);


%                 errorbar(beta(wave_idx,phase_idx),alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})(1,1),beta_error_mean(wave_idx,phase_idx),'horizontal','black','linewidth',1);

%                 plot(vol_bulge_re,y_fitVV + 2*deltaVV,'m',vol_bulge_re,y_fitVV - 2*deltaVV,'m');
                
                set(gca,'Fontsize',sizefont1);
                set(gca,'Color',[1 1 1]);
                ylabel('$ \boldmath{{V}_{max,tot}}$ ($\mathrm{mm^{3}}$)','Fontsize',sizefont,'interpreter','latex');
%                 xlabel('$ \boldmath{V_{cyl}}$ $\mathrm{(mm^3}$)','Fontsize', sizefont,'interpreter','latex');
%                 xlabel('$ \mathrm{\epsilon_{l}}$','Fontsize', sizefont,'interpreter','latex');
%                 xlabel('$ \mathrm{S}$','Fontsize', sizefont,'interpreter','latex');
                xlabel('$ \boldmath{V_{blg}}$ ($\mathrm{mm^{3}}$)','Fontsize', sizefont,'interpreter','latex');
%                 ylim([0 11500]);
%                 xlim([0 800]);
                ylim([0 2*10^4]);
%                 ylim([0. 0.4]);
%                 xlim([0. 0.5]);
                box on;
                grid on;
            end
        end
%--------------------------------------------------------------------------
        export_fig([dir_save,'\',fig_title],'-png','-nocrop');
%-------------------------------------------------------------------------- 
        


        
        
        
else
end




























% %% Hinze scale calculation
%     title_plot_P04 = [];
%     P04 = figure('Name',[title_plot_P04,comments],'NumberTitle','off');
%     sizefont = 26;
%     sizefont1 = 26;
%     left = 0.25;
%     bottom = 0.17;
%     width = 0.60;
%     height = 0.75;
%     limmax_B = 10^5;
%     limmax_text =9.*10^4;
%     set(P04,'units','normalized','outerposition',[0.1 0. 0.45 0.95]);
%     set(P04,'PaperSize',[70 70]);
%     set(P04,'PaperUnits', 'centimeters');
%     set(P04,'Color',[1 1 1]);
%     diff_idx = [2];
%     diff_idx_p = [1,3,4];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P014 = axes('Position',[left bottom width height]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clear cc
% %     aa = [winter(4),summer(4),gray(4),autumn(4)];
% %     bb = summer(4);
% %     cc = gray(4);
% %     dd = autumn(4);
%        colors.(wave{1})= winter(4);  colors.(wave{2}) = summer(4);  colors.(wave{3}) = gray(4);  colors.(wave{4}) = autumn(4);
%        colors.(wave{wave_idx})(phase_idx,:,:)
%       for wave_idx = 1:4
%         cc= jet(16);
%         cc255 = 255.*cc;
%       end
%     title_save = ['hinze_test']; 
%     c_idx = 1;
%     for wave_idx = 4
%               if wave_idx == 2
%                 c_idx = 5;
%               elseif wave_idx == 3
%                 c_idx = 9;
%               elseif wave_idx ==4
%                 c_idx = 13;
%               end
%         for phase_idx = 1
% %             for ROI_idx = 3
% %                 for i = 1: 201
% %                 title_save = ['time_sel_1_norm_',wave{wave_idx},'_',phase{phase_idx},'_',ROI{ROI_idx}]; 
%                 plot(BSDX_median.(phase{phase_idx}),(BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(alpha_tot_max_smooth_ROI_idx.(wave{wave_idx}).(phase{phase_idx})-200:alpha_tot_max_smooth_ROI_idx.(wave{wave_idx}).(phase{phase_idx})+200,:).*...
%                     vol_factor.*norm_bin_size)...
%                    ,'^','color','blue','linewidth',0.5);
%                                hold on;
% 
% %                 plot(BSDX_median.(phase{phase_idx}),(BSDROIsum.(wave{wave_idx}).(phase{phase_idx})(alpha_tot_max_smooth_ROI_idx.(wave{wave_idx}).(phase{phase_idx}),:).*...
% %                     vol_factor.*norm_bin_size)...
% %                    ,'-','color','black','linewidth',3);
% % %                cc(c_idx,:,:
% %                 c_idx = c_idx + 1;
% 
%                 
% 
% %                ./hsq_vol(wave_idx,phase_idx)                      
% %                ./gamma(wave_idx,phase_idx)
% %                ./alpha_tot_max_smooth_ROI.(wave{wave_idx}).(phase{phase_idx})
% %                ./vol_bulge.(wave{wave_idx}).(phase{phase_idx})
%                               
%           
% % 
% %                  plot(BSDX_median.(phase{2})(Hz_down:Hz_up),(N_experimental_e.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)...
% %                      ./alpha_avg_total_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1),'--','color','red','linewidth',2);  
% %                  plot(BSDX_median.(phase{2})(Hz_down:Hz_up),(N_experimental_e_half.(wave{wave_idx}).(phase{phase_idx}).(ROI{ROI_idx}).*vol_factor.*norm_bin_size)...
% %                     ./alpha_avg_total_mtx.(wave{wave_idx}).(phase{phase_idx})(ROI_idx,1),'--','color','red','linewidth',2);  
% 
% %                 end
% %             end
%         end
%     end
%     set(gca,'Fontsize',sizefont);
%     grid on;
%     box on;
%     ylabel('$\boldmath{\mathrm{\bar{N}(r)}}$ ($\mathrm{mm^{-1}}$ $\mathrm{m{^{-3}}}$)','Fontsize', sizefont1,'interpreter', 'latex');
% 
%     xlim([0.2 10]);
% %     ylim([10^-5 10^1]);
% %     yticks([10^-5 10^-3 10^-1 10^1]);
% 
% %     ylim([10^-4 10^2]);
%     ylim([10^0 10^6]);
% %     yticks([10^-1 10^0 10^1 10^2 10^3 10^4 10^5]);
%     xticks([0.2 0.5 1 5 10]);
%     xlabel('r(mm)','Fontsize',sizefont1); 
%     set(gca,'YScale','log','Fontsize',sizefont);
%     set(gca,'XScale','log','Fontsize',sizefont);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % save result
%     export_fig([dir_save,'\',title_save],'-png','-nocrop');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










