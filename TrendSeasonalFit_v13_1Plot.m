function rec_cg = TrendSeasonalFit_v13_1Plot(N_row,N_col,mini,T_cg,conse,B_detect)
% function rec_cg = TrendSeasonalFit_v12Plot(3845,2918,100,0.99,6,2:6);
% CCDC 13.01 version - Zhe Zhu, EROS, USGS
% It is based on 7 bands fitting for Iterative Seasonal, Linear, and Break Models
% This function works for analyzing and plotting one time series pixel
%
%% Revisions: $ Date: 03/14/2015 $ Copyright: Zhe Zhu
%  Version 13.01  Add more categories and update i_start in the end (03/14/2015)
%% Version 13.0   Fit curve for disturbed peroid (03/13/2015)
%  Version 12.20  Convert BT to from K to C before analysis (03/12/2015)
%  Version 12.19  Fit for permanent snow if is more than 75% (03/12/2015)
%  Version 12.18  No change detection if clear observation less than 25% (03/12/2015)
%  Version 12.17  Use median value for very simple model & change magnitude (02/24/2015)
%  Version 12.16  Finding changes in all water pixels (02/24/2015)
%  Version 12.15  Use the original multitemporal cloud mask (02/15/2015)
%  Version 12.14  Do not need example_img in images folder (02/09/2015)
%  Version 12.13: More infromation in "category" (11/10/2014)
%  This version (12.13) is used for the third round of the LCMS project.
%  Command: TrendSeasonalFit_v12Plot(N_row,N_col,min=0.5,T_cg=0.99,Tmax_cg=1-1e-6,conse=6,B_detect=2:6)
%  Version 12.12: Fit for pixels where Fmask fails (11/09/2014)
%  Version 12.11: Bug fixed in num_fc (11/09/2014)
%  Version 12.10: Better multietmporal cloud detection at the beginning (11/06/2014)
%  Version 12.9:  Detect change only for land pixel (water/snow speical case) (10/31/2014)
%  Version 12.8:  Speed up by reducing time for RMSE and model computing (10/17/2014)
%  Version 12.7:  mini rmse should be larger than 10% of the mean (10/13/2014)
%  Version 12.6:  Fit model again when there are a 33.3% more data (10/08/2014)
%  Version 12.5:  Use subset of bands (2-6) for detecting surface change (10/01/2014)
%  Version 12.4:  Only apply multitemporal cloud masking during model initialization (09/29/2014)
%  Version 12.3:  Use subset of bands (3-5) to balance change in diferent dimensions (09/01/2014)
%  This version (12.3) is used for the second round of the LCMS project.
%  Command: TrendSeasonalFit_v12Plot(N_row,N_col,min=1,T_cg=0.99,n_times=3,conse=5,B_detect=3:6)
%  Version 12.2:  Bug fixed in model intialization (08/14/2014)
%  Version 12.1:  Use subset of bands (3-6) to avoid atmosphere influences (08/04/2014)
%% Version 12.0   Detecting change based on probability (07/19/2014)
%  Version 11.6:  No need to change folder name & faster in speed (by Christ Holden 06/06/2014)
%  Version 11.5:  Improved calculation of temporally adjusted RMSE (04/23/2014)
%  Version 11.4:  Revise "rec_cg.category" to better seperate different fit processes (04/01/2014)
%  This version (11.4) is used for generating synthetic data for ACRE project and
%  detecting change for LCMS project.
%  Command: TrendSeasonalFit_v11Plot(N_row,N_col,min=1,T_cg=2,n_times=3,conse=6,B_detect=1:6)
%  Version 11.3:  Add "rec_cg.magnitude" as indicator of change magnitude (04/01/2014)
%  Version 11.2:  Change very simple fit with mean value for start and end of timeseries (04/01/2014)
%  Version 11.1:  Do not need metadata in the image folder to run CCDC (03/25/2014)
%% Version 11.0:  Use change vector magnitude as threshold for detecting change (03/25/2014)
%  Version 10.13: Use subset of bands (1-6) to avoid atmosphere influences (01/31/2014)
%  Version 10.12: More accurate number of days per year "num_yrs" (01/30/2014)
%  Version 10.11: RMSE updates with time series fit (01/26/2014)
%  Version 10.10: Update temperature extreme in recent studies (01/16/2014)
%  Version 10.9:  Find break in max value in any of the band (01/08/2014)
%  Version 10.8:  Add very simple fit with median value for start and end of timeseries (10/21/2013)
%  This version (10.8) is used for generating synthetic data for the LCMS project.
%  Command: TrendSeasonalFit_v10Plot('stack',N_row,N_col,mini=0.5,T_cg=3,n_times=3,conse=6,B_detect=2:6)
%  Version 10.7:  Better multitemporal cloud detection (10/19/2013)
%  Version 10.6:  Add "Tmax_cg" for last step noise removal (10/18/2013)
%  Version 10.5:  Use subset of bands (2-6) to avoid atmosphere influences (10/18/2013)
%  Version 10.4:  Let dynamic fitting for pixels at the beginning (09/23/2013)
%  Version 10.3:  Able to detect change at the verying beginning (09/06/2013)
%  Version 10.2:  Add mini years "mini_yrs" in model intialization (09/03/2013)
%  Version 10.1:  Reduce time for calcuating "v_dif" (09/02/2013)
%% Version 10.0:  Fit for beginning and end of the time series (08/31/2013)
%  Version 9.9:   Only fit more than 50% of Landat images overlap area (08/28/2013)
%  Version 9.8:   Force model fit for persistent snow pixels (08/27/2013)
%  Version 9.7:   Add "rec_cg.category" as indicator of fitting procudure (08/20/2013)
%                 Add rec_cg.change_prob as indicator of change probability (08/20/2013)
%                 Add rec_cg.num_obs ad indicator of number of observations (08/20/2013)
%  Version 9.6:   Remove mininum rmse "mini" and minimum years "mini_yrs" (08/16/2013)
%  Version 9.5:   Model gets more coefficients with more observations (08/16/2013)
%  Version 9.4:   Bug fixed in calculating temporally adjusted rmse (08/01/2013)
%  Version 9.3:   Fit curve again after one year (03/28/2013)
%  This version (9.3) is used for mapping land cover for the IDS project.
%  Command: TrendSeasonalFit_v9Plot('stack',N_row,N_col,T_cg=2,n_times=3,conse=4)
%  Version 9.2:   Use "mini = T_const/T_cg" for small rmse cases (03/26/2013)
%  Version 9.1:   Remove out of range pixels before time series analysis (02/09/2013)
%% Version 9.0:   Using 8 coefficients and lasso fit (02/01/2013)
%  Version 8.4:   Use "max v_slope" instead of "average v_slope" (01/16/2013)
%  Version 8.3:   Start initialization when "time_span" > 1 year (01/16/2013)
%  Version 8.2:   Bug fixed in not fitting models at the begining (01/16/2013)
%  Version 8.1:   Bug fixed in counting "i" and "i_span"(01/13/2013)
%% Version 8.0:   Temporally changing RMSE (01/09/2013)
%% Version 7.3:   Continuous Change Detection and Classification (CCDC) (07/11/2012)
%  This version (7.3) is explained by Zhu, Z. & Woodcock, C.E., Continuous Change
%  Detection and Classification (CCDC) of land cover using all available
%  Landsat data, Remote Sensing of Environment (2014).
%  Command: TrendSeasonalFit_v7Plot('stack',N_row,N_col,T_cg=3,n_times=3,conse=3)
%% Version 1.0:   Continous Monitoring of Forest Disturbance Algorithm (CMFDA) (07/13/2010)
%  This version (1.0) is explained by Zhu, Z., Woodcock, C.E., Olofsson, P.,
%  Continuous monitoring of forest disturbance using all available Landsat
%  data, Remote Sensing of Environment (2012).
%
% Inputs:
% N_row: number of rows
% N_col: number of cols
%
% Outputs:
% rec_cg RECord information about all curves between ChanGes
% rec_cg(i).t_start record the start of the ith curve fitting (julian_date)
% rec_cg(i).t_end record the end of the ith curve fitting (julian_date)
% rec_cg(i).t_break record the first observed break time (julian_date)
% rec_cg(i).coefs record the coefficients of the ith curve
% rec_cg(i).pos record the position of the ith pixel (pixel id)
% rec_cg(i).magnitude record the change vector of all spectral bands
%
% rec_cg(i).category record what fitting procudure and model is used
% cateogry category 5x: persistent snow    4x: Fmask fails
% cateogry category 3x: distubed fit       2x: end fit
% category category 1x: start fit           x: normal procedure
% cateogry category x1: mean value         x4: simple model 
% category category x6: advanced model     x8: full model 

% Basic Tools
addpath('/usr3/graduate/zhuzhe/Algorithms/CCDC/Tools');
% % % Robustfit Tools (no need anymore because it is within this folder)
% addpath('/usr3/graduate/zhuzhe/Algorithms/CCDC/Tools/fast_robust_fit');
% Tools of fast lasso fit
addpath('/usr3/graduate/zhuzhe/Algorithms/CCDC/Tools/glmnet_fit');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  defining variables & constants
%% Constants
% maximum number of coefficient required
% 2 for tri-modal; 2 for bi-modal; 2 for seasonality; 2 for linear;
min_num_c = 4;
mid_num_c = 6;
max_num_c = 8;
% max number of coefficients for the model
num_c = 8;
% number of clear observation / number of coefficients
n_times = 3;
% intialize NUM of Functional Curves
num_fc = 0;
% number of days per year
num_yrs = 365.25;
% number of bytes: int16
num_byte = 2;
% total bands (1-5,7,6,Fmask)
nbands = 8;
% Band for multitemporal cloud/snow detection (green)
num_B1 = 2;
% Band for multitemporal shadow/snow shadow detection (SWIR)
num_B2 = 5;
% Threshold for cloud, shadow, and snow detection.
T_const = 400;
% minimum year for model intialization
mini_yrs = 1;
% number of bands for detection
num_detect = length(B_detect);
% percetn of overall ref for mini rmse
p_min = 0.;
% % no change detection for permanent water pixels
% t_ws = 0.95;
% no change detection for permanent snow pixels
t_sn = 0.75;
% Fmask fails threshold
t_clr = 0.25;
% % robust fit type
% fit_type = 'rlowess';
% % change probability % error area less than 10 pixels per scene
% cg_prob = 2.e-11; % e-8=>1000 e-9=>100 e-10=>10 e-11=>1
% conse is able to calculated
% conse = ceil(log(cg_prob)/log(1-T_cg));
% fprintf('At least %d consecutive observations are needed!\n',conse);
Tmax_cg = chi2inv(1-1e-6,num_detect);
% change T_cg from p & v to cdf
T_cg = chi2inv(T_cg,num_detect); % T_cg = 0.99

% get num of total folders start with "L"
imf=dir('L*'); % folder names
% filter for Landsat folders
imf = regexpi({imf.name}, 'L(T5|T4|E7|C8|ND)(\w*)', 'match');
imf = [imf{:}];
imf = vertcat(imf{:});
% sort according to yeardoy
yeardoy = str2num(imf(:, 10:16));
[~, sort_order] = sort(yeardoy);
imf = imf(sort_order, :);
% number of folders start with "L"
num_t=size(imf,1);
% name of the first stacked image
filename = dir([imf(1,:),'/','L*MTLstack']); %%% Revision v12.16
% read in dimension and zone number of the data
dim = envihdrread([imf(1,:),'/',filename.name]);
% initialize the struct data by RECording of ChanGe (rec_cg)
rec_cg=struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'rmse',[],...
    'pos',[],'change_prob',[],'num_obs',[],'category',[],'magnitude',[]);

pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get ready for Xs & Ys
%% Read in Xs & Ys
% transforming to serial date number (0000 year)
sdate=zeros(num_t,1); % Xs
line_t=zeros(num_t,nbands); %Ys

for i=1:num_t
    im_dir = dir(imf(i, :));
    im = '';
    for f = 1:size(im_dir, 1)
        % use regular expression to match:
        %   'L(\w*)'    Any word begining with L that has any following chars
        %   stk_n       includes stack name somewhere after L
        %   '$'           ends with the stack name (e.g., no .hdr, .aux.xml)
        if regexp(im_dir(f).name, ['L(\w*)', 'stack', '$']) == 1
            im = [imf(i, :), '/', im_dir(f).name];
            break
        end
    end
    
    if i == 1
        dim = envihdrread(im);
        % get row, col
        dim=fliplr(dim);
    end
    
    % Find date for folder imf(i)
    yr = str2num(imf(i, 10:13));
    doy = str2num(imf(i, 14:16));
    sdate(i) = datenum(yr, 1, 0) + doy;
    dummy_name=im;
    fid_t=fopen(dummy_name,'r'); % get file ids
    fseek(fid_t,num_byte*((N_row-1)*dim(2)+N_col-1)*nbands,'bof');
    line_t(i,:)=fread(fid_t,nbands,'int16=>double','ieee-le'); % get Ys
end
fclose('all'); % close all files
% profile on;
% tic
% for iii = 1:N
% mask data
line_m=line_t(:,nbands);

% Only run CCDC for places where more than 50% of images has data
idexist = line_m < 255;
overlap_pct = sum(idexist)/num_t;
if overlap_pct < 0.5
    fprintf('Less than 50%% overlap area (%.2f%%)!\n',100*overlap_pct);
    return;
else
    fprintf('More than 50%% overlap area (%.2f%%)!\n',100*overlap_pct);
end

% convert Kelvin to Celsius
line_t(:,7) = line_t(:,7)*10 - 27315;

% pixel value ranges should follow physical rules
idrange=line_t(:,1)>0&line_t(:,1)<10000&...
    line_t(:,2)>0&line_t(:,2)<10000&...
    line_t(:,3)>0&line_t(:,3)<10000&...
    line_t(:,4)>0&line_t(:,4)<10000&...
    line_t(:,5)>0&line_t(:,5)<10000&...
    line_t(:,6)>0&line_t(:,6)<10000&...
    line_t(:,7)>-9320&line_t(:,7)<7070;

% line_t(:,6) = 10000*(line_t(:,2)-line_t(:,5))./(line_t(:,5)+line_t(:,2));
% line_t(:,4) = 10000*(line_t(:,4)-line_t(:,5))./(line_t(:,4)+line_t(:,5));

% # of clear observatons
idclr = line_m < 2;
% # of all available observations
idall = line_m < 255;
% clear observation percentage
clr_pct = sum(idclr)/sum(idall);
% snow pixels
idsn = line_m == 3;
% percent of snow observations
sn_pct = sum(idsn)/sum(idclr);

% not enough clear observations for change detection
if clr_pct < t_clr    
    % permanent snow pixels
    if sn_pct > t_sn        
        % number of snow pixel within range
        n_sn = sum(idsn);
        
        if n_sn < n_times*min_num_c % not enough snow pixels
            fprintf('Not enough "good" snow observations!\n');
            return;
        else
            % start model fit for snow persistent pixels
            fprintf('Fit permanent snow observations (%.0f%%)!\n',100*sn_pct);
            
            % snow observations are "good" now
            idgood = idsn;
            % Xs & Ys for computation
            clrx = sdate(idgood);
            clry = line_t(idgood,1:nbands-1);
            
            % copy Surface reflectance & Brightness Temperature for plot
            Xs = clrx;
            Ys = clry;
            
            % the first observation for TSFit
            i_start = 1;
            % identified and move on for the next curve
            num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
            
            % defining computed variables
            fit_cft = zeros(max_num_c,nbands-1);
            % rmse for each band
            rmse = zeros(nbands-1,1);
            
            for i_B = 1:nbands-1
                if i_B < nbands-1 % treat saturated and unsaturated pixels differently
                    idgood = clry(:,i_B) < 10000 & clry(:,i_B) > 0; % saturate if ref > 1
                    i_span = sum(idgood);
                    if i_span < min_num_c*n_times % fill value for frequently saturated snow pixels
                        fit_cft(1,i_B) = 10000; % give constant value
                    else % fit for enough unsaturat snow pixels
                        [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx(idgood),clry(idgood,i_B),min_num_c);
                    end
                else % fit for temperature band
                    idgood = clry(:,i_B)>-9320&clry(:,i_B)<7070;
                    [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx(idgood),clry(idgood,i_B),min_num_c);
                end
            end
            
            % updating information at each iteration
            % record time of curve start
            rec_cg(num_fc).t_start = clrx(i_start);
            % record time of curve end
            rec_cg(num_fc).t_end = clrx(end);
            % record break time
            rec_cg(num_fc).t_break = 0; % no break at the moment
            % record postion of the pixel
            rec_cg(num_fc).pos = [N_row,N_col];
            % record fitted coefficients
            rec_cg(num_fc).coefs = fit_cft;
            % record rmse of the pixel
            rec_cg(num_fc).rmse = rmse;
            % record change probability
            rec_cg(num_fc).change_prob = 0;
            % record number of observations
            rec_cg(num_fc).num_obs = n_sn;
            % record fit category
            rec_cg(num_fc).category = 50 + min_num_c; % snow pixel
            % record change magnitude
            rec_cg(num_fc).magnitude = zeros(1,nbands-1);
        end
    else % no change detection for clear observations
        n_clr = sum(idclr);
        if n_clr < n_times*min_num_c % not enough snow pixels
            fprintf('Not enough "good" clear observations!\n');
            return;
        else            
            % start model fit for snow persistent pixels
            fprintf('Fmask failed (%.0f%%)!\n',100*clr_pct);
            % clear and within physical range pixels
            idgood = idclr & idrange;
            
            % Xs & Ys for computation
            clrx=sdate(idgood);
            clry=line_t(idgood,1:nbands-1);
            
            % copy Surface reflectance & Brightness Temperature for plot
            Xs=clrx;
            Ys=clry;
            
            % the first observation for TSFit
            i_start = 1;
            % identified and move on for the next curve
            num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
            
            % defining computed variables
            fit_cft = zeros(max_num_c,nbands-1);
            % rmse for each band
            
            rmse = zeros(nbands-1,1);
            for i_B = 1:nbands-1
                [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx,clry(:,i_B),min_num_c);
            end
            
            % record time of curve start
            rec_cg(num_fc).t_start = clrx(i_start);
            % record time of curve end
            rec_cg(num_fc).t_end = clrx(end);
            % record break time
            rec_cg(num_fc).t_break = 0;
            % record postion of the pixel
            rec_cg(num_fc).pos = [N_row,N_col];
            % record fitted coefficients
            rec_cg(num_fc).coefs = fit_cft;
            % record rmse of the pixel
            rec_cg(num_fc).rmse = rmse;
            % record change probability
            rec_cg(num_fc).change_prob = 0;
            % record number of observations
            rec_cg(num_fc).num_obs = length(clrx);
            % record fit category
            rec_cg(num_fc).category = 40 + min_num_c; % simple model fit at the end
            % record change magnitude
            rec_cg(num_fc).magnitude = zeros(1,nbands-1);
        end
    end
    % clear land or water pixels
else
    fprintf('Seasonal Snow (Snow < %.0f%%)!\n',100*sn_pct);
    fprintf('Fmask works (Clear observation > %.0f%%)!\n',clr_pct*100);
    
    % clear and within physical range pixels
    idgood = idclr & idrange;
    
    % Xs & Ys for computation
    clrx = sdate(idgood);
    clry = line_t(idgood,1:nbands-1);
    
    % caculate median variogram
    var_clry = clry(2:end,:)-clry(1:end-1,:);
    var_med = median(abs(var_clry),1);
    
    
    % copy Surface reflectance & Brightness Temperature for plot
    Xs = clrx;
    Ys = clry;
    
    % start with mininum requirement of clear obs
    i = n_times*min_num_c;
    
    % initializing variables
    % the first observation for TSFit
    i_start = 1;
    % record the start of the model initialization (0=>initial;1=>done)
    BL_train = 0;
    % identified and move on for the next curve
    num_fc = num_fc+1; % NUM of Fitted Curves (num_fc)
    % record the num_fc at the beginning of each pixel
    rec_fc = num_fc;
    
    % while loop - process till the last clear observation - conse
    while i<= length(clrx)-conse
        % span of "i"
        i_span = i-i_start+1;
        % span of time (num of years)
        time_span = (clrx(i)-clrx(i_start))/num_yrs;        
        
        % basic requrirements: 1) enough observations; 2) enough time
        if i_span >= n_times*min_num_c && time_span >= mini_yrs
            
            % initializing model
            if BL_train == 0
                % Tmask: noise removal (good => 0 & noise => 1)
                blIDs = autoMask(clrx(i_start:i+conse),clry(i_start:i+conse,[num_B1,num_B2]),...
                    (clrx(i+conse)-clrx(i_start))/num_yrs,T_const);   
                
                % IDs to be removed
                IDs = i_start:i+conse;
                rmIDs = IDs(blIDs(1:end-conse) == 1);  
                                
                % update i_span after noise removal
                i_span = sum(~blIDs(1:end-conse));
                                
                % check if there is enough observation
                if i_span < n_times*min_num_c
                    % move forward to the i+1th clear observation
                    i = i+1;
                    % not enough clear observations
                    continue;
                    % check if there is enough time
                else
                    % copy x & y
                    cpx = clrx;
                    cpy = clry;
                    
                    % remove noise pixels between i_start & i
                    cpx(rmIDs) = [];
                    cpy(rmIDs,:) = [];
                    
                    % record i before noise removal
                    % This is very important as if model is not initialized
                    % the multitemporal masking shall be done again instead
                    % of removing outliers in every masking
                    i_rec = i;
                    
                    % update i afer noise removal (i_start stays the same)
                    i = i_start+i_span-1;
                    % update span of time (num of years)
                    time_span = (cpx(i)-cpx(i_start))/num_yrs;
                    
                    % check if there is enough time
                    if time_span < mini_yrs
                        % keep the original i
                        i = i_rec;
                        % move forward to the i+1th clear observation
                        i = i+1;
                        % not enough time
                        continue;
                        % Step 2: model fitting
                    else
                        % remove noise
                        clrx = cpx;
                        clry = cpy;
                        
                        % Step 2: model fitting
                        % initialize model testing variables
                        % defining computed variables
                        fit_cft = zeros(max_num_c,nbands-1);
                        % rmse for each band
                        rmse = zeros(nbands-1,1);
                        % value of differnce
                        v_dif = zeros(nbands-1,1);
                        % record the diference in all bands
                        rec_v_dif = zeros(i-i_start+1,nbands-1);
                        
                        for i_B = 1:nbands-1
                            % initial model fit
                            [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = autoTSFit(clrx(i_start:i),clry(i_start:i,i_B),min_num_c);
                        end
                        
                        % adjusted mini RMSE based on p_min % of the mean ref
                        adj_rmse = zeros(nbands-1,1);
                        for i_B = B_detect
                            adj_rmse(i_B) = p_min*(fit_cft(1,i_B)+fit_cft(2,i_B)*(clrx(i_start)+clrx(i))/2);
                        end
                        % adjust agaist a constant
                        adj_rmse(adj_rmse < mini) = mini;
                        
                        adj_rmse = var_med;
                        adj_rmse(adj_rmse < mini) = mini;
                        % normalized to z-score
                        for i_B = B_detect
                            % minimum rmse
                            mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                            
                            % compare the first clear obs
                            v_start = rec_v_dif(1,i_B)/mini_rmse;
                            % compare the last clear observation
                            v_end = rec_v_dif(end,i_B)/mini_rmse;
                            % anormalized slope values
                            v_slope = fit_cft(2,i_B)*(clrx(i)-clrx(i_start))/mini_rmse;
                            
                            % differece in model intialization
                            v_dif(i_B) = abs(v_slope) + abs(v_start) + abs(v_end);
                            % v_dif(i_B) = max(abs(v_slope),max(abs(v_start),abs(v_end)));
                            % v_dif(i_B) = abs(v_slope + v_end - v_start);
                        end
                        v_dif = norm(v_dif(B_detect))^2;

                        % find stable start for each curve
                        if v_dif > T_cg
                            % start from next clear obs
                            i_start = i_start + 1;
                            % move forward to the i+1th clear observation
                            i = i + 1;
                            % keep all data and move to the next obs
                            continue;
                        else
                            % model ready!
                            BL_train = 1;
                            % count difference of i for each iteration
                            i_count = 0;
                            
                            % find the previous break point
                            if num_fc == rec_fc
                                % first curve
                                i_break = 1;
                            else
                                % after the first curve
                                i_break = find(clrx > rec_cg(num_fc-1).t_end);
                                % i_break = find(clrx >= rec_cg(num_fc-1).t_start);
                                i_break = i_break(1);
                            end
                            
                            if i_start > i_break
                                % model fit at the beginning of the time series
                                for i_ini = i_start-1:-1:i_break
                                    if i_start - i_break < conse
                                        ini_conse = i_start - i_break;
                                    else
                                        ini_conse = conse;
                                    end
                                    % value of difference for conse obs
                                    v_dif = zeros(ini_conse,nbands-1);
                                    % record the magnitude of change
                                    v_dif_mag = v_dif;
                                    % chagne vector magnitude
                                    vec_mag = zeros(ini_conse,1);
                                    
                                    for i_conse = 1:ini_conse
                                        for i_B = 1:nbands-1
                                            % absolute difference
                                            v_dif_mag(i_conse,i_B) = clry(i_ini-i_conse+1,i_B)-autoTSPred(clrx(i_ini-i_conse+1),fit_cft(:,i_B));
                                            % normalized to z-scores
                                            if sum(i_B == B_detect)
                                                % minimum rmse
                                                mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                                
                                                % z-scores
                                                v_dif(i_conse,i_B) = (v_dif_mag(i_conse,i_B))/mini_rmse;
                                            end
                                        end
                                        vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                                    end
                                    
                                    % change detection
                                    if min(vec_mag) > T_cg % change detected
                                        break;
                                    elseif vec_mag(1) > Tmax_cg % false change
                                        % remove noise
                                        clrx(i_ini,:) = [];
                                        clry(i_ini,:) = [];
                                        i=i-1; % stay & check again after noise removal
                                    end
                                    
                                    % update new_i_start if i_ini is not a confirmed break
                                    i_start = i_ini;
                                    % update curves
                                    for i_B = 1:nbands-1
                                        [fit_cft(:,i_B),rmse(i_B)] = autoTSFit(clrx(i_start:i),clry(i_start:i,i_B),min_num_c);
                                    end
                                end
                            end                           
                            
                            % enough to fit simple model and confirm a break
                            if i_start - i_break >= conse
                                % defining computed variables
                                fit_cft = zeros(max_num_c,nbands-1);
                                % rmse for each band
                                rmse = zeros(nbands-1,1);
                                for i_B=1:nbands-1
                                    [fit_cft(:,i_B),rmse(i_B)] = autoTSFit(clrx(i_break:i_start-1),clry(i_break:i_start-1,i_B),min_num_c);
                                end
                                
                                % record time of curve start
                                rec_cg(num_fc).t_start = clrx(i_break);
                                % record time of curve end
                                rec_cg(num_fc).t_end = clrx(i_start-1);
                                % record postion of the pixel
                                rec_cg(num_fc).pos = [N_row,N_col];
                                % record fitted coefficients
                                rec_cg(num_fc).coefs = fit_cft;
                                % record rmse of the pixel
                                rec_cg(num_fc).rmse = rmse;
                                % treat first curve differently
                                if num_fc == rec_fc
                                    % record fit category 
                                    rec_cg(num_fc).category = 10 + min_num_c;
                                    % record break time
                                    rec_cg(num_fc).t_break = clrx(i_start);
                                    % record change probability
                                    rec_cg(num_fc).change_prob = 1;
                                else
                                    % record fit category 
                                    rec_cg(num_fc).category = 30 + min_num_c;
                                    % record break time
                                    rec_cg(num_fc).t_break = 0;
                                    % record change probability
                                    rec_cg(num_fc).change_prob = 0;
                                end
                                % record number of observations
                                rec_cg(num_fc).num_obs = i_start - i_break;
                                % record change magnitude
                                rec_cg(num_fc).magnitude = - median(v_dif_mag,1);
                                
                                % identified and move on for the next functional curve
                                num_fc = num_fc + 1;
                            
                            % enough to use median value
                            elseif i_start > i_break
                                % median value fit for the rest of the pixels < conse
                                % defining computed variables
                                fit_cft = zeros(max_num_c,nbands-1);
                                % rmse for each band
                                rmse = zeros(nbands-1,1);
                                
                                % median of the last < conse pixels
                                for i_B = 1:nbands-1
                                    fit_cft(1,i_B) = median(clry(i_break:i_start-1,i_B));
                                    rmse(i_B) = sqrt(mean((clry(i_break:i_start-1,i_B)-fit_cft(1,i_B)).^2));
                                end
                                
                                % record time of curve start
                                rec_cg(num_fc).t_start = clrx(i_break);
                                % record time of curve end
                                rec_cg(num_fc).t_end = clrx(i_start-1);
                                % record postion of the pixel
                                rec_cg(num_fc).pos = [N_row,N_col];
                                % record fitted coefficients
                                rec_cg(num_fc).coefs = fit_cft;
                                % record rmse of the pixel
                                rec_cg(num_fc).rmse = rmse;
                                % treat first curve differently
                                if num_fc == rec_fc
                                    % record fit category 
                                    rec_cg(num_fc).category = 10 + 1;
                                    % record break time
                                    rec_cg(num_fc).t_break = clrx(i_start);
                                    % record change probability
                                    rec_cg(num_fc).change_prob = (i_start - i_break)/conse;                                    
                                else
                                    % record fit category 
                                    rec_cg(num_fc).category = 30 + 1;
                                    % record break time
                                    rec_cg(num_fc).t_break = 0;
                                    % record change probability
                                    rec_cg(num_fc).change_prob = 0;
                                end
                                % record number of observations
                                rec_cg(num_fc).num_obs = i_start - i_break;
                                % record change magnitude
                                rec_cg(num_fc).magnitude = - median(v_dif_mag,1);
                                
                                % identified and move on for the next functional curves
                                num_fc = num_fc + 1;
                            end
                        end
                    end % end of "if num_c"
                end % end of "if v_dif > T_cg"
            end % end of initializing model
            
            % continuous monitoring started!!!
            if BL_train == 1 
                % all IDs
                IDs = i_start:i;
                i_span = i-i_start+1;
                
                % determine the time series model
                update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                
                % dynamic model fit when there are not many obs
                if  i_count == 0 || update_num_c < max_num_c
                    % update i_count at each interation
                    i_count = i-i_start+1;
                    
                    % defining computed variables
                    fit_cft=zeros(max_num_c,nbands-1);
                    % rmse for each band
                    rmse=zeros(nbands-1,1);
                    % record the diference in all bands
                    rec_v_dif = zeros(length(IDs),nbands-1);
                    
                    for i_B = 1:nbands-1
                        [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                            autoTSFit(clrx(IDs),clry(IDs,i_B),update_num_c);
                    end
                    
                    % adjusted mini RMSE based on p_min % of the mean ref
                    adj_rmse = zeros(nbands-1,1);
                    for i_B = B_detect
                        adj_rmse(i_B) = p_min*(fit_cft(1,i_B)+fit_cft(2,i_B)*(clrx(i_start)+clrx(i))/2);
                    end
                    % adjust agaist a constant
                    adj_rmse(adj_rmse < mini) = mini;
                    
                    adj_rmse = var_med;
                    adj_rmse(adj_rmse < mini) = mini;
                    % updating information for the first iteration
                    % record time of curve start
                    rec_cg(num_fc).t_start = clrx(i_start);
                    % record time of curve end
                    rec_cg(num_fc).t_end = clrx(i);
                    % record break time
                    rec_cg(num_fc).t_break = 0; % no break at the moment
                    % record postion of the pixel
                    rec_cg(num_fc).pos = [N_row,N_col];
                    % record fitted coefficients
                    rec_cg(num_fc).coefs = fit_cft;
                    % record rmse of the pixel
                    rec_cg(num_fc).rmse = rmse;
                    % record change probability
                    rec_cg(num_fc).change_prob = 0;
                    % record number of observations
                    rec_cg(num_fc).num_obs = i-i_start+1;
                    % record fit category
                    rec_cg(num_fc).category = 0 + update_num_c;
                    % record change magnitude
                    rec_cg(num_fc).magnitude = zeros(1,nbands-1);
                    
                    % detect change
                    % value of difference for conse obs
                    v_dif = zeros(conse,nbands-1);
                    % record the magnitude of change
                    v_dif_mag = v_dif;
                    vec_mag = zeros(conse,1);

                    for i_conse = 1:conse
                        for i_B = 1:nbands-1
                            % absolute difference
                            v_dif_mag(i_conse,i_B) = clry(i+i_conse,i_B)-autoTSPred(clrx(i+i_conse),fit_cft(:,i_B));
                            % normalized to z-scores
                            if sum(i_B == B_detect)
                                % minimum rmse
                                mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                
                                % z-scores
                                v_dif(i_conse,i_B) = (v_dif_mag(i_conse,i_B))/mini_rmse;
                            end
                        end
                        vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                    end
                    % IDs that haven't updated
                    IDsOld = IDs;
                else
                    if i-i_start+1 >= i_count + i_count/3
                        % update i_count at each interation year
                        i_count = i-i_start+1;
                        
                        % defining computed variables
                        fit_cft=zeros(max_num_c,nbands-1);
                        % rmse for each band
                        rmse=zeros(nbands-1,1);
                        % record the diference in all bands
                        rec_v_dif = zeros(length(IDs),nbands-1);
                        
                        for i_B = 1:nbands-1
                            [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                autoTSFit(clrx(IDs),clry(IDs,i_B),update_num_c);
                        end
                        
                        % adjusted mini RMSE based on p_min % of the mean ref
                        adj_rmse = zeros(nbands-1,1);
                        for i_B = B_detect
                            adj_rmse(i_B) = p_min*(fit_cft(1,i_B)+fit_cft(2,i_B)*(clrx(i_start)+clrx(i))/2);
                        end
                        % adjust agaist a constant
                        adj_rmse(adj_rmse < mini) = mini;
                        
                        adj_rmse = var_med;
                        adj_rmse(adj_rmse < mini) = mini;
                        
                        % record fitted coefficients
                        rec_cg(num_fc).coefs = fit_cft;
                        % record rmse of the pixel
                        rec_cg(num_fc).rmse = rmse;
                        % record number of observations
                        rec_cg(num_fc).num_obs = i-i_start+1;
                        % record fit category
                        rec_cg(num_fc).category = 0 + update_num_c;
                        
                        % IDs that haven't updated
                        IDsOld = IDs;
                    end
                    
                    % record time of curve end
                    rec_cg(num_fc).t_end=clrx(i);
                    
                    % use temporally-adjusted RMSE
                    if length(IDsOld) <= n_times*max_num_c
                        % number of observations for calculating RMSE
                        n_rmse = length(IDsOld);
                        tmpcg_rmse = rmse;
                    else
                        % use fixed number for RMSE computing
                        n_rmse = n_times*max_num_c;
                        tmpcg = zeros(nbands-1,1);
                        % better days counting for RMSE calculating
                        % relative days distance
                        d_rt = clrx(IDsOld) - clrx(i+conse);
                        d_yr = abs(round(d_rt/num_yrs)*num_yrs-d_rt);
                        
                        [~,sorted_indx] = sort(d_yr);
                        sorted_indx = sorted_indx(1:n_rmse);
                        
                        for i_B = B_detect
                            % temporally changing RMSE
                            tmpcg_rmse(i_B) = norm(rec_v_dif(IDsOld(sorted_indx)-IDsOld(1)+1,i_B))/...
                                sqrt(n_rmse-update_num_c);
                        end
                    end
                    
                    % move the ith col to i-1th col
                    v_dif(1:conse-1,:) = v_dif(2:conse,:);
                    % only compute the difference of last consecutive obs
                    v_dif(conse,:) = 0;
                    % move the ith col to i-1th col
                    v_dif_mag(1:conse-1,:) = v_dif_mag(2:conse,:);
                    % record the magnitude of change of the last conse obs
                    v_dif_mag(conse,:) = 0;
                    % move the ith col to i-1th col
                    vec_mag(1:conse-1) = vec_mag(2:conse);
                    % change vector magnitude
                    vec_mag(conse) = 0;
                    
                    for i_B = 1:nbands-1
                        % absolute difference for all bands
                        v_dif_mag(conse,i_B) = clry(i+conse,i_B)-autoTSPred(clrx(i+conse),fit_cft(:,i_B));
                        % normalized to z-scores
                        if sum(i_B == B_detect)
                            % minimum rmse
                            mini_rmse = max(adj_rmse(i_B),tmpcg_rmse(i_B));
                            
                            % z-scores
                            v_dif(end,i_B) = (v_dif_mag(end,i_B))/mini_rmse;
                        end
                    end
                    vec_mag(conse) = norm(v_dif(end,B_detect))^2;
                end

                % change detection
                if min(vec_mag) > T_cg % change detected
                    fprintf('Change Magnitude = %.2f\n',min(vec_mag) - T_cg);                    
                    % record break time
                    rec_cg(num_fc).t_break = clrx(i+1);
                    % record change probability
                    rec_cg(num_fc).change_prob = 1;
                    % record change magnitude
                    rec_cg(num_fc).magnitude = median(v_dif_mag,1);
                    
                    % identified and move on for the next functional curve
                    num_fc = num_fc + 1;
                    % start from i+1 for the next functional curve
                    i_start = i + 1;
                    % start training again
                    BL_train = 0;
                    
                elseif vec_mag(1) > Tmax_cg % false change
                    % remove noise
                    clrx(i+1,:)=[];
                    clry(i+1,:)=[];
                    i=i-1; % stay & check again after noise removal
                end
            end % end of continuous monitoring
        end % end of checking basic requrirements
        
        % move forward to the i+1th clear observation
        i=i+1;
    end % end of while iterative
    
    % Two ways for processing the end of the time series
    if BL_train == 1
        % 1) if no break find at the end of the time series
        % define probability of change based on conse
        for i_conse = conse:-1:1
            if vec_mag(i_conse) <= T_cg
                % the last stable id
                id_last = i_conse;
                break;
            end
        end
        
        % update change probability
        rec_cg(num_fc).change_prob = (conse-id_last)/conse;
        % update end time of the curve
        rec_cg(num_fc).t_end=clrx(end-conse+id_last);
        
        if conse > id_last % > 1
            % update time of the probable change
            rec_cg(num_fc).t_break = clrx(end-conse+id_last+1);
            % update magnitude of change
            rec_cg(num_fc).magnitude = median(v_dif_mag(id_last+1:conse,:),1);
            
            % defining computed variables
            fit_cft=zeros(max_num_c,nbands-1);
            % rmse for each band
            rmse=zeros(nbands-1,1);
            % median of the last < conse pixels
            for i_B=1:nbands-1
                fit_cft(1,i_B) = median(clry(end-conse+id_last+1:end,i_B));
                rmse(i_B) = sqrt(mean((clry(end-conse+id_last+1:end,i_B)-fit_cft(1,i_B)).^2));
            end
            
            % identified and move on for the next functional curve
            num_fc = num_fc+1;
            
            % record time of curve start
            rec_cg(num_fc).t_start = clrx(end-conse+id_last+1);
            % record time of curve end
            rec_cg(num_fc).t_end = clrx(end);
            % record break time
            rec_cg(num_fc).t_break = 0;
            % record postion of the pixel
            rec_cg(num_fc).pos = [N_row,N_col];
            % record fitted coefficients
            rec_cg(num_fc).coefs = fit_cft;
            % record rmse of the pixel
            rec_cg(num_fc).rmse = rmse;
            % record change probability
            rec_cg(num_fc).change_prob = 0;
            % record number of observations
            rec_cg(num_fc).num_obs = conse-id_last;
            % record fit category
            rec_cg(num_fc).category = 20 + 1; 
            % record change magnitude
            rec_cg(num_fc).magnitude = zeros(1,nbands-1);
        end
        
    elseif BL_train == 0
        % 2) if break find close to the end of the time series
        % Use [conse,min_num_c*n_times+conse) to fit curve
        % update i_start
        if num_fc == rec_fc
            % first curve
            i_start = 1;
        else
            i_start = find(clrx > rec_cg(num_fc-1).t_end);
            i_start = i_start(1);
        end
                            
        % multitemporal cloud mask
        blIDs = autoMask(clrx(i_start:end),clry(i_start:end,[num_B1,num_B2]),...
            (clrx(end)-clrx(i_start))/num_yrs,T_const);
        
        % update i_span after noise removal
        i_span = sum(~blIDs);
        
        IDs = i_start:length(clrx); % all IDs
        rmIDs = IDs(blIDs(1:end-conse) == 1); % IDs to be removed

        % remove noise pixels between i_start & i
        clrx(rmIDs) = [];
        clry(rmIDs,:) = [];
        
        % defining computed variables
        fit_cft = zeros(max_num_c,nbands-1);
        % rmse for each band
        rmse = zeros(nbands-1,1);
        for i_B = 1:nbands-1
            [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx(i_start:end),clry(i_start:end,i_B),min_num_c);
        end
        
        % record time of curve start
        rec_cg(num_fc).t_start = clrx(i_start);
        % record time of curve end
        rec_cg(num_fc).t_end = clrx(end);
        % record break time
        rec_cg(num_fc).t_break = 0;
        % record postion of the pixel
        rec_cg(num_fc).pos = [N_row,N_col];
        % record fitted coefficients
        rec_cg(num_fc).coefs = fit_cft;
        % record rmse of the pixel
        rec_cg(num_fc).rmse = rmse;
        % record change probability
        rec_cg(num_fc).change_prob = 0;
        % record number of observations
        rec_cg(num_fc).num_obs = i_span;
        % record fit category
        rec_cg(num_fc).category = 20 + min_num_c; % simple model fit at the end
        % record change magnitude
        rec_cg(num_fc).magnitude = zeros(1,nbands-1);
    end
end
% toc
% profile viewer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start of plot
% Plot raw data (after Fmask screening)
figure;
set(gcf,'Position',[0 300*nbands 3000 300]);
set(gca,'FontSize',13)
plot(clrx,clry)
legend('Band 1','Band 2','Band 3','Band 4','Band 5','Band 7','Band 6','location','NorthEastOutside');
datetick('x',29,'keeplimits');

Band_Plot=B_detect;
for B_plot=Band_Plot
    figure;
    set(gcf,'Position',[0 300*(B_plot-2) 3000 300]);
    set(gca,'FontSize',13)
    
    plot(Xs,Ys(:,B_plot),'b.');
    hold on;
    
    % plot clear pixels
    plot(clrx,clry(:,B_plot),'k.');
    
    hold on;
    % plot curves
    num_fit=size(rec_cg,2);
    
    m=5;
    color_tab=jet(m);
    for i=1:num_fit
        x_plot=rec_cg(i).t_start:rec_cg(i).t_end;
        pred_y=autoTSPred(x_plot',rec_cg(i).coefs(:,B_plot));
        icolor=i;
        if icolor > m
            icolor = mod(icolor,m)+1;
        end
        
        if rec_cg(i).category < 60
            plot(x_plot,pred_y,'-','color',color_tab(icolor,:));
        else
            % plot(x_plot,pred_y,'k-');
        end
        
        hold on;
        if rec_cg(i).change_prob == 1
%              if rec_cg(i).category < 10
                plot(rec_cg(i).t_break,Ys(Xs==rec_cg(i).t_break,B_plot),'ro','Markersize',10);
%              else
%                  plot(rec_cg(i).t_break,Ys(Xs==rec_cg(i).t_break,B_plot),'g^','Markersize',10);
%              end
            %             x_time = datenum(1996,7,31);
            %             plot(x_time,Ys(Xs==x_time,B_plot),'bo','Markersize',20);
            ntime=datevec(rec_cg(i).t_break);
            ndays = datenum(rec_cg(i).t_break)-datenum(ntime(1),0,0);
            hold on;
            plot(rec_cg(i).t_break,Ys(Xs==rec_cg(i).t_break,B_plot),'k.');
            if B_plot==Band_Plot(end)
                fprintf('The %dth change occurred at %d %d\n',i,ntime(1),ndays);
            end
%         elseif rec_cg(i).change_prob ~= 0
%             plot(rec_cg(i).t_break,Ys(Xs==rec_cg(i).t_break,B_plot),'b^','Markersize',10);
        end
    end
    B_num=B_plot;
    if B_plot==6
        B_num=7;
    elseif B_plot==7
        B_num=6;
    end
    title(['Landsat Band ',num2str(B_num),' at row/col=',num2str(N_row),'/',num2str(N_col)]);
    if B_num == 6
        ylabel(['Band ',num2str(B_num),' Brightness Temperature (^oCX10^2)']);
    else
        ylabel(['Band ',num2str(B_num),' Surface Reflectance (X10^4)']);
    end
    % axis([datenum(1985,1,1) datenum(2014,12,31) 0 4000]);
    datetick('x',29,'keeplimits');
    if rec_cg(1).change_prob < 1
        legend('Ephermal Change','Clear Observations','Model Fitting','location','NorthEastOutside');
    else
        legend('Ephermal Change','Clear Observations','Model Fitting','Change','location','NorthEastOutside');
    end
end
% length(Xs)
% length(clrx)
pause;
close all;
fprintf('CCDC v12.21 Done!\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of plot
% end of function
end