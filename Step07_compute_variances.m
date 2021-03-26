% Tests
%==========================================================================

% Switches
Switch     = zeros(9,1);
Switch(01) = 1;     % Do

% Directories
addpath('Functions');
dir_data   = '.\Data';
dir_mdl    = 'F:\Datasets\HYCOM\PhaseMatch';
str_mdl    = {'NonDA_2016','DA_2016','DA_Modern'};
str_altm   = {'j2','j3'};
str_mthd   = {'track','lin'};
str_tide   = {'M2','S2','N2','K1','O1'};

% Options
do_debug   = 1;  
ref_m_to_a = [1 1 2];
ref_a_to_m = [2 3];
tref       = datenum('2016-01-01');

% Binning
bin_dx = 5;
bin_dy = 5;
bin_lon_edge = 0:bin_dx:360;
bin_lat_edge = -65:bin_dy:65;
bin_lon_cntr = (bin_lon_edge(2:end)+bin_lon_edge(1:end-1))./2;
bin_lat_cntr = (bin_lat_edge(2:end)+bin_lat_edge(1:end-1))./2;

% Constants
n_track = 254;                      % # altimeter tracks
n_mdl   = numel(str_mdl);           % # model simulations
n_altm  = numel(str_altm);          % # altimeter missions
n_mthd  = numel(str_mthd);
n_binx  = numel(bin_lon_cntr);
n_biny  = numel(bin_lat_cntr);

%==========================================================================
if(Switch(01))
    
  % Loop through model outputs and sampling methods
  for im=2:n_mdl
  for id=2 %1:n_mthd
      
    % Data
    load([dir_data '\' str_mthd{id} '_info_' str_altm{ref_m_to_a(im)} '.mat']);
    n_point = size(lon,2);
    
    % Altimeter data
    load([dir_data '\' str_mthd{id} '_data_' str_altm{ref_m_to_a(im)} '.mat'],'time','sla');
    altm_time = time;   clear time;
    
    % Tide fits
    load([dir_data '\' str_mthd{id} '_' str_mdl{im} '_tidefit.mat']);
    mdl_tide = tide_fit;    clear tide_fit;  
    
    % What will be computed
    v_z_t	= NaN(n_track,n_point,8);   n_z_t	= NaN(n_track,n_point);
    v_m_t   = NaN(n_track,n_point,8);   n_m_t   = NaN(n_track,n_point);
    v_m_f   = NaN(n_track,n_point,4);   n_m_f   = NaN(n_track,n_point);
    s2      = NaN(n_track,n_point,7);   d2      = NaN(n_track,n_point,7);
    v_a     = NaN(n_track,n_point);     n_a     = NaN(n_track,n_point);
    v_n_z_t = NaN(n_track,n_point,8);   n_n_z_t = NaN(n_track,n_point);
    v_n_m_t = NaN(n_track,n_point,8);   n_n_m_t = NaN(n_track,n_point);
    v_n_m_f = NaN(n_track,n_point,4);   n_n_m_f = NaN(n_track,n_point);

    % Loop through tracks
    for it=1:n_track
    clc; disp(['On model output #' num2str(im) ', sampling method #' num2str(id) ', track #' num2str(it) '...']);
    
      % This track's model time series file
      f_mdl = [dir_mdl '\' str_mdl{im} '\_Processed\MATLAB\Track' sprintf('%0.3d',it) '.mat'];
      load(f_mdl,'time',[str_mthd{id} '_*']);
      mdl_time = time;      clear time;
      
      % Renaming data
      eval(['dat_filt = ' str_mthd{id} '_dat_filt; clear ' str_mthd{id} '_dat_filt;']);
      eval(['dat_tide_timeseries = ' str_mthd{id} '_tide_timeseries; clear ' str_mthd{id} '_tide_timeseries;']);

      % Loop through points
      parfor ip=1:n_point
      if(is_valid(it,ip))
          
        %-----------------------------------------------------------------%
        % PARALLELIZABLE PIECE
        
        [v_z_t(it,ip,:),   n_z_t(it,ip),   ...
         v_m_t(it,ip,:),   n_m_t(it,ip),   ...
         v_m_f(it,ip,:),   n_m_f(it,ip),   ...
         d2(it,ip,:),      s2(it,ip,:),  ...
         v_a(it,ip),       n_a(it,ip),     ...
         v_n_z_t(it,ip,:), n_n_z_t(it,ip), ...
         v_n_m_t(it,ip,:), n_n_m_t(it,ip), ...
         v_n_m_f(it,ip,:), n_n_m_f(it,ip)  ...
        ] = var_reduc_tests_parpiece( squeeze(tide_zaron_re(it,ip,:)), ...
                                      squeeze(tide_zaron_im(it,ip,:)), ...
                                       altm_time{it}(:,ip), ...
                                       tref , ...
                                       mdl_time, ...
                                       squeeze(dat_tide_timeseries(:,ip,:)), ...
                                       squeeze(dat_filt(:,ip,:)), ...
                                       squeeze(mdl_tide(it,ip,:,:)), ...
                                       sla{it}(:,ip) ...
                                     );
        
        % END PARALLELIZABLE PIECE
        %-----------------------------------------------------------------%
      
      end
      end
      clear ip;
      
      % Debug check
      
      
      % Clean-up
      clear f_mdl mdl_time dat_filt dat_tide_timeseries;
      
    end
    clear it;
    
    % Debug plot

    
    % Save
    save(['Data/Variances_' str_mdl{im} '_' str_mthd{id} '.mat'],...
         'v_*','s2','d2','n_z_t','n_m_*','n_a','n_n_*');
    
    % Clean-up
    clear var_* num_*;
    
  end
  end
  clear im id;
    
end
%==========================================================================