% Collect altimeter data and interpolate to standard and linearly-spaced tracks
% Arin Nelson
% Rewritten 12/12/2020
%==========================================================================
clc; close all;

% Switches
Switch     = zeros(9,1);
Switch(01) = 1;     % Collect raw altimeter data
Switch(02) = 1;     % Interpolate to standard-track and linear-track points

% Directories
dir_data = '.\Data';
dir_altm = 'F:\Datasets\Altimetry\Jason_StandardCorrections_fromEd';
str_altm = {'j2','j3'};         % Altimeter missions
str_smpl = {'track','lin'};     % Sampling points

% Options
do_debug    = 1;                        % Show debug plots and information
rE          = 6371;                     % Radius of Earth
time_offset = datenum('1985-01-01');    % Time offset of altimeter data
bath_min    = 500;               % Only allow depths >= this (m) (Zaron 2019: 500)
dshore_min  = 12;                % Only allow distances from shore >= this (km) (Zaron 2019: 12)
basins      = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21];
min_gap     = 30;
min_length  = 150;
num_min     = [18 111];

% Constants
n_track = 254;
n_altm  = numel(str_altm);
n_smpl  = numel(str_smpl);

% NOTES
% Min num to separate semidiurnal and diurnal constituents: 18
% Min num to separate all constituents: 111

%==========================================================================
if(Switch(01))
    
  % Loop through altimeter missions
  for ia=1:n_altm
 
    % Initialize containers
     lon = cell(n_track,1);     % Longitudes (^oE)
     lat = cell(n_track,1);     % Latitudes (^oN)
    time = cell(n_track,1);     % Times (datenum)
     num = cell(n_track,1);     % Total # of samples at this pt
     sla = cell(n_track,1);     % sea level anomaly (cm)
     
    % Load data
    for it=1:n_track
    clc; disp(['Loading altimeter data for mission ' num2str(ia) ', track ' num2str(it) '...']);
    
      % Construct file name  
      path_data = [dir_altm '\' str_altm{ia} 'acolin\' str_altm{ia} 'acolin_p' sprintf('%0.4d',it) '.nc'];
  
      % Read in data
       lon{it} = ncread(path_data,'lon');  
       lat{it} = ncread(path_data,'lat');
      time{it} = ncread(path_data,'time');
       sla{it} = ncread(path_data,'sla');
      clear path_data;

      % Compute total number of cycles with valid observations
      num{it} = sum(~isnan(sla{it}),1);
      
      % Set time units to datenum
      time{it} = (time{it}./86400) + time_offset;
     
      % Convert sea level anomaly from units of m to cm
      sla{it} = sla{it}.*100;
      
    end
    clear it;
    
    % Debug check
    if(do_debug==1)
      figure('units','normalized','outerposition',[0.2 0 0.6 1]);
        xx = [];    for it=1:n_track;   xx = [xx(:); nanmean(lon{it},1)'];   end
        yy = [];    for it=1:n_track;   yy = [yy(:); nanmean(lat{it},1)'];   end
        nn = [num{:}];
        zm = [];    for it=1:n_track;   zm = [zm(:); nanmean(sla{it},1)'];   end
        zs = [];    for it=1:n_track;   zs = [zs(:); nanstd(sla{it},0,1)'];   end
        ii = 5:10:numel(xx);
        subplot(3,1,1); scatter(xx(ii),yy(ii),2,zm(ii),'o','filled'); axis tight; box on; colormap(jet); colorbar; title('sla mean');
        subplot(3,1,2); scatter(xx(ii),yy(ii),2,zs(ii),'o','filled'); axis tight; box on; colormap(jet); colorbar; title('sla std'); 
        subplot(3,1,3); scatter(xx(ii),yy(ii),2,nn(ii),'o','filled'); axis tight; box on; colormap(jet); colorbar; title('num');
      print(gcf,'-dpng',['Plots/Debug/Debug_Step02-01_' str_altm{ia} '_raw.png']);
      close all; clear xx yy nn zm zs ii
    end
    
    % Save
    save([dir_data '\raw_data_' str_altm{ia} '.mat'],'lon','lat','time','sla','num');
     
    % Clean-up
    clear lon lat time sla num;
    
  end
  clear ia;

end
%==========================================================================
if(Switch(02))
    
  % Loop through altimeter missions
  for ia=1:n_altm
  for is=1:n_smpl
      
    % Load info
     raw_data = load([dir_data '\raw_data_'   str_altm{ia} '.mat']);
    smpl_info = load([dir_data '\' str_smpl{is} '_info_' str_altm{ia} '.mat']);
    n_point   = size(smpl_info.lon,2);  
    
    % Init new variables
    time = cell(n_track,1);
    sla  = cell(n_track,1);
    %num  = cell(n_track,1);
    num  = zeros(n_track,n_point);
    
    % Variables that will be appended to raw data file
    dist_raw = cell(n_track,1);
    
    % Loop through tracks
    for it=1:n_track
    clc; disp(['On altimeter mission ' num2str(ia) ', Sampling ' num2str(is) ', Track #' num2str(it) '...']);
    
      % This track's data
       lon_raw = raw_data.lon{it};
       lat_raw = raw_data.lat{it};
      time_raw = raw_data.time{it}; % NOTE: Units of datenum
       sla_raw = raw_data.sla{it};  % NOTE: Units of cm
       num_raw = raw_data.num{it};
        
      % Distance vectors on new grids
      lon_0     = smpl_info.lon(it,1);
      lat_0     = smpl_info.lat(it,1);
      if(strcmp(str_smpl{is},'track'))
        dist_smpl = smpl_info.dist(it,:); 
      elseif(strcmp(str_smpl{is},'lin'))
        dist_smpl = smpl_info.dist(:); 
      else
        error('unkonwn str_smpl!');
      end
      
      % And land-sea masks
      mask_smpl = smpl_info.mask(it,:);   imask = find(mask_smpl==1);
      
      % Some useful constants
      nc = size(sla_raw,1);
      nr = size(sla_raw,2);
      np = numel(dist_smpl);
      
      % Compute along-track distances
      dist_raw{it} = NaN(nc,nr);
      for ic=1:nc
      for ir=1:nr
        dist_raw{it}(ic,ir) = rE.*circledist(lon_0,lat_0,lon_raw(ic,ir),lat_raw(ic,ir));  
      end
      end
      clear ic ir;
      
      % Interpolate to new track points
      time{it} = NaN(nc,np);
       sla{it} = NaN(nc,np);
       %num{it} = NaN(1,np);
      for ic=1:nc
        if(sum(~isnan(time_raw(ic,:)))>1)	
          time{it}(ic,imask) = naninterp1(dist_raw{it}(ic,:),time_raw(ic,:),dist_smpl(imask),'linear');
        end
        if(sum(~isnan( sla_raw(ic,:)))>1)	
           sla{it}(ic,imask) = naninterp1(dist_raw{it}(ic,:),sla_raw(ic,:), dist_smpl(imask),'linear');
        end
        if(ic==1)
          %num{it}(1,imask) = naninterp1(nanmean(dist_raw{it},1),num_raw, dist_smpl(imask),'linear');
          num(it,imask) = naninterp1(nanmean(dist_raw{it},1),num_raw,dist_smpl(imask),'linear');
        end
      end
      clear ic;
      
      % Clean-up
      clear lon_raw lat_raw time_raw sla_raw num_raw lon_0 lat_0 dist_smpl mask_smpl imask nc nr np;
         
    end
    clear it;
    
    % Debug plot
    if(do_debug)
        
      % Figure #1: track   
      figure('units','normalized','outerposition',[0 0 1 1]);
        xx = smpl_info.lon';   xx = xx(:);
        yy = smpl_info.lat';   yy = yy(:);
        ii = 5:10:numel(xx);
        subplot(2,3,1); % mask
          zz = smpl_info.mask';      
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('mask');
        subplot(2,3,2); % bath
          zz = smpl_info.bath';      
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('bath');
        subplot(2,3,3); % dshore
          zz = smpl_info.dshore';	
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('dshore');
        subplot(2,3,4); % num
          zz = num';  zz = zz(:);    
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('num');
        subplot(2,3,5); % sla mean
          zz = []; for it=1:n_track; zz = [zz(:); nanmean(sla{it},1)']; end
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('sla mean');
        subplot(2,3,6); % sla std
          zz = []; for it=1:n_track; zz = [zz(:); nanstd(sla{it},0,1)']; end
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('sla std');
      print(gcf,'-dpng',['Plots/Debug/Debug_Step02-02_' str_altm{ia} '_' str_smpl{is} '.png']);
      close all; clear xx yy zz ii; close all;  
    end

    % Save (track)
    save([dir_data '\' str_smpl{is} '_data_' str_altm{ia} '.mat'],'time','sla'); %,'num');
    
    % Save (append)
    dist = dist_raw;
    save([dir_data '\raw_data_' str_altm{ia} '.mat'],'dist','-append');
    save([dir_data '\' str_smpl{is} '_info_' str_altm{ia} '.mat'],'num','-append');
    
    % Clean-up
    clear raw_data smpl_info time sla num dist dist_raw;
    
  end
  end
  clear ia is;
  
end
