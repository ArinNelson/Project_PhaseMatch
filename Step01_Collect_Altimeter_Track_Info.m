% Collect altimeter standard tracks and information along them
% Arin Nelson
% 12/12/2020
%==========================================================================
clc; close all;

% Switches
Switch     = zeros(9,1);
Switch(01) = 0;         % Gather altimeter track & sampling point info
Switch(02) = 0;         % Generate linearly-spaced-points tracks
Switch(03) = 0;         % Interpolate bathymetry & land-sea mask from GEBCO
Switch(04) = 0;         % Interpolate distance-from-coast from GEBCO
Switch(05) = 0;         % Separate info's into separate files per mission
Switch(06) = 0;         % Interpolate basin indices to track points
Switch(07) = 0;         % Link model points to track points

% Directories & Files
dir_info     = '\Datasets\PODAAC\';
file_bath    = '\Datasets\HYCOM\depth_GLBc0.04_27.a';
file_dshore  = '\Datasets\GEBCO\GEBCO_gridone_distance_to_coastline_v2.0m.nc';
file_mdlgrid = '\Datasets\HYCOM/HYCOM25_Grid.mat';
file_basins  = '\Datasets\NOAA\range_area.msk';
dir_mdl      = '\Datasets\HYCOM\PhaseMatch';
str_mdl      = {'DA_2016','DA_Modern'};

% Options
do_debug      = 1;              % Show debug plots
rE            = 6371;           % Radius of Earth (km)
str_altm      = {'j2','j3'};	  % Altimeter labels (Jason-2 & -3)
dist_lin_step = 6.6;            % Step size of linearly-spaced-points (km)

% Constants
n_altm = numel(str_altm);
  
%==========================================================================
if(Switch(01))
    
  % Loading reference data
  eqcrosses2   = importdata([dir_info 'eqxl_tp_then_j1.txt']);
  eqcrosses3   = importdata([dir_info 'eqxl_j3.txt']);
  reftrack_asc = importdata([dir_info 'reftrk_ascending.txt']);
  reftrack_dsc = importdata([dir_info 'reftrk_descending.txt']);

  % Total number of tracks
  n_track = numel(eqcrosses2);
  
  % Parse data
  lat_asc = reftrack_asc(:,1);
  lon_asc = reftrack_asc(:,2);
  lat_dsc = reftrack_dsc(:,1);
  lon_dsc = reftrack_dsc(:,2);
  
  % Total number of sample points per track
  n_point = numel(lat_asc);
  
  % Compute track longitudes and latitudes
  lon = zeros(n_track,n_point,2);
  lat = zeros(n_track,n_point);
  for it=1:n_track
    if(mod(it,2)~=0)
      lon(it,:,1) = lon_asc + eqcrosses2(it);
      lon(it,:,2) = lon_asc + eqcrosses3(it);
      lat(it,:)   = lat_asc;
    else
      lon(it,:,1) = lon_dsc + eqcrosses2(it);
      lon(it,:,2) = lon_dsc + eqcrosses3(it);
      lat(it,:)   = lat_dsc;
    end
  end
  clear it;
  
  % Limit points to 0-360
  lon(lon < 0  ) = lon(lon < 0  ) + 360;
  lon(lon > 360) = lon(lon > 360) - 360;
  
  % Along-track distance
  dist = zeros(n_track,n_point,2);
  for it=1:n_track
  for ia=1:n_altm
  for ip=1:n_point    
    dist(it,ip,ia) = rE.*circledist(lon(it,1,ia),lat(it,1),lon(it,ip,ia),lat(it,ip));
  end
  end
  end
  clear it ia ip;
  
  % Debug plot: map of distance-along-track's
  if(do_debug)
    figure('units','normalized','outerposition',[0 0 1 1]);
      for ia=1:n_altm
        subplot(1,2,ia); hold on;
          for it=1:n_track
            scatter(lon(it,:,ia),lat(it,:),2,dist(it,:,ia),'o','filled'); 
          end
        hold off; axis tight; title(['Altimeter sampling info (' str_altm{ia} ')']);
        colormap(jet); cb=colorbar; set(get(cb,'ylabel'),'string','along-track distance (km)');
      end
    print(gcf,'-dpng','Plots/Debug/Debug_Step01-01_AltimSamplingInfo.png');
    close all;  clear ia it cb;
  end

  % Save?
  save('Data/track_info_raw.mat','lon','lat','dist','str_altm');
  
  % Clean-up
  clear eqcrosses* reftrack_* lon* lat* dist

end
%==========================================================================
if(Switch(02))
  
  % Raw track info
  track_info = load('Data/track_info_raw.mat');
    
  % Generate evenly-spaced-points distance-along-track values
  dist_start = 0;
  dist_end   = max(max(track_info.dist(:,end,:)));
  dist       = dist_start : dist_lin_step : dist_end;
  n_lin      = numel(dist);
  
  % Interpolate to find the lons and lats of these tracks
  lon = zeros(n_track,n_lin,2);
  lat = zeros(n_track,n_lin);
  for it=1:n_track
  for ia=1:n_altm
      
    % This track's points
    xx = track_info.lon(it,:,ia);
    yy = track_info.lat(it,:);
    
    % Unwrap looping longitudes
    xx = unwrap( xx.*(pi/180) ).*(180/pi);
      
    % Complex position
    pos_track = xx + sqrt(-1).*yy;
    
    % Interpolate to linearly-spaced-points
    pos_lin = interp1(track_info.dist(it,:,ia),pos_track,dist,'spline');
    
    % Re-loop longitudes
    x_lin = real(pos_lin);  x_lin(x_lin<0) = x_lin(x_lin<0)+360;    x_lin(x_lin>360) = x_lin(x_lin>360)-360;
    y_lin = imag(pos_lin);
    
    % Save the points
    lon(it,:,ia) = x_lin;
    lat(it,:)    = y_lin;
    
    % Clean-up
    clear xx yy pos_track pos_lin x_lin y_lin;
    
  end
  end
  clear it ia;
  
  % Debug plots
  if(do_debug)
    figure('units','normalized','outerposition',[0 0 1 1]);
      hold on;
        for it=1 %:n_track
        for ia=1 %:n_altm
          scatter(track_info.lon(it,:,ia),track_info.lat(it,:),8,track_info.dist(it,:,ia),'d','filled'); 
          scatter(lon(it,:,ia),lat(it,:),8,dist,'s','filled'); 
        end
        end
      hold off; colormap(jet); cb=colorbar; set(get(cb,'ylabel'),'string','along-track distance (km)');
    print(gcf,'-dpng','Plots/Debug/Debug_Step01-02_LinearPointInterp-j2-track001.png');
    close all;  clear it ia cb;
  end
  
  % Save?
  save('Data/lin_info_raw.mat','lon','lat','dist','str_altm');
    
  % Clean-up
  clear track_info dist_start dist_end dist lon lat
  
end
%==========================================================================
if(Switch(03))

  % Load in needed data
  track_info = load('Data/track_info_raw.mat');
    lin_info = load('Data/lin_info_raw.mat');
    
  % Dimensions  
  n_track    = size(track_info.lon,1);
  n_point    = size(track_info.lon,2);
  n_lin      = size(  lin_info.lon,2);
    
  % Read in HYCOM lon and lat
  tmp = load(file_mdlgrid,'lon','lat');
  lon_mdl = tmp.lon';
  lat_mdl = tmp.lat';
  clear tmp

  % Read in HYCOM-bathymetry-derived-from-GEBCO file
  nt  = numel(lon_mdl);
  fid = fopen(file_bath,'r','b');
  tmp = fread(fid,[1,nt],'real*4');
  bath_gebco = reshape(tmp,size(lon_mdl)); 
  fclose(fid); clear fid tmp nt;
  
%   % Check
%   zz = bath_gebco(5:10:end,5:10:end);   zz(zz>1e10) = NaN;
%   surf(lon_mdl(5:10:end,5:10:end),lat_mdl(5:10:end,5:10:end),zz,'edgecolor','none'); view(2); colormap(jet); colorbar;

  % Don't need so many points?
  xx = lon_mdl; %(2:3:end,2:3:end);
  yy = lat_mdl; %(2:3:end,2:3:end);
  zz = bath_gebco; %(2:3:end,2:3:end);
  ntrplnt = scatteredInterpolant(double(xx(:)),double(yy(:)),zz(:),'linear');
  clear xx yy zz;
  
  % Ensure lon's are within HYCOM's range (74.14 - 434.14)
  xa = track_info.lon; xa(xa<74.14) = xa(xa<74.14)+360;
  xl =   lin_info.lon; xl(xl<74.14) = xl(xl<74.14)+360;
  
  % Interpolate bathymetry to track points
  bath_track = NaN(n_track,n_point,2);
  bath_lin   = NaN(n_track,n_lin,2);
  for ia=1:n_altm
    bath_track(:,:,ia) = ntrplnt(xa(:,:,ia),track_info.lat);
    bath_lin(:,:,ia)   = ntrplnt(xl(:,:,ia),  lin_info.lat);
  end
  clear ia xa xl;
  
  % Set land values to NaN
  maxb = max(bath_gebco(bath_gebco < 1e10));
  bath_track(bath_track > maxb) = NaN;
  bath_lin(  bath_lin > maxb)   = NaN;
  clear maxb;
  
  % Generate land-sea masks
  mask_track = ~isnan(bath_track);
  mask_lin   = ~isnan(bath_lin  );
  
  % Debug plots
  if(do_debug)
    figure('units','normalized','outerposition',[0 0 1 1]);
      subplot(2,1,1);
        xx = track_info.lon(:,:,1);  xx = xx(:);
        yy = track_info.lat(:);
        ii = 2:3:numel(xx);
        zz = bath_track(:,:,1); zz = zz(:);
        scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); caxis([0 7000]); title('bath track');
      subplot(2,1,2);
        xx = lin_info.lon(:,:,1);  xx = xx(:);
        yy = lin_info.lat(:);
        ii = 2:3:numel(xx);
        zz = bath_lin(:,:,1); zz = zz(:);
        scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); caxis([0 7000]); title('bath lin');
     cb=colorbar; set(get(cb,'ylabel'),'string','bathymetry (m)');
                  print(gcf,'-dpng','Plots/Debug/Debug_Step01-03_AltimPointBathymetry.png');
    close all;  clear xx yy ii zz cb;
  end
      
  % Save
  bath = bath_track;    mask = mask_track;  save('Data/track_info_raw.mat','bath','mask','-append');    clear bath mask;
  bath = bath_lin;      mask = mask_lin;    save('Data/lin_info_raw.mat',  'bath','mask','-append');    clear bath mask;
    
  % Clean-up
  clear *_info ntrplnt bath_* mask_* *_mdl;
  
end
%==========================================================================
if(Switch(04))
    
  % Load track info's
  track_info = load('Data/track_info_raw.mat');
    lin_info = load('Data/lin_info_raw.mat');
    
  % Dimensions  
  n_track = size(track_info.lon,1);
  n_point = size(track_info.lon,2);
  n_lin   = size(  lin_info.lon,2);
    
  % Load GEBCO distance-from-shore data
     lon_gebco = ncread(file_dshore,'longitude');
     lat_gebco = ncread(file_dshore,'latitude');
  dshore_gebco = ncread(file_dshore,'dist_to_coastline');
  
  % Convert dshore_gebco to kilometers
  dshore_gebco = dshore_gebco./1000;
  
%   % Don't need so many points?
%      lon_gebco =    lon_gebco(2:3:end,2:3:end);
%      lat_gebco =    lat_gebco(2:3:end,2:3:end);
%   dshore_gebco = dshore_gebco(2:3:end,2:3:end);
    
  % Interpolate dshore to track points
  dshore_track = NaN(n_track,n_point,2);
  dshore_lin   = NaN(n_track,n_lin,2);
  for ia=1:n_altm
    dshore_track(:,:,ia) = interp2(lon_gebco,lat_gebco,dshore_gebco',track_info.lon(:,:,ia),track_info.lat,'linear');
    dshore_lin(:,:,ia)   = interp2(lon_gebco,lat_gebco,dshore_gebco',  lin_info.lon(:,:,ia),  lin_info.lat,'linear');
  end
  clear ia xa xl;
  
  % Debug plot
  if(do_debug)
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1);
      xx = track_info.lon(:,:,1);    xx = xx(:);
      yy = track_info.lat(:);
      mm = track_info.mask(:,:,1);
      ii = find(mm(:)==1);
      zz = dshore_track(:,:,1); zz = zz(:);
      scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); caxis([0 2000]); title('dshore track');
    subplot(2,1,2);
      xx = lin_info.lon(:,:,1);    xx = xx(:);
      yy = lin_info.lat(:);
      mm = lin_info.mask(:,:,1);
      ii = find(mm(:)==1);
      zz = dshore_lin(:,:,1); zz = zz(:);
      scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); caxis([0 2000]); title('dshore lin');
    cb=colorbar; set(get(cb,'ylabel'),'string','distance from shore (km)');
    print(gcf,'-dpng','Plots/Debug/Debug_Step01-04_AltimPointDistFromShore.png');
    close all; clear xx yy zz mm ii;
  end
      
  % Save
  dshore = dshore_track;    save('Data/track_info_raw.mat','dshore','-append');     clear dshore dshore_track;
  dshore = dshore_lin;      save('Data/lin_info_raw.mat',  'dshore','-append');     clear dshore dshore_lin;
  
  % Clean-up
  clear *_info *_gebco
  
end
%==========================================================================
if(Switch(05))
    
  % Get basin data
  basin_info  = importdata(file_basins);
  basin_lat   = basin_info.data(:,1);
  basin_lon   = basin_info.data(:,2);   basin_lon(basin_lon<0) = basin_lon(basin_lon<0)+360;
  basin_index = basin_info.data(:,3);
  clear basin_info basin_file;
  
  % Create interpolant
  ntrplnt = scatteredInterpolant(basin_lon,basin_lat,basin_index,'nearest');
  
  % Loop through track info's
  str_info = {'track','lin'};
  for ii=1:numel(str_info)
      
    % Load track info
    load(['Data/' str_info{ii} '_info_raw.mat'],'lon','lat');
    
    % Interpolate
    basin = NaN(size(lon));
    for jj=1:2
      basin(:,:,jj) = ntrplnt(lon(:,:,jj),lat);
    end
    
    % Save
    save(['Data/' str_info{ii} '_info_raw.mat'],'basin','-append');
    
    % Clean-up
    clear lon lat basin;
    
  end
  clear ii jj basin_* ntrplnt str_info
    
end
%==========================================================================
if(Switch(06))
   
  % Load info's  
  info_track = load('Data/track_info_raw.mat');
  info_lin   = load('Data/lin_info_raw.mat');
    
  % Save via mission
  for ia=1:n_altm
      
    % Track info
       lon = info_track.lon(:,:,ia);
       lat = info_track.lat;
      dist = info_track.dist(:,:,ia);
      bath = info_track.bath(:,:,ia);
      mask = info_track.mask(:,:,ia);
    dshore = info_track.dshore(:,:,ia);
    basin  = info_track.basin(:,:,ia);
      
    % Save
    save(['Data/track_info_' str_altm{ia} '.mat'],'lon','lat','dist','bath','mask','dshore','basin');
    clear lon lat dist bath mask dshore basin;
    
    % Lin info
       lon = info_lin.lon(:,:,ia);
       lat = info_lin.lat;
      dist = info_lin.dist;
      bath = info_lin.bath(:,:,ia);
      mask = info_lin.mask(:,:,ia);
    dshore = info_lin.dshore(:,:,ia);
    basin  = info_lin.basin(:,:,ia);
      
    % Save
    save(['Data/lin_info_' str_altm{ia} '.mat'],'lon','lat','dist','bath','mask','dshore','basin');
    clear lon lat dist bath mask dshore basin; 
      
  end
  clear ia;
  
  % Delete original files
  %delete('Data/track_info.mat');
  %delete('Data/lin_info.mat');
      
end
%==========================================================================
if(Switch(07))
    
  % Loop through altimeter missions
  for ia=1:n_altm
  
    % Load track info
    track_info = load(['Data\track_info_' str_altm{ia} '.mat'],'lon','lat');
    n_track    = size(track_info.lon,1);
    
    % Variable to determine: index_link
    index_link = cell(n_track,1);
    for it=1:n_track
      
      % Model file
      fmdl = [dir_mdl '\' str_mdl{ia} '\_Processed\Track' sprintf('%0.3d',it) '.nc']; 
        
      % Model points
      xm = ncread(fmdl,'x');    xm(xm>360) = xm(xm>360)-360;
      ym = ncread(fmdl,'y');
      
      % Find track point closest to each model point
      xt = track_info.lon(it,:);
      yt = track_info.lat(it,:);
      
      % Loop x if need be
      xm = unwrap(xm.*(pi/180)).*(180/pi);
      xt = unwrap(xt.*(pi/180)).*(180/pi);
      
      % Ensure x's unwrapped around the same boundary
      if(min(xm)>max(xt));      xm = xm-360;    end
      if(max(xm)<min(xt));      xm = xm+360;    end
       
      % Positions
      pm = xm + sqrt(-1).*ym;
      pt = xt + sqrt(-1).*yt;
      
      % Loop though model positions and link to track positions
      np = numel(pm);
      index_link{it} = NaN(np,1);
      for ip=1:np
          
        % Distances
        dd = abs(pm(ip)-pt);
        
        % Min distance
        index_link{it}(ip) = find(dd==min(dd),1);
        
        % Clean-up
        clear dd;
          
      end
      clear ip;
        
      % Clean-up
      clear xm ym xt yt pm pt np;
        
    end
    clear it;
    
    % Debug plot
    if(do_debug)
      figure;
        x = track_info.lon + sqrt(-1).*track_info.lat;     
        y = x;
        for it=1:n_track
          y(it,~ismember(1:size(y,2),index_link{it})) = NaN;
        end
        plot(real(x(:)),imag(x(:)),'.k',real(y(:)),imag(y(:)),'.b');
      axis tight; box on; title('index link');
      print(gcf,'-dpng',['Debug_step01-6_index-link_' str_altm{ia} '.png']);
      clear x y it; close all;
    end
    
    % Save (appends)
    save(['Data\track_info_' str_altm{ia} '.mat'],'index_link','-append');
    
  end
  clear ia;
    
end
%==========================================================================
