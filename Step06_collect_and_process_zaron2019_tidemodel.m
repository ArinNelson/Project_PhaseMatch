% Load and process Ed's baroclinic tide model
%==========================================================================

% Switches
Switch     = zeros(9,1);
Switch(01) = 1;     % Load and process raw data
Switch(02) = 1;     % Interpolate to method points

% Strings
dir_zaron  = '\Datasets\Altimetry\Zaron2019\';
file_zaron = 'HRET_v8.1_compressed.nc';
dir_data = '.\Data';
str_altm = {'j2','j3'};
str_mthd = {'track','lin'};

% Options
do_debug   = 1;                 % Show debug plots and information

% Constants
n_track = 254;                  % # altimeter tracks
n_altm  = numel(str_altm);      % # altimeter missions
n_mthd  = numel(str_mthd);

%==========================================================================
if(Switch(01))

  % Data location  
  file_ed   = [dir_zaron file_zaron];
  str_ed    = {'M2','S2','K1','O1'};    
  ne        = numel(str_ed);
  
  % Load data
  lon = ncread(file_ed,'longitude');    nx = numel(lon);
  lat = ncread(file_ed,'latitude');     ny = numel(lat);
  tide_ed = NaN(nx,ny,ne);
  for i=1:ne
    tide_ed(:,:,i) = ncread(file_ed,[str_ed{i} 're']) + sqrt(-1).*ncread(file_ed,[str_ed{i} 'im']);
  end
  clear i;
  tide_ed(tide_ed==0) = NaN;
  
  % Ed's data is in units of m, convert to cm
  tide_ed = tide_ed.*100;
  
  % Save raw data?
  tide_re = real(tide_ed);
  tide_im = imag(tide_ed);
  tide_name = str_ed;
  save([dir_data '/Zaron2019.mat'],'lon','lat','tide_re','tide_im','tide_name');
  
  % Clean-up
  clear file_data file_ed str_ed ne lon lat nx ny tide_*;
  
end
%==========================================================================
if(Switch(02))
    
  % Load Ed's data
  zaron = load([dir_data '\Zaron2019.mat']);
  
  % Prepare interpolant
  [xx,yy] = meshgrid(zaron.lon,zaron.lat);  %xx=xx'; yy=yy';
  zr = zaron.tide_re;
  zi = zaron.tide_im;
    
  % Loop through methods and altimeter missions
  for ia=2 %1:n_altm
  for id=1:n_mthd
  clc; disp(['Performing interpolation for mission #' num2str(ia) ', sampling method #' num2str(id) '...']);
  
    % Info
    info = load([dir_data '\' str_mthd{id} '_info_' str_altm{ia} '.mat']);
    
    % Valid interp points
    ii = find(~isnan(info.segnum));
    
    % Tide info
    tide_zaron_re = NaN(size(info.lon,1),size(info.lon,2),4);
    tide_zaron_im = NaN(size(info.lon,1),size(info.lon,2),4);
    
    % Loop through tides
    for iw=1:4
        
      % Real piece  
      tmp = tide_zaron_re(:,:,iw);
      tmp(ii) = interp2(xx,yy,zr(:,:,iw)',info.lon(ii),info.lat(ii),'linear');
      tide_zaron_re(:,:,iw) = tmp;
      
      % Imag piece
      tmp = tide_zaron_im(:,:,iw);
      tmp(ii) = interp2(xx,yy,zi(:,:,iw)',info.lon(ii),info.lat(ii),'linear');
      tide_zaron_im(:,:,iw) = tmp;
      
    end
    clear iw tmp;
    
    % Debug plot
    if(do_debug)
      figure('units','normalized','outerposition',[0 0.25 1 0.5]);
        for iw=1:4
          zA = log10( sqrt( tide_zaron_re(:,:,iw).^2 + tide_zaron_im(:,:,iw).^2 ) );
          zP = cos(atan2(tide_zaron_im(:,:,iw),tide_zaron_re(:,:,iw)));  
          subplot(2,4,iw);
            scatter(info.lon(ii), info.lat(ii), 2, zA(ii), 'o', 'filled'); 
            axis tight; box on; set(gca,'color','k');
            colormap(jet); colorbar; caxis([-3 1]);
            title(['Zaron 2019 ' zaron.tide_name{iw} ' Amplitude ( Log_{10}(cm) )']);
          subplot(2,4,iw+4)
            scatter(info.lon(ii), info.lat(ii), 2, zP(ii), 'o', 'filled'); 
            axis tight; box on; set(gca,'color','k');
            colormap(jet); colorbar; caxis([-1 1]);
            title(['Zaron 2019 ' zaron.tide_name{iw} ' cos(Phase)']);
        end
      print(gcf,'-dpng',['Plots\Debug\Debug_Step06-2_' str_mthd{id} '_' str_altm{ia} '.png']);
      close all; clear x iq;
    end
    
    % Save
    save([dir_data '\' str_mthd{id} '_info_' str_altm{ia} '.mat'],'tide_zaron_re','tide_zaron_im','-append');
    
    %Clean-up
    clear info ii tide_zaron;
    
  end
  end
  clear ia id;
    
end
