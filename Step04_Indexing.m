% Determine segment and box points
% Arin Nelson
%==========================================================================

% Switches
Switch     = zeros(9,1);
Switch(01) = 0;     % Determine valid segments
Switch(02) = 1;     % Determine bin indices
Switch(03) = 1;     % Determine box indices

% Directories
dir_data   = '.\Data';
dir_mdl    = '\Datasets\HYCOM\PhaseMatch';
str_mdl    = {'NonDA_2016','DA_2016','DA_Modern'};
str_altm   = {'j2','j3'};
str_mthd   = {'track','lin'};

% Options
do_debug   = 1;  
ref_m_to_a = [1 1 2];
ref_a_to_m = [2 3];

% Options for segments
num_min    = [7,111];
dshore_min = 12;
bath_min   = 500;
gap_max    = 30;
length_min = 1200;
basin_incl = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21];

% Options for bins
dlat   = 5;
dlon   = 5;
bin_lat_edge = -65 : dlat : 65;
bin_lat_cntr = (bin_lat_edge(2:end)+bin_lat_edge(1:end-1))./2;
bin_lon_edge = 0 : dlon : 360;
bin_lon_cntr = (bin_lon_edge(2:end)+bin_lon_edge(1:end-1))./2;

% Options for boxes (lon range, lat range, basin indices, name/label)
box_info = {[140  240],     [ 40  60],  [8 9 10 11 12 13], 'N  PACIFIC';   ...
            [120  170],     [ 05  40],  [8 9 10 11 12 13], 'LUZON';        ...
            [170  210],     [ 05  40],	[8 9 10 11 12 13], 'HAWAII';       ...
            [210  280],     [ 05  40],  [8 9 10 11 12 13], 'NW PACIFIC';   ...
            [160  280],     [-05  05],  [8 9 10 11 12 13], 'EQ PACIFIC';  ...
            [150  180],     [-30 -05],  [8 9 10 11 12 13], 'TAHITI';       ...
            [180  230],     [-30 -05],	[8 9 10 11 12 13], 'SC PACIFIC';	...
            [230  290],     [-30 -05],  [8 9 10 11 12 13], 'SW PACIFIC';   ...
            [150  290],     [-60 -30],  [8 9 10 11 12 13], 'S PACIFIC';    ...
            [ 40  100],     [-05  30],  [14 15 16 17 18 19], 'N INDIAN';     ...
            [ 30   70],     [-30 -05],  [14 15 16 17 18 19], 'MADAGASCAR';   ...
            [ 70  130],     [-30 -05],  [14 15 16 17 18 19], 'W INDIAN';     ...
            [ 20  150],     [-60 -30],  [14 15 16 17 18 19], 'S INDIAN';     ...
            [300  360],     [ 50  66],  [2 3 4 5 6 7], 'N ATLANTIC';   ...
            [280  320],     [ 20  50],  [2 3 4 5 6 7], 'GULF STREAM';  ...
            [320  360],     [ 20  50],  [2 3 4 5 6 7], 'NW ATLANTIC';  ...
            [280  380],     [-20  20],  [2 3 4 5 6 7], 'EQ ATLANTIC';  ...
            [290  380],     [-60 -20],  [2 3 4 5 6 7], 'S ATLANTIC';   ...
           };   

% Constants
n_track = 254;                      % # altimeter tracks
n_mdl   = numel(str_mdl);           % # model simulations
n_altm  = numel(str_altm);          % # altimeter missions
n_mthd  = numel(str_mthd);
n_binx  = numel(bin_lon_cntr);
n_biny  = numel(bin_lat_cntr);
n_box   = size(box_info,1);


%==========================================================================
if(Switch(01))
    
  % Loop through simulations
  for ia=1:n_altm
  for id=1:n_mthd
      
    % Load track info
    info = load([dir_data '\' str_mthd{id} '_info_' str_altm{ia} '.mat']);
    n_track    = size(info.lon,1); 
    n_point    = size(info.lon,2);
    
    % Load model info
    mdl_info = load([dir_data '\mdl_info_' str_mdl{ref_a_to_m(ia)} '.mat'],[str_mthd{id} '_is_avail']);
    
    % Load altimeter data
    altm_data = load([dir_data '\' str_mthd{id} '_data_' str_altm{ia} '.mat']);
    
    % New variables
    is_valid = NaN(n_track,n_point);
    segnum   = NaN(n_track,n_point);
    
    % Loop through tracks
    for it=1:n_track
    clc; disp(['Altim mission ' num2str(ia) ', sampling method ' num2str(id) ', Track ' num2str(it) '...']);    
    %try
        
      % Data for this track
%       f_mdl = [dir_mdl '\' str_mdl{ref_a_to_m(ia)} '\_Processed\MATLAB\Track' sprintf('%0.3d',it) 'mat.mat'];
%       mdl_data = load(f_mdl,[str_mthd{id} '_dat']);
%       eval(['mdl_data = mdl_data.' str_mthd{id } '_dat;']);
      eval(['mdl_data = squeeze(mdl_info.' str_mthd{id} '_is_avail(it,:,:))'';']);
      
      % Info for this track
      altm_test = sum(~isnan(altm_data.sla{it}) & altm_data.sla{it}~=0);
      mdl_test  = sum(~isnan(mdl_data) & mdl_data~=0);
      
      % Valid points
      is_valid(it,:) = (  altm_test>=num_min(ia) ...
                        & info.dshore(it,:)>=dshore_min ...
                        & info.bath(it,:)>=bath_min ...
                        & info.mask(it,:)==1 ...
                        & mdl_test>=num_min(ia) ...
                        & ismember(info.basin(it,:),basin_incl) ...
                       );
      
      % Separate data into segments
      i1 = find( ~is_valid(it,1:end-1) &  is_valid(it,2:end) );
      i2 = find(  is_valid(it,1:end-1) & ~is_valid(it,2:end) );

      % Some fixes
      if(isempty(i1) && isempty(i2))
            
        % If only one segment, segment = full track!
        i1 = 0;
        i2 = n_point;
          
      else
            
        % Ensure end points are counted in segments  
        if(is_valid(it,1));       i1 = [0; i1(:)];            end
        if(is_valid(it,end));     i2 = [i2(:); n_point];      end
          
        % Remove last i1 if no more non-NaN points afterwards
        if(i1(end)>i2(end));	i1(end) = [];               end
          
      end
      
%       % Check
%       plot(1:n_point,is_valid,'-k',i1+1,is_valid(i1+1),'or',i2,is_valid(i2),'db');

%       % Length of each segment
%       try          
%         seg_lngth = info.dist(it,i2)-info.dist(it,i1+1); 
%         seg_gap  = info.dist(it,i1(2:end)+1) - info.dist(it,i2(1:end-1));
%       catch err 
%         seg_lngth = info.dist(i2)-info.dist(i1+1);
%         seg_gap  = info.dist(i1(2:end)+1) - info.dist(i2(1:end-1));
%       end
      
      % Determine size of gaps.  If gap <length_gap, combine the tracks.
      try d=info.dist(it,:); catch err; d=info.dist(:); end
      ns_old = numel(i1);
      ns_new = 0;
      while(ns_new~=ns_old)
        ii = 1;
        while ii<numel(i1)
          
          % Perform gap check
          if( d(i1(ii+1)+1)-d((i2(ii))) < gap_max )
            i2(ii)   = [];
            i1(ii+1) = [];
          end
          ii = ii + 1;
          
        end
        ns_old = ns_new;
        ns_new = numel(i1); 
      end
      clear ii ns_old ns_new;
      
      % Remove segments with length < length_min
      min_seg = min(d(i2)-d(i1+1));
      while(min_seg < length_min)
            
        % Remove short segment(s)  
        ii = 1;
        while ii<=numel(i1)
            
          % Perform length check
          if( d(i2(ii))-d(i1(ii)+1) < length_min )
            i2(ii) = [];
            i1(ii) = [];
          end
          ii = ii + 1;
            
        end
        clear ii;
          
        % Check again
        min_seg = min(d(i2)-d(i1+1));
          
      end
      clear min_seg
      
      % Save the segments
      for ii=1:numel(i1)
        segnum(it,(i1(ii)+1):i2(ii)) = ii;
      end
      clear ii i1 i2;
      
      % Clean-up
      clear f_mdl mdl_data;
     
    %catch err; end  
    end
    clear it;
    
    % Debug plot
    if(do_debug)
      figure;
        xx = info.lon(~isnan(segnum));
        yy = info.lat(~isnan(segnum));
        zz = segnum(~isnan(segnum));
        scatter(xx,yy,2,zz,'o','filled'); colorbar;
        title(['Segment Numbering (' str_altm{ia} ')']);
      print(gcf,'-dpng',['Plots/Debug/Debug_Step03-02_' str_mthd{id} '_' str_altm{ia} '_segmenting.png']);
      close all; clear xx yy zz;
    end

    % Save
    save([dir_data '\' str_mthd{id} '_info_' str_altm{ia} '.mat'],'segnum','is_valid','-append');
    
    % Clean-up
    clear info n_track n_point altm_data segnum is_valid
    
  end
  end
  clear ia id;
    
end
%==========================================================================
if(Switch(02))
    
  % Loop through info's and methods
  for ia=1:n_altm
  for id=1:n_mthd
   
    % Load info  
    info = load([dir_data '\' str_mthd{id} '_info_' str_altm{ia} '.mat'],'lon','lat','segnum','is_valid');
    
    % The bins
    bin_index = NaN(size(info.lon));
    for ix=1:n_binx
    for iy=1:n_biny
    
      % Shouldn't be too bad
      ii = find(  info.lon  > bin_lon_edge(ix)   ...
                & info.lon <= bin_lon_edge(ix+1) ...
                & info.lat  > bin_lat_edge(iy)   ...
                & info.lat <= bin_lat_edge(iy+1) ...
                & info.is_valid==1 );
            
      % Save
      bin_index(ii) = ix + sqrt(-1).*iy;
            
    end
    end
    clear ix iy ii;
    
    % Debug plot
    if(do_debug)
      figure;
        subplot(2,1,1);
          scatter(info.lon(:),info.lat(:),2,real(bin_index(:)),'o','filled'); 
          axis tight; colorbar; colormap(jet); title('bin index (x)');
        subplot(2,1,1);
          scatter(info.lon(:),info.lat(:),2,real(bin_index(:)),'o','filled');
          axis tight; colorbar; colormap(jet); title('bin index (y)');
      print(gcf,'-dpng',['Plots\Debug\Debug_Step04-2_bin-index_' str_mthd{id} '_' str_altm{ia} '.png']);
      close all;
    end
        
    % Save
    save([dir_data '\' str_mthd{id} '_info_' str_altm{ia} '.mat'],'bin_index','-append');
    
  end
  end
  clear ia id;
    
end
%==========================================================================
if(Switch(03))
    
  % Loop through info's and methods
  for ia=1:n_altm
  for id=1:n_mthd
   
    % Load info  
    info = load([dir_data '\' str_mthd{id} '_info_' str_altm{ia} '.mat'],'lon','lat','basin','segnum','is_valid');
    
    % The bins
    box_index = NaN(size(info.lon));
    for ib=1:n_box
    
      % Shouldn't be too bad
      ii = find(  info.lon  > box_info{ib,1}(1) ...
                & info.lon <= box_info{ib,1}(2) ...
                & info.lat  > box_info{ib,2}(1) ...
                & info.lat <= box_info{ib,2}(2) ...
                & ismember(info.basin, box_info{ib,3}) ...
                & info.is_valid == 1 );
            
      % Save
      box_index(ii) = ib;
            
    end
    clear ii ib;
    
    % Debug plot
    if(do_debug)
      figure;
        scatter(info.lon(:),info.lat(:),2,box_index(:),'o','filled'); 
        axis tight; colorbar; colormap(hsv(n_box)); title('box index');
      print(gcf,'-dpng',['Plots\Debug\Debug_Step04-3_box-index_' str_mthd{id} '_' str_altm{ia} '.png']);
      close all;
    end
    
    % Save
    save([dir_data '\' str_mthd{id} '_info_' str_altm{ia} '.mat'],'box_index','-append');
    
  end
  end
  clear ia id;
    
end
