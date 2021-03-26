function [v_z_t,   n_z_t,   ...
          v_m_t,   n_m_t,   ...
          v_m_f,   n_m_f,   ...
          s2,      d2,    ...
          v_a,     n_a,     ...
          v_n_z_t, n_n_z_t, ...
          v_n_m_t, n_n_m_t, ...
          v_n_m_f, n_n_m_f ...
         ] = var_reduc_tests_parpiece(tide_zaron_re, tide_zaron_im, ...
                                      altm_time, tref, ...
                                      mdl_time, dat_tide_timeseries, ...
                                      dat_filt, mdl_tide, sla ...
                                     )
                                 
    %----------------------------------------------------------------------
    % CONSTANTS
    
    % Much more convienient...
    [omega,chi,f,nu] = tideinfo(tref);
    
    %----------------------------------------------------------------------
    % INIT
    
    % De-mean altimeter data
    if(sum(~isnan(sla))>1)
      sla = nandetrend(sla,0);
    end                           
                                 
	% Construct Zaron 2019 baroclinic tide time series at this point
    ii = [1 2 4 5];
	zaron_tide_timeseries = NaN(numel(altm_time),8);
	for iw=1:4
	if(~isnan(tide_zaron_re(iw)))                          
                                 
      % Compute the stationary tide time series
      % VERIFIED TO BE RIGHT! SEE OLDER CODES...
      trigarg = omega(ii(iw)).*(altm_time-tref) + (pi/180)*(chi(ii(iw))+nu(ii(iw)));
      R = tide_zaron_re(iw);
      I = tide_zaron_im(iw);
      A = sqrt(R^2 + I^2);
      P = atan2(I,R);
      zaron_tide_timeseries(:,ii(iw)) = f(ii(iw))*A*cos(trigarg - P);
          
	end
	end
    clear iw ii;                         
                                 

	% Interpolate model time series to altimeter times
	mdl_tide_timeseries = NaN(numel(altm_time),8);
	for iw=1:5
    if(~all(isnan(dat_tide_timeseries(:,iw))))    
        mdl_tide_timeseries(:,iw) = naninterp1(mdl_time,dat_tide_timeseries(:,iw),altm_time(:),'linear');
    end
	end
    clear iw;
    
	% Same for filtered time series
	mdl_filt_timeseries = NaN(numel(altm_time),4);
	for iw=1:4
	if(~all(isnan(dat_filt(:,iw))))    
        mdl_filt_timeseries(:,iw) = naninterp1(mdl_time,dat_filt(:,iw),altm_time(:),'linear');
	end
	end
    clear iw;
    
	% Some totals
	zaron_tide_timeseries(:,6) = nansum(zaron_tide_timeseries(:,1:3),2); 
	zaron_tide_timeseries(:,7) = nansum(zaron_tide_timeseries(:,4:5),2); 
	zaron_tide_timeseries(:,8) = nansum(zaron_tide_timeseries(:,1:5),2); 
	  mdl_tide_timeseries(:,6) = nansum(  mdl_tide_timeseries(:,1:3),2); 
	  mdl_tide_timeseries(:,7) = nansum(  mdl_tide_timeseries(:,4:5),2); 
	  mdl_tide_timeseries(:,8) = nansum(  mdl_tide_timeseries(:,1:5),2); 
      
    % Initial Variances
    v_z_t = nanvar(zaron_tide_timeseries,0,1);    n_z_t = max( sum(~isnan(zaron_tide_timeseries),1) );
    v_m_t = nanvar(mdl_tide_timeseries,  0,1);    n_m_t = max( sum(~isnan(mdl_tide_timeseries)  ,1) );
    v_m_f = nanvar(mdl_filt_timeseries,  0,1);    n_m_f = max( sum(~isnan(mdl_filt_timeseries)  ,1) );
    
    %----------------------------------------------------------------------  
    % Test #0: debug plots
    
    % Plot #1: time series
    %figure; plot(altm_time,sla,'-k',mdl_time,dat_filt(:,3),'-r',mdl_time,sum(dat_tide_timeseries,2),'-b',altm_time,zaron_tide_timeseries(:,end),'-g');
    %plot(altm_time,mdl_tide_timeseries(:,1),'-k',altm_time,zaron_tide_timeseries(:,1),'-b',altm_time,mdl_tide_timeseries(:,1)-zaron_tide_timeseries(:,1),'-r');
    
    %----------------------------------------------------------------------  
    % Test #1: variance in Zaron 2019 explained by HYCOM
    
    % Init
    ii = [1 2 4 5];
    S2 = NaN(7,1);
    D2 = NaN(7,1);
    
    % Amps and phases
    a_m = abs(mdl_tide(ii,1));       a_m(a_m<mdl_tide(ii,2)) = NaN;
    p_m = angle(mdl_tide(ii,1));
    a_z = sqrt( tide_zaron_re(:).^2 + tide_zaron_im(:).^2 );   
    p_z = atan2( tide_zaron_im(:), tide_zaron_re(:) ); 
    
    % S2
    S2(1:4) = 0.5 .* ( a_z.^2 );

    % D2
    D2(1:4) = 0.5 .* (a_m.^2 + a_z.^2) - a_m.*a_z.*cos(p_m - p_z);
    
    % Sums
    S2(5,:) = nansum(S2(1:2,:),1);
    S2(6,:) = nansum(S2(3:4,:),1);
    S2(7,:) = nansum(S2(1:4,:),1);
    D2(5,:) = nansum(D2(1:2,:),1);
    D2(6,:) = nansum(D2(3:4,:),1);
    D2(7,:) = nansum(D2(1:4,:),1);
    
    % Save
    s2 = S2;
    d2 = D2;
    
    %----------------------------------------------------------------------  
    % Test #2: Variance Reductions
    if(~all(isnan(sla)))
    
      % Altimeter variance
      tmp = sla;    tmp = tmp(~isnan(tmp));     v_a = nanvar(tmp);  n_a = numel(tmp);
    
      % Anomaly time series
      anom_zaron_tides = repmat(sla(:),[1 8]) - zaron_tide_timeseries;
      anom_mdl_tides   = repmat(sla(:),[1 8]) - mdl_tide_timeseries;
      anom_mdl_filts   = repmat(sla(:),[1 4]) - mdl_filt_timeseries;
    
      % Variances
      v_n_z_t = nanvar(anom_zaron_tides,0,1);   n_n_z_t = max( sum(~isnan(anom_zaron_tides),1) );
      v_n_m_t = nanvar(anom_mdl_tides,0,1);     n_n_m_t = max( sum(~isnan(anom_mdl_tides),1) );
      v_n_m_f = nanvar(anom_mdl_filts,0,1);     n_n_m_f = max( sum(~isnan(anom_mdl_filts),1) );
      
      % V and N should be enough to use for statistical tests 
      % (assuming means are 0, which they oughta be by construction...)
      
    else
      
      % No altimeter data, no altimeter variance or anomaly variances
      v_a = NaN;            n_a = 0;
      v_n_z_t = NaN(8,1);   n_n_z_t = zeros(8,1);
      v_n_m_t = NaN(8,1);   n_n_m_t = zeros(8,1);
      v_n_m_f = NaN(8,1);   n_n_m_f = zeros(8,1);
      
    end

%--------------------------------------------------------------------------                                 
end