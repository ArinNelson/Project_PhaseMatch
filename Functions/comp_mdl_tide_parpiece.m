function [tidefit, tidets] = comp_mdl_tide_parpiece(time,dat_filt,lat)

  % Constants
  str_tide = {'M2','S2','N2','K1','O1'};
  n_tide   = numel(str_tide);

  % Tide fitting
  coef = ut_solv(time(:),dat_filt(:),[],lat,str_tide,...
                 'OLS','LinCI','NoDiagn','RunTimeDisp','nnn','NodsatLinT','GwchLinT','NoTrend'); 
             
  % Save results
  tidefit      = NaN(n_tide,2);
  for it=1:n_tide
    ii = find(strcmp(coef.name,str_tide{it})==1);
    tidefit(it,1) = coef.A(ii) .* exp(sqrt(-1).*coef.g(ii).*(pi/180));
    tidefit(it,2) = coef.A_ci(ii); %.* ( squeeze(tidefit(:,1))./coef.A )) + (sqrt(-1)*coef.g_ci.*(pi/180));
  end
  clear it 

  % Compute stationary tide timeseries
  tidets = NaN(numel(time),n_tide);
  for iw=1:n_tide  
  if(abs(tidefit(iw,1))>tidefit(iw,2)) 
    tidets(:,iw) = ut_reconstr(time,coef,'Cnstit',str_tide{iw});
  end
  end
  clear iw;
  
%   % Test: plot
%   var_filt = nanvar(dat_filt);
%   var_tide = nanvar(nansum(tidets,2));
%   anom     = dat_filt - nansum(tidets,2);
%   var_anom = nanvar(anom);
%   test = [var_filt var_tide var_anom];
%   plot(time,dat_filt,'-k',time,nansum(tidets,2),'-b',time,dat_filt-nansum(tidets,2),'-r');
%   title(num2str(test));
%   
%   % Test: ensure tide time series removes variance, but not (much) more 
%   % variance than the full (stationary + nonstationary) time series
%   if(var_anom > 0.25*var_filt)
%     pause(1e-9);
%   end

end