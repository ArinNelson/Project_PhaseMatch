function yp = naninterp1(x,y,xp,mthd)
% yp = naninterp1(x,y,xp,mthd)
% Works like interp1 but works with NaNs in y (but NOT in x!).

  % The interpolated values
  yp = NaN(numel(xp),1);
  
  % Set NaN values to some huge number
  yy = y;
  ymax = max(y(~isnan(y)));
  ymin = min(y(~isnan(y)));
  if(any(isnan(yy)))
    yy(isnan(yy)) = ymax*1e30;
  end
  
  % Do interpolation
  yp = interp1(x(~isnan(x)),yy(~isnan(x)),xp,mthd);

  % Large values are NaNs
  yp(yp>ymax) = NaN;
  yp(yp<ymin) = NaN;
  
end