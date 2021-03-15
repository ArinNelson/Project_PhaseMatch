function y = nandetrend(x,n)
  x = x(:);
  t = 1:numel(x);   t = t(:);
  i = find(~isnan(x));
  if(~isempty(i))
    tt = t(i);        tt = tt(:);
    xx = x(i);        xx = xx(:);
    p = polyfit(tt,xx,n);
    z = polyval(p,t); z = z(:);
    y = x (:) - z(:);
  end
end