function y = nanfill(x)
% Fill in NaNs using linear interpolation

t = 1:numel(x);
y = x;
y(isnan(x)) = interp1(t(~isnan(x)),x(~isnan(x)),t(isnan(x)));

end