function [h,s,p,N,T,time_mjd]=ray_arguments(y,d);

% NOTE ADDED APRIL 17, 2019:  THIS CODE DOES NOT ACTUALLY USE THE "DAY" ARGUMENT; SO IT IS WRITTEN FOR EVERYTHING TO REFERENCED TO "OOZ" ON JANUARY 1 OF THE YEAR IN QUESTION


% y = year number (i.e. 2003)
% d = day number (1-366)

inty=floor((y-1857)/4)-1;

time_mjd=365*(y-1858)+inty-(31+28+31+30+31+30+31+31+30+31+17) + 1;

T=time_mjd-51544.4993;

s=218.3164+13.17639648*T;

h=280.4661+0.98564736*T;

p=83.3535+0.11140353*T;

N=125.0445 - 0.05295377*T;









