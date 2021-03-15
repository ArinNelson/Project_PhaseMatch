function D = circledist(o1,p1,o2,p2)
% function D = circledist(o1,p1,o2,p2)
% computes the normalized distance on a sphere between lon/lat points o1,p1
% and o2,p2 utilizing the haversine formula to compute the central angle
% btween the points

% Some strange error occurs sometimes... for now let's try x,y
%D = sqrt( (o2-o1)^2 + (p2-p1)^2 );

% % radians?
% o1 = o1*(pi/180);
% o2 = o2*(pi/180);
% p1 = p1*(pi/180);
% p2 = p2*(pi/180);

% cartesian?
x1 = sind(90-p1)*cosd(o1);        x2 = sind(90-p2)*cosd(o2);
y1 = sind(90-p1)*sind(o1);        y2 = sind(90-p2)*sind(o2);
z1 = cosd(90-p1);                 z2 = cosd(90-p2);

% cartesian distance?
X = [x1,y1,z1];
Y = [x2,y2,z2];
D = atan2(norm(cross(X,Y)),dot(X,Y));

% % angle differences
% do = o2-o1;
% dp = p2-p1;
% 
% % various arguments
% cp1 = cos(p1);
% sp1 = sin(p1);
% cp2 = cos(p2);
% sp2 = sin(p2);
% cdo = cos(do);
% sdo = sin(do);
% 
% % numerator
% a = cp2*sdo;
% b = cp1*sp2 - sp1*sp2*cdo;
% n = sqrt( a^2 + b^2 );
% 
% % denominator
% d = sp1*sp2 + cp1*cp2*cdo;
% 
% % the distance
% D =atan2(n,d);

end