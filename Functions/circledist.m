function D = circledist(o1,p1,o2,p2)
% function D = circledist(o1,p1,o2,p2)
% computes the normalized distance on a sphere between lon/lat points o1,p1
% and o2,p2 utilizing the haversine formula to compute the central angle
% between the points

    % cartesian points
    x1 = sind(90-p1)*cosd(o1);        x2 = sind(90-p2)*cosd(o2);
    y1 = sind(90-p1)*sind(o1);        y2 = sind(90-p2)*sind(o2);
    z1 = cosd(90-p1);                 z2 = cosd(90-p2);

    % cartesian distance
    X = [x1,y1,z1];
    Y = [x2,y2,z2];
    D = atan2(norm(cross(X,Y)),dot(X,Y));

end
