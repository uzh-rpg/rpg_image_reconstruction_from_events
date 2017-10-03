function [point_2d, d2d3] = project_EquirectangularProjection(point_3d, pano_width, pano_height)
% project_EquirectangularProjection Project a 3D point according to equirectangular model
% Used for 360 degrees panoramic cameras that output a full panorama frame
%
% Input:
% -point_3d: a 3D point
% -pano_width: width of the panorama
% -pano_height: height of the panorama
%
% Output:
% -point_2d: projected point (coordinates in the panorama image)

rho = sqrt(sum(point_3d.*point_3d,1)); % norm of each 3D point

fx = pano_width / (2.0 * pi);
fy = pano_height / pi;
principal_point = 0.5 * [pano_width; pano_height];

phi = atan2(point_3d(1,:),point_3d(3,:));
theta = asin(-point_3d(2,:) ./ rho);
point_2d = [phi * fx; -theta * fy];
point_2d(1,:) = point_2d(1,:) + principal_point(1);
point_2d(2,:) = point_2d(2,:) + principal_point(2);

if nargout > 1
    num_points = size(point_3d,2);
    d2d3 = zeros(2,3,num_points);
    X = point_3d(1,:);
    Y = point_3d(2,:);
    Z = point_3d(3,:);
    tmp1 = fx ./((1 + (X./Z).^2) .* Z);
    tmp2 = -fy ./sqrt(1 - (Y./rho).^2);
    tmp3 = Y./(rho.^3);
    d2d3(1,1,:) = tmp1;
    d2d3(1,3,:) = -tmp1 .* (X./Z);
    d2d3(2,1,:) = tmp2 .* tmp3 .* X;
    d2d3(2,2,:) = tmp2 .* (tmp3.*Y - 1./rho);
    d2d3(2,3,:) = tmp2 .* tmp3 .* Z;
end
