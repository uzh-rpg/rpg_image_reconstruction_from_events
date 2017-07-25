function R = AA2R(AA)
% function R = AA2R(AA)
% Convert Axis Angle (AA) to Rotation Matrix (R)
%
%  Input
%     AA: 1x4 axis angle rotation.
%
% Ouput:
%     R: a 3x3 rotation matrix
%
% Written by Garrick Orchard July 2017
% Based on:
% http://mathworld.wolfram.com/RodriguesRotationFormula.html

assert(numel(AA)==4,'Input must be 1x4 or 4x1');

% Axis
norm_ax = norm(AA(1:3));
if (norm_ax < 1e-6)
    R = eye(3);
    return
end
ax = AA(1:3) ./ norm_ax; % Unit norm, avoid division by zero

% Cross-product matrix
omega = [   0, -ax(3),  ax(2);
        ax(3),      0, -ax(1);
       -ax(2),  ax(1),  0];

% Rotation angle
theta = AA(4);

% Rotation matrix, using Rodrigues formula
R = eye(3) + omega*sin(theta) + omega*omega*(1-cos(theta));