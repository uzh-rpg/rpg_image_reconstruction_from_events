function R = q2R( q )
% function R = q2R( q )
% Convert quaternion (q) to 3x3 rotation matrix (R)
%
% Input:
%     q: [qw, qx, qy, qz]*[1 i j k]';
%
% Ouput:
%     R: a 3x3 Rotation Matrix
%
% Written by Garrick Orchard July 2017
% Based on:
% http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/

R = zeros(3);

R(1,1) = 1 - 2*(q(3)^2 + q(4)^2);
R(1,2) = 2*(q(2)*q(3) - q(4)*q(1));
R(1,3) = 2*(q(2)*q(4) + q(3)*q(1));

R(2,1) = 2*(q(2)*q(3)+q(4)*q(1));
R(2,2) = 1 - 2*(q(2)^2 + q(4)^2);
R(2,3) = 2*(q(3)*q(4) - q(2)*q(1));

R(3,1) = 2*(q(2)*q(4) - q(3)*q(1));
R(3,2) = 2*(q(3)*q(4) + q(2)*q(1));
R(3,3) = 1 - 2*(q(2)^2 + q(3)^2);

% Numerically improve result by projecting on the space of rotation matrices
[u,~,v] = svd(R);
R = u * diag([1,1,det(u*v.')]) * v.'; % Closest rot matrix