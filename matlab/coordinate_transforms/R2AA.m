function AA = R2AA(R)
% function AA = R2AA(R)
% Convert Rotation Matrix (R) to Axis Angle (AA)
%
% Input
%     R: a 3x3 rotation matrix
%
% Ouput:
%     AA: a 1x4 axis angle. The axis is not normalized
%
% Written by Garrick Orchard July 2017
% Based on:
% https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis.E2.80.93angle

% Check that input is a rotation matrix
assert(all(size(R)==3) && numel(R)==9,'Input must be a 3x3 matrix');
assert(norm(R.'*R - eye(3),'fro') < 1e-7,'Input must be a 3D rotation matrix');
assert(det(R) > 0,'Input must be a 3D rotation matrix');

% Get rotation angle
theta = acos((trace(R)-1)/2);
theta = real( theta ); % in case cosine is slightly out of the range [-1,1]

% Get rotation axis
if (abs(theta-pi) < 1e-3)
    % Rotations with an angle close to pi
    
    % Obtain the axis from the quadratic term in [u]_x in Rodrigues formula
    % Get the vector that generates the rank-1 matrix
    [U,~,~] = svd( 0.5*(R + eye(3)) );
    ax = U(:,1).';
    
    % Adjust the sign of the axis
    if (norm(AA2R([ax, theta])-R,'fro') > norm(AA2R([-ax, theta])-R,'fro'))
        ax = - ax;
    end
    
else
    % Most rotations
    
    % Obtain the axis from the linear term in [u]_x in Rodrigues formula
    ax = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)];
    
    norm_ax = norm(ax);
    if (norm_ax > 1e-8)
        ax = ax./norm_ax; % Ensure a unit length axis direction
    else
        % Rotation close to zero degrees. Axis is undetermined
        ax = [0 0 1]; % this is what rotm2axang outputs
    end
end

% Ouput 4-vector: [axis, angle]
AA = [ax, theta];