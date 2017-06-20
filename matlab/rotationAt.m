function rot_interp = rotationAt(t_ctrl, rotmats_ctrl, t_query)
% rotationAt Function for linear interpolation of rotations
% Interpolate the orientation (rotation) at a given time.
%
% Input:
% -t_ctrl: timestamps of discrete set of orientations ("control poses")
% -rotmats_ctrl: discrete set of rotation matrices
% -t_query: time of the requested rotation matrix
%
% Output:
% -rot_interp: interpolated rotation matrix
%
% I/O in the style of interp1


if (t_query < t_ctrl(1) || t_query > t_ctrl(end) )
    rot_interp = nan * rotmats_ctrl(1,:);
    
elseif (t_query == t_ctrl(end))
    rot_interp = rotmats_ctrl(end,:);
    
else
    idx_0 = find(t_query > t_ctrl,1,'last');
    % two rotations and their times
    t_0 = t_ctrl(idx_0);
    t_1 = t_ctrl(idx_0+1);
    rot_0 = rotmats_ctrl(:,:,idx_0);
    rot_1 = rotmats_ctrl(:,:,idx_0+1);
    % interpolation parameter in [0,1]
    dt = (t_query - t_0) / (t_1 - t_0);
    % Linear interpolation, Lie group formulation
    axang_increm = rotm2axang(rot_0.'*rot_1);
    axang_increm(4) = axang_increm(4)*dt;
    rot_interp = rot_0 * axang2rotm( axang_increm );
end
