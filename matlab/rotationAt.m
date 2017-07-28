function rot_interp = rotationAt(t_ctrl, rotmats_ctrl, t_query, f_r2a, f_a2r)
% rotationAt Function for linear interpolation of rotations
% Interpolate the orientation (rotation) at a given time.
%
% Input:
% -t_ctrl: timestamps of discrete set of orientations ("control poses")
% -rotmats_ctrl: discrete set of rotation matrices
% -t_query: time of the requested rotation matrix
% -f_r2a: (Optional) handle to function that converts rotation to axis-angle
% -f_a2r: (Optional) handle to function that converts axis-angle to rotation
%
% Output:
% -rot_interp: interpolated rotation matrix
%
% I/O in the style of interp1

% Check function handles
if (nargin < 5) || ~isa(f_r2a,'function_handle') || ~isa(f_a2r,'function_handle')
    
    % If no conversion functions are provided, specify them
    
    if license('test', 'Robotics_System_Toolbox')
        % Requires Robotic System toolbox
        f_r2a = @rotm2axang;
        f_a2r = @axang2rotm;
        
    else
        % Use VR toolbox if Robotics System is not available
        toolboxes_info = ver;
        if any(strcmp('Simulink 3D Animation', {toolboxes_info.Name}))
            % Requires VR toolbox
            f_r2a = @vrrotmat2vec;
            f_a2r = @vrrotvec2mat;
            
        else
            % Use custom conversion functions, which should be in the path
            f_r2a = @R2AA;
            f_a2r = @AA2R;
        end
    end
end


% Rotation interpolation
if (t_query < t_ctrl(1) || t_query > t_ctrl(end) )
    rot_interp = nan * rotmats_ctrl(1,:);
    
elseif (t_query == t_ctrl(end))
    rot_interp = rotmats_ctrl(end,:);
    
else
    idx_0 = find(t_query >= t_ctrl,1,'last');
    % Two rotations and their times
    t_0 = t_ctrl(idx_0);
    t_1 = t_ctrl(idx_0+1);
    rot_0 = rotmats_ctrl(:,:,idx_0);
    rot_1 = rotmats_ctrl(:,:,idx_0+1);
    % Interpolation parameter in [0,1]
    dt = (t_query - t_0) / (t_1 - t_0);
    % Linear interpolation, Lie group formulation
    axang_increm = f_r2a( rot_0.'*rot_1 );
    axang_increm(4) = axang_increm(4)*dt;
    rot_interp = rot_0 * f_a2r( axang_increm );
end