% DVS image reconstruction demo
%
% Given DVS events and camera poses (rotational motion), reconstruct the
% gradient map that caused the events.
%
% Guillermo Gallego
% Robotics and Perception Group
% University of Zurich

clc; clear; close ALL


%% Check for toolboxes and integration methods
% Select the integration method. Options are:
%   'poisson_dirichlet'    % Requires the PDE toolbox
%   'poisson_neumann'      % Requires the Image Processing toolbox
%   'frankotchellappa'     % Does not require a toolbox
integration_method = 'frankotchellappa';

% Check which toolboxes are available
% NOTE:  Returning '1' means that a license for the toolbox is present, but
%        it does guarantee that a license server will allow the license to
%        be checked out
use_PDE_toolbox = license('test', 'PDE_toolbox');
use_image_processing_toolbox = license('test', 'Image_Toolbox');
use_robotics_system_toolbox = license('test', 'Robotics_System_Toolbox');
% Use VR toolbox if Robotics System is not available
toolboxes_info = ver;
use_VR_toolbox = any(strcmp('Simulink 3D Animation', {toolboxes_info.Name}));

addpath('integration_methods') % Path to different integration methods

% Provide function handles to convert from rotation matrix to axis-angle
% and vice-versa
if use_robotics_system_toolbox
    % Use Robotic System toolbox
    f_r2a = @rotm2axang;
    f_a2r = @axang2rotm;

elseif use_VR_toolbox
    % Use VR toolbox
    f_r2a = @vrrotmat2vec;
    f_a2r = @vrrotvec2mat;
    addpath('coordinate_transforms') % Need q2R function
    
else
    % Use custom functions
    addpath('coordinate_transforms')
    f_r2a = @R2AA;
    f_a2r = @AA2R;
end

% Check that the integration method is compatible with the available
% toolboxes and give some meaningful feedback if it isn't
if strcmp(integration_method, 'poisson_dirichlet')
    % Requires the Partial Differential Equation toolbox
    if ~use_PDE_toolbox
        error(['Cannot use option "PDE_poicalc" without the PDE toolbox. ' ...
            'Please select another integration method']);
    end
elseif strcmp(integration_method, 'poisson_neumann')
    % Requires the Image Processing toolbox
    if ~use_image_processing_toolbox
        error(['Cannot use option "poisson_neumann" without the image processing toolbox. ' ...
            'Please select another integration method']);
    end
elseif strcmp(integration_method, 'frankotchellappa')
    % Does not require a toolbox.
    % Do nothing (this method does not require a toolbox)
else
    error([integration_method ' is not a valid integration method. ' ...
        'Valid options are: \n\t poisson_dirichlet  \n\t poisson_neumann \n\t frankotchellappa']);
end


%% Dataset

% Synthetic dataset
% Load DVS intrinsic calibration and undistortion map
load DVS_synth_undistorted_pixels.mat
sensor_height = 128; sensor_width = 128;
folder = '../data/synth1';
C_th = 0.45; % DVS Contrast threshold


%% Load events
disp('Loading events');
filename_events = fullfile(folder,'events.txt');
events = load(filename_events);
num_events = size(events,1);
disp(['Number of events in file: ' num2str(num_events)]);

% Remove time of offset
first_event_sec = events(1,1);
first_event_nsec = events(1,2);
events_time = (events(:,1)-first_event_sec) + 1e-9 * (events(:,2) - first_event_nsec);
events = [events_time, events(:,3:5)]; % t x y pol


%% Load camera poses
disp('Loading camera orientations');
filename_poses = fullfile(folder,'poses.txt');
gt = load(filename_poses);
time_ctrl = (gt(:,1) - first_event_sec) + 1e-9*(gt(:,2) - first_event_nsec);
poses_ctrl = gt(:,3:end);

num_poses_ctrl = size(poses_ctrl,1);
quats = poses_ctrl(:,4:7); % keep only the rotational part (quaterions): qx qy qz qw
quats = quats(:,[4,1:3]); % qw,qx,qy,qz same as Matlab's order

% Convert quaternions to rotation matrices
rotmats_ctrl = zeros(3,3,num_poses_ctrl);
for k = 1:num_poses_ctrl
    if use_robotics_system_toolbox
        % Requires Robotics System Toolbox
        rotmats_ctrl(:,:,k) = quat2rotm(quats(k,:));
    else
        % Removes the dependency on the Robotics System Toolbox
        rotmats_ctrl(:,:,k) = q2R(quats(k,:));
    end
end




%% Image reconstruction using pixel-wise EKF
% Input: events, contrast threshold and camera orientation (discrete poses + interpolation)
% Output: reconstructed intensity image

% profile('on','-detail','builtin','-timer','performance')
tic  % to measure execution time

num_events_batch = 300;
num_events_display = 100000;
num_batches_display = floor(num_events_display / num_events_batch);

% Extended Kalman Filter (EKF) parameters
% CHOOSE one measurement function: "contrast" or "event rate"
measurement_criterion = 'contrast';
% measurement_criterion = 'event rate';

if ( strcmpi(measurement_criterion,'contrast') )
    var_R = (0.17).^2; % Measurement variance using constrast criterion (Gallego et al. arXiv 2015)
else
    var_R = 1e4; % Measurement variance using event rate criterion (Kim et al. BMVC 2014)
end
grad_initial_variance = 10;


% Variables related to the reconstructed image mosaic (panorama)
pano_height = 1024;
pano_width = 2 * pano_height;
% Gradient map
grad_map.x = zeros(pano_height, pano_width);
grad_map.y = zeros(pano_height, pano_width);
% Covariance matrix of each gradient pixel
grad_map_covar.xx =  ones(pano_height, pano_width) * grad_initial_variance;
grad_map_covar.xy = zeros(pano_height, pano_width);
grad_map_covar.yx = zeros(pano_height, pano_width);
grad_map_covar.yy =  ones(pano_height, pano_width) * grad_initial_variance;

% For efficiency, a structure, called event map, contains for every pixel
% the time of the last event and its rotation at that time.
s.sae = -1e-6;
s.rotation = eye(3)*nan;
event_map = repmat(s, sensor_height, sensor_width);

Rot0 = rotmats_ctrl(:,:,1); % to center the map around the first pose
one_vec = ones(num_events_batch,1);
Id = eye(2);

disp('Processing events');

close ALL
fig_show_evol = figure('Color','w', 'units','normalized','outerposition',[0 0 1 1]);
first_plot = true; % for efficient plotting

iEv = 1; % event counter
iBatch = 1; % packet-of-events counter
while true
    
    if (iEv + num_events_batch > num_events)
        break; % There are no more events
    end
    
    % Get batch of events
    events_batch = events(iEv + (0:num_events_batch-1),:);
    iEv = iEv + num_events_batch;
    
    t_events_batch = events_batch(:,1);
    x_events_batch = events_batch(:,2);
    y_events_batch = events_batch(:,3);
    pol_events_batch = 2 * (events_batch(:,4) - 0.5);
    
    
    %----------------------------------------------------------------------
    % Get the two map points corresponding to each event,
    % and update the event map (time and rotation of last event)
    
    % Get time of previous event at same DVS pixel
    idx_to_mat = x_events_batch*sensor_height + y_events_batch + 1;
    t_prev_batch = [event_map(idx_to_mat).sae].';
    
    % Get (interpolated) rotation of current event
    t_ev_mean = (t_events_batch(1) + t_events_batch(end)) * 0.5;
    if (t_ev_mean > time_ctrl(end))
        break; % event later than last known pose
    end
    Rot = rotationAt(time_ctrl, rotmats_ctrl, t_ev_mean, f_r2a, f_a2r);
    
    % Get bearing vector of the event
    bearing_vec = [undist_pix_calibrated(idx_to_mat,:), one_vec].'; % 3xN
    
    % Get map point corresponding to current event
    rotated_vec = Rot0.' * Rot * bearing_vec;
    pm = project_EquirectangularProjection(rotated_vec, pano_width, pano_height);
    
    % Get map point corresponding to previous event at same pixel
    rotated_vec_prev = zeros(size(rotated_vec));
    for ii = 1:num_events_batch
        Rot_prev = event_map(idx_to_mat(ii)).rotation;
        rotated_vec_prev(:,ii) = Rot0.' * Rot_prev * bearing_vec(:,ii);
        % Update last rotation and time of event (SAE)
        event_map(idx_to_mat(ii)).sae = t_events_batch(ii);
        event_map(idx_to_mat(ii)).rotation = Rot;
    end
    pm_prev = project_EquirectangularProjection(rotated_vec_prev, pano_width, pano_height);
    
    if (t_prev_batch(end) < 0) || (t_prev_batch(end) < time_ctrl(1))
        continue; % initialization phase. Fill in event_map
    end
    
    % Discard nan values
    mask_uninitialized = isnan(pm_prev(1,:)) | isnan(pm_prev(2,:));
    num_uninitialized = sum(mask_uninitialized);
    if (num_uninitialized > 0)
        % Delete uninitialized events
        % disp(['deleting ' num2str(num_uninitialized) ' points'])
        t_events_batch(mask_uninitialized) = [];
        t_prev_batch(mask_uninitialized) = [];
        pol_events_batch(mask_uninitialized) = [];
        pm(:,mask_uninitialized) = [];
        pm_prev(:,mask_uninitialized) = [];
    end
    
    % Get time since previous event at same pixel
    tc = t_events_batch - t_prev_batch;
    event_rate = 1 ./ (tc + 1e-12); % measurement or observation (z)
    
    % Get velocity vector
    vel = (pm - pm_prev) .* ([1;1]*event_rate.');
    
    
    %----------------------------------------------------------------------
    % Extended Kalman Filter (EKF) for the intensity gradient map.
    % Get gradient and covariance at current map points pm
    ir = floor(pm(2,:)); % row is y coordinate
    ic = floor(pm(1,:)); % col is x coordinate
    idx_map = (ic*pano_height + ir + 1).';
    gm = [grad_map.x(idx_map), grad_map.y(idx_map)];
    Pg = [grad_map_covar.xx(idx_map), grad_map_covar.xy(idx_map), ...
        grad_map_covar.yx(idx_map), grad_map_covar.yy(idx_map)];
    
    % EKF update
    if ( strcmpi(measurement_criterion,'contrast') )
        % Use contrast as measurement function
        dhdg = vel.' .* ((tc .* pol_events_batch)*[1,1]); % deriv. of measurement function
        nu_innovation = C_th - sum(dhdg .* gm,2);
    else
        % Use the event rate as measurement function
        dhdg = vel.' ./ ((C_th * pol_events_batch)*[1,1]); % deriv. of measurement function
        nu_innovation = event_rate - sum(dhdg .* gm,2);
    end
    Pg_dhdg = [Pg(:,1).*dhdg(:,1) + Pg(:,2).*dhdg(:,2), ...
        Pg(:,3).*dhdg(:,1) + Pg(:,4).*dhdg(:,2)];
    S_covar_innovation = (dhdg(:,1).*Pg_dhdg(:,1) + dhdg(:,2).*Pg_dhdg(:,2)) + var_R;
    Kalman_gain = Pg_dhdg ./ (S_covar_innovation * [1,1]);
    % Update gradient and covariance
    gm = gm + Kalman_gain .* (nu_innovation * [1,1]);
    Pg = Pg - [Pg_dhdg(:,1).*Kalman_gain(:,1), Pg_dhdg(:,1).*Kalman_gain(:,2), ...
        Pg_dhdg(:,2).*Kalman_gain(:,1), Pg_dhdg(:,2).*Kalman_gain(:,2)];
    
    % Store updated values
    grad_map.x(idx_map) = gm(:,1);
    grad_map.y(idx_map) = gm(:,2);
    grad_map_covar.xx(idx_map) = Pg(:,1);
    grad_map_covar.xy(idx_map) = Pg(:,2);
    grad_map_covar.yx(idx_map) = Pg(:,3);
    grad_map_covar.yy(idx_map) = Pg(:,4);
    
    
    %----------------------------------------------------------------------
    % Visualization every so many batches of events
    iBatch = iBatch + 1;
    if ( mod(iBatch, num_batches_display) == 0 )
        disp(['Update display: Event # ' num2str(iEv)]); % display current event number
        figure(fig_show_evol);
        
        % Compute trace image (estimated uncertainty)
        trace_map = grad_map_covar.xx + grad_map_covar.yy;
        
        % Compute reconstructed image (Poisson solver)
        % high resolution and high dynamic range
        % (MATLAB code uses Dirichlet b.c.)
        grad_map_clip.x = grad_map.x;
        grad_map_clip.y = grad_map.y;
        mask = trace_map > 0.05; % reconstruct only gradients with small covariance
        grad_map_clip.x(mask) = 0;
        grad_map_clip.y(mask) = 0;
        
        % The lines below are different integration methods
        if strcmp(integration_method, 'poisson_dirichlet')
            % Method 1: Dirichlet boundary conditions
            % Requires the Partial Differential Equation toolbox
            div = divergence(grad_map_clip.x,grad_map_clip.y);
            rec_image = reshape( poicalc(-div(:),1,1,pano_height,pano_width), pano_height,pano_width);
            
        elseif strcmp(integration_method, 'poisson_neumann')
            % Method 2: Neumann boundary conditions
            % Requires the Image Processing toolbox
            % The code is from http://www.cs.cmu.edu/~ILIM/projects/IM/aagrawal/software.html
            rec_image = poisson_solver_function_neumann(grad_map_clip.x, grad_map_clip.y);
            
        elseif strcmp(integration_method, 'frankotchellappa')
            % Method 3: Neumann boundary conditions
            % Does not require a toolbox
            % The code is from http://www.cs.cmu.edu/~ILIM/projects/IM/aagrawal/software.html
            rec_image = frankotchellappa(grad_map_clip.x, grad_map_clip.y);
            rec_image = rec_image - mean(rec_image(:));
        end
        
        if first_plot
            subplot(2,2,1),
            % Plot points on panoramic image, colored according to polarity
            idx_pos = pol_events_batch > 0;
            idx_neg = pol_events_batch < 0;
            h_map_pts = plot(pm(1,idx_pos),pm(2,idx_pos),'.b', ...
                pm(1,idx_neg),pm(2,idx_neg),'.r');
            grid on, axis ij equal
            axis([0 pano_width 0 pano_height])
            title('Map points from events')
            set(gca,'FontSize',20);
            
            % Display reconstructed image
            subplot(2,2,2), h_img = imshow(rec_image / max(abs(rec_image(:))),[-1,1]);
            title('Reconstructed image')
            set(gca,'FontSize',20);
            
            % Display one of the gradient images
            subplot(2,2,3), h_gx = imshow(grad_map.x / std(grad_map.x(:)),[-5,5]);
            title('Gradient in X direction')
            set(gca,'FontSize',20);
            
            % Display trace of the error covariance
            subplot(2,2,4), h_trace = imshow(trace_map / max(abs(trace_map(:))),[0 1]);
            title('Trace of covariance')
            set(gca,'FontSize',20);
            
            first_plot = false;
        else
            subplot(2,2,1),
            idx_pos = pol_events_batch > 0;
            set(h_map_pts(1),'XData',pm(1,idx_pos),'YData',pm(2,idx_pos));
            idx_neg = pol_events_batch < 0;
            set(h_map_pts(2),'XData',pm(1,idx_neg),'YData',pm(2,idx_neg));
            
            subplot(2,2,2), set(h_img,'CData',rec_image / max(abs(rec_image(:))));
            subplot(2,2,3), set(h_gx,'CData',grad_map.x / std(abs(grad_map.x(:))));
            subplot(2,2,4), set(h_trace,'CData',trace_map / max(abs(trace_map(:))));
        end
        drawnow
    end
    
end

toc  % measure execution time
% profile viewer


%% Display in separate figure
rec_image_normalized = rec_image / max(abs(rec_image(:)));
figure('Name','Reconstructed image (log scale)');
imshow(rec_image_normalized,[])
filename_out = fullfile(folder,'reconstructed_image_log_scale.jpg');
imwrite(rec_image_normalized+0.5,filename_out,'Quality',90);

rec_image_exp = exp(0.001 + rec_image);
figure('Name','Reconstructed image (linear scale)');
imshow(rec_image_exp,[0,5])

trace_map_normalized = trace_map / max(trace_map(:));
fig_trace = figure('Name','Trace of covariance');
imshow(trace_map,[]), colorbar
filename_out = fullfile(folder,'covariance_trace.jpg');
imwrite(trace_map_normalized,filename_out,'Quality',90);


%% Save gradient images
num_sigmas = 5; % clipping

% x-component of the gradient
gx_normalized = grad_map.x / std(abs(grad_map.x(:)));
gx_normalized(gx_normalized > num_sigmas) = num_sigmas;
gx_normalized(gx_normalized < -num_sigmas) = -num_sigmas;
gx_normalized = (gx_normalized + num_sigmas)/(2*num_sigmas);
filename_out = fullfile(folder,'mosaicing_gx.jpg');
imwrite(gx_normalized,filename_out,'Quality',90);
% y-component of the gradient
gy_normalized = grad_map.y / std(abs(grad_map.y(:)));
gy_normalized(gy_normalized > num_sigmas) = num_sigmas;
gy_normalized(gy_normalized < -num_sigmas) = -num_sigmas;
gy_normalized = (gy_normalized + num_sigmas)/(2*num_sigmas);
filename_out = fullfile(folder,'mosaicing_gy.jpg');
imwrite(gy_normalized,filename_out,'Quality',90);


%% Colored gradient using HSV according to magnitude and direction
% Combine gradient magnitude and direction into a single image.
% Normalize for color coding

if use_image_processing_toolbox
    % Requires the Image Processing Toolbox
    [g_grad,g_ang] = imgradient(grad_map.x,grad_map.y);
else
    % Same as above, but remove the need for the Image Processing Toolbox
    g_ang = -atan2d(grad_map.y,grad_map.x);
    g_grad = sqrt(grad_map.x.^2+grad_map.y.^2);
end

g_grad_unit = g_grad/1;
g_grad_unit(g_grad_unit > 1.0) = 1.0;
g_ang_unit = g_ang/360 + 0.5;
hsv_img = cat(3, g_ang_unit, g_grad_unit, ones(size(g_ang)));
rgb_edges = hsv2rgb(hsv_img);
figure('Name','Gradient (magnitude and direction)')
imshow(rgb_edges,[])
filename_out = fullfile(folder,'mosaicing_grad_map_hsv.png');
imwrite(rgb_edges,filename_out);


%% Pseudo-colored trace ("confidence map")
figure('Color','w','Name','Trace of covariance (in log10 scale)');
cmap = flipud(colormap('jet'));
imshow(log10(trace_map + 1e-3),[], 'Colormap',cmap); colorbar
saveas(gcf, fullfile(folder,'covariance_trace_colored_log_scale'), 'jpg');
