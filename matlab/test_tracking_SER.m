% DVS rotation tracking (localization) demo with an event camera, using an EKF
%
% Given DVS events and a panoramic image (i.e., map), recover the
% rotational motion of the camera that caused the events.
%
% Guillermo Gallego
% TU Berlin

clc; clear; close ALL


%% Dataset ________________________________________________________________

% % Synthetic dataset
% % Load DVS intrinsic calibration and undistortion map
% load DVS_synth_undistorted_pixels.mat
% sensor_height = 128; sensor_width = 128;
% folder = '../data/synth1';
% C_th = 0.45; % DVS Contrast threshold

% Input: given panoramic image (i.e., map)
load rec_image_Neumann.mat
map = rec_image;
fig_show_evol = figure();
imshow(map,[]); title('Map'); colorbar; hold on;

% Compute derivatives of the map (if not provided)
% [grad_map.x,grad_map.y] = imgradientxy(map); % Scaled derivative?
kernel_1 = 0.5*[-1 0 1]; % gradient kernel
grad_map.x = filter2(kernel_1 ,map); % vertical edges
grad_map.y = filter2(kernel_1',map); % horizontal edges


%% Main algorithmic parameters ____________________________________________
% Parameters of the EKF
sigma_covar = 1e-3;         % process noise covariance
scale_factor = 1e-3;        % for process noise covariance propagation
var_meas_noise = (0.17)^2;  % variance of measurement noise; units [C_th]^2
num_events_batch = 1000;    % for speed-up we process measurements in packets

%__________________________________________________________________________
% Output: array with the "trajectory" (history of rotations)
traj.time = 0; % tuples (t, rotation(t), covar(t))
traj.rotmat = zeros(9,1); % state mean (rotation matrix)
traj.covar = sigma_covar * reshape(eye(3),9,1);


%% Ground truth poses (for comparison)
Rot0 = rotmats_ctrl(:,:,1); % to center the map around the first pose

% Debugging: start somewhere in the middle, using ground truth
t0 = 0.001; % fails if t0=0. Need to DEBUG
idx_t = find(time_ctrl < t0, 1,'last');
t0 = time_ctrl(idx_t);
traj.time = time_ctrl(1:idx_t)';
traj.rotmat = zeros(9,idx_t);
for k = 1:idx_t
    Rot_mat = Rot0' * rotmats_ctrl(:,:,k);
    traj.rotmat(:,k) = Rot_mat(:);
    traj.covar(:,k) = sigma_covar * reshape(eye(3),9,1); % Debug: some random initialization
end


%% Tracking of camera orientation using an EKF
% Input: events, panoramic image (i.e., map), contrast threshold and camera calibration
% Output: camera orientation (at discrete times)

% profile('on','-detail','builtin','-timer','performance')
tic  % to measure execution time

num_events_display = 50000;
num_batches_display = floor(num_events_display / num_events_batch);

% For efficiency, a structure, called event map, contains for every sensor pixel
% the time of the last event and its rotation at that time.
%s.sae = -1e-6;    % not used for tracking
s.rotation = eye(3)*nan;
event_map = repmat(s, sensor_height, sensor_width);

first_plot = true; % for efficient plotting

iEv = 1; % event counter
iBatch = 1; % packet-of-events counter
iEv_last = 1;
while true
    
    if (iEv + num_events_batch > num_events)
        break; % There are no more events
    end
    
    % Get a batch of events
    events_batch = events(iEv + (0:num_events_batch-1),:);
    iEv = iEv + num_events_batch;
    iBatch = iBatch + 1;
    
    t_events_batch = events_batch(:,1);
    x_events_batch = events_batch(:,2);
    y_events_batch = events_batch(:,3);
    pol_events_batch = 2 * (events_batch(:,4) - 0.5);
    
    % Get index into 1-D array of pixel grid
    idx_to_mat = x_events_batch*sensor_height + y_events_batch + 1;

    % Assume all events in a batch share the same rotation (for speed-up)
    t_ev_mean = (t_events_batch(1) + t_events_batch(end)) * 0.5;
    
    if ( any(t_events_batch(:) < t0) )
        iEv_last = iEv;
        
        % Events in a batch are assigned the same rotation
        idx_0 = find( traj.time <= t_ev_mean, 1, 'last');
        Rot_prev0 = reshape(traj.rotmat(:,idx_0),3,3);
        if (idx_0 == numel(traj.time))
            Rot_prev = Rot_prev0;
        else
            % Linear interpolation of rotation
            idx_1 = idx_0 + 1;
            Rot_prev1 = reshape(traj.rotmat(:,idx_1),3,3);
            t_ctrl = traj.time([idx_0,idx_1]);
            Rot_prev = rotationAt(t_ctrl, cat(3, Rot_prev0, Rot_prev1), t_ev_mean);
        end
            
        % Debug: skip the first events but do some bookkeeping
        for ii=1:num_events_batch
            % Set rotations for each event
            event_map(idx_to_mat(ii)).rotation = Rot_prev;
        end
        
        continue;
    end
        
    % State and error covariance are updated using each event.
    % They are copied to the trajectory array every so many events or time 
    % (to avoid storing a rotation per event)
    if (iEv_last + num_events_batch == iEv)
        % Debug: provide values to test
        rotmat_cur = reshape(traj.rotmat(:,end),3,3);
        covar_cur = reshape(traj.covar(:,end),3,3);
        disp(['Last event before debugging: ' num2str(iEv)]);
    end
    
    % EKF
    % 1. Prediction / Propagation
    t_cur_state = t_ev_mean;
    t_last_update = traj.time(end);
    delta_t_state = t_cur_state - t_last_update; % time elapsed since last state update
    rotmat_pred = rotmat_cur; % update using motion model? Arbitrary motion
    % The covariance of the process noise should depend on the time elapsed
    % since the last measurement update: the longer, the larger the covariance.
    covar_process_noise = scale_factor * delta_t_state * eye(3);
    covar_pred = covar_cur + covar_process_noise;
    
    % 2. Correction / Update
    % Compute innovation and its covariance
    
    % Discard NaN values due to uninitialized rotation at previous time
    idx_notnan = idx_to_mat;
    pol_events_batch_notnan = pol_events_batch;
    mask_uninitialized = false(1,num_events_batch);
    for ii=1:num_events_batch
        if any(isnan(event_map(idx_to_mat(ii)).rotation(:)))
            mask_uninitialized(ii) = 1;
        end
    end
    num_uninitialized = sum(mask_uninitialized);
    if (num_uninitialized > 0)
        % Delete uninitialized events
        %disp(['deleting ' num2str(num_uninitialized) ' points'])
        idx_notnan(mask_uninitialized) = [];
        pol_events_batch_notnan(mask_uninitialized) = [];
    end
    
    % Compute innovation and Kalman gain
    % First, compute predicted contrast and analytical Jacobian
    f_contrast = @(x) compute_contrast_SO3(idx_notnan, event_map, x, map, ...
        undist_pix_calibrated, grad_map);
    [contrast_pred, df] = f_contrast(rotmat_pred);
    innov = C_th - pol_events_batch_notnan .* contrast_pred;
    G_jac = (pol_events_batch_notnan * [1 1 1]) .* df;
    
    % Discard NaN. Points that may fall out of the image 
    % Ideally one would pad the image before extrapolating pixel values
    % Just discard the event; there are many more. Padding is more expensive
    idx_nan = isnan(innov) | isnan(sum(G_jac,2));
    if any(idx_nan)
        innov(idx_nan) = [];
        G_jac(idx_nan,:) = [];
    end
    
    % Using the matrix inversion lemma we invert only a 3x3 matrix, 
    % and inverting the s*Id matrix is fast, too.
    num_innov = numel(innov);
    ivar_meas_noise = 1/var_meas_noise; % precision is inverse of variance
    S_inv_covar_innov = -(ivar_meas_noise*ivar_meas_noise) * ...
        G_jac * inv(inv(covar_pred) + ivar_meas_noise*(G_jac'*G_jac)) * G_jac';
    idx_diag = 1:(num_innov+1):(num_innov*num_innov);
    S_inv_covar_innov(idx_diag) = S_inv_covar_innov(idx_diag) + ivar_meas_noise;
    % Finally, compute Kalman gain
    Kalman_gain = covar_pred * G_jac' * S_inv_covar_innov;
    
    % Update rotation mean and covariance
    % See Barfoot's SER book, Section 8.2.4.
    rotmat_cur = expm_v_SO3(Kalman_gain * innov)*rotmat_pred;
    if any(isnan(rotmat_cur))
        break;
    end
    % Joseph's form covariance update requires computing S_covar_innovation;
    % the usual update is faster (and seems to preserve symmetry)
    covar_cur = covar_pred - Kalman_gain * G_jac * covar_pred;
    
    % Prepare for next events
    for ii = 1:num_events_batch
        % Update last rotation
        event_map(idx_to_mat(ii)).rotation = rotmat_cur;
    end
    
    % Store current rotation in state history array (trajectory)
    if (delta_t_state > 1e-3)
        disp(['Saving pose: Event # ' num2str(iEv)]);
        traj.time = [traj.time, t_cur_state];
        traj.rotmat = [traj.rotmat, rotmat_cur(:)];
        traj.covar = [traj.covar, covar_cur(:)];
    end
    
    if ( mod(iBatch, num_batches_display) == 0 )
        % Print current event number
        disp(['Update display: Event # ' num2str(iEv)]);
        
        %%% <!-- DEBUG. Visualization
        % Plot points projected on panoramic image using ground truth pose
        num_events_batch_plot = max([400,num_events_batch]);
        events_batch_plot = events(iEv + (0:num_events_batch_plot-1),:);
        t_events_batch_plot = events_batch_plot(:,1);
        x_events_batch_plot = events_batch_plot(:,2);
        y_events_batch_plot = events_batch_plot(:,3);
        t_ev_mean = (t_events_batch_plot(1) + t_events_batch_plot(end)) * 0.5;
        % Using ground truth rotation
        Rot = rotationAt(time_ctrl, rotmats_ctrl, t_ev_mean, f_r2a, f_a2r);
        % Get bearing vector of the event
        idx_to_mat_plot = x_events_batch_plot*sensor_height + y_events_batch_plot + 1;
        one_vec_plot = ones(num_events_batch_plot,1);
        bearing_vec = [undist_pix_calibrated(idx_to_mat_plot,:), one_vec_plot]'; % 3xN
        % Get map point corresponding to current event
        rotated_vec = Rot0' * Rot * bearing_vec;
        pm_gt = project_EquirectangularProjection(rotated_vec, pano_width, pano_height);

        % Using currently estimated rotation:
        Rot = rotmat_cur;
        % Get map point corresponding to current event
        rotated_vec = Rot * bearing_vec;
        pm = project_EquirectangularProjection(rotated_vec, pano_width, pano_height);

        figure(fig_show_evol);
        if first_plot 
            h_map_pts_gt = plot(pm_gt(1,:),pm_gt(2,:),'go');
            h_map_pts = plot(pm(1,:),pm(2,:),'r.');
            legend ('ground truth','estimated')
            first_plot = false;
        else
            set(h_map_pts_gt,'XData',pm_gt(1,:),'YData',pm_gt(2,:));
            set(h_map_pts,'XData',pm(1,:),'YData',pm(2,:));
        end
        drawnow
        %%% DEBUG -->
    end
    
end

toc  % measure execution time
% profile viewer


%% Visualize state and covariance

% Plot evolution of trace of the state covariance
figure, 
semilogy(traj.time, sqrt(sum(traj.covar([1,5,9],:))) * (180/pi) )
title('sqrt(Trace of the state covariance)')
xlabel('time'), ylabel('[deg]')
grid on

% Plot evolution of the state
figure, 
% First, get angles in the ball of radius pi before unwrapping
num = numel(traj.time);
rotvec_ballpi = zeros(3,num);
for k = 1:num
    axang = rotm2axang(reshape(traj.rotmat(:,k),3,3));
    rotvec_ballpi(:,k) = axang(1:3)' * axang(4);
end
% Actual plot
ax = gca;
ax.ColorOrderIndex = 1;
hold on, plot(traj.time, unwrap(rotvec_ballpi,[],2)'*(180/pi),'LineWidth',1)
title('Unwrapped angles vs time')
axis tight, grid on
xlabel('time'), ylabel('angle [deg]')
legend('\theta_1','\theta_2','\theta_3')

% Plot ground truth
gt_rotvec = zeros(3,num_poses_ctrl);
for k = 1:num_poses_ctrl
    ax = rotm2axang( Rot0' * rotmats_ctrl(:,:,k));
    gt_rotvec(:,k) = ax(1:3)' * ax(4);
end
ax = gca;
ax.ColorOrderIndex = 1;
plot(time_ctrl, unwrap(gt_rotvec,[],2)'*(180/pi),'--','LineWidth',2)

legend({'\theta_1','\theta_2','\theta_3',...
    'true \theta_1','true \theta_2','true \theta_3'},...
    'Location','southwest');