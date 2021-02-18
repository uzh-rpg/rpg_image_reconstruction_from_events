% Rotation tracking with an event camera, using an EKF approach

clear; close ALL; clc

% Prototyping numerical derivatives
addpath('/home/ggb/pCloudDrive/MATLAB/MVGToolbox/General/Derivatives');


% % Synthetic dataset
% % Load DVS intrinsic calibration and undistortion map
% load DVS_synth_undistorted_pixels.mat
% sensor_height = 128; sensor_width = 128;
% folder = '../data/synth1';
% C_th = 0.45; % DVS Contrast threshold

% Given: map
load rec_image_Neumann.mat
map = rec_image;
fig_show_evol = figure();
imshow(map,[]); title('Map'); colorbar; hold on;

% Compute derivatives of the map (if not provided)
% [Gx,Gy] = imgradientxy(map); % Scaled derivative?
% grad_map.x = Gx;
% grad_map.y = Gy;
kernel_1 = 0.5*[-1 0 1]; % gradient kernel
grad_map.x = filter2(kernel_1  ,map); % vertical edges
grad_map.y = filter2(kernel_1.',map); % horizontal edges

% Array with the "trajectory" (history of rotations)
traj.time = 0; % tuples (t, rotation(t), covar(t))
traj.rotmat = zeros(9,1); % state mean (rotation matrix)
sigma_covar = 1e-3;
traj.covar = sigma_covar * reshape(eye(3),9,1);
var_meas_noise = (0.17)^2; % variance of measurement noise
scale_factor = 1e-3; % scale factor for process noise covariance propagation


num_events_batch = 1000;
num_events_display = 50000;
num_batches_display = floor(num_events_display / num_events_batch);


%%
Rot0 = rotmats_ctrl(:,:,1); % to center the map around the first pose

s.sae = -1e-6;
s.rotation = eye(3)*nan;
event_map = repmat(s, sensor_height, sensor_width);

% Debugging: start somewhere in the middle, using ground truth
t0 = 0.001; % fails of t0=0. Need to DEBUG
idx_t = find(time_ctrl < t0, 1,'last');
t0 = time_ctrl(idx_t);
traj.time = time_ctrl(1:idx_t)';
traj.rotmat = zeros(9,idx_t);
for k = 1:idx_t
    Rot_mat = Rot0' * rotmats_ctrl(:,:,k);
    traj.rotmat(:,k) = Rot_mat(:);
    traj.covar(:,k) = sigma_covar * reshape(eye(3),9,1); % Debug: some random initialization
end


%%

% profile('on','-detail','builtin','-timer','performance')
tic  % to measure execution time

first_plot = true; % for efficient plotting

iEv = 1; % event counter
iBatch = 1; % packet-of-events counter
iEv_last = 1;
while true
%while iEv < 500000 
    
    %if (iEv + num_events_batch > num_events/6)
    if (iEv + num_events_batch > num_events)
        break; % There are no more events
    end
    
    % Get batch of events
    events_batch = events(iEv + (0:num_events_batch-1),:);
    iEv = iEv + num_events_batch;
    iBatch = iBatch + 1;
    
    t_events_batch = events_batch(:,1);
    x_events_batch = events_batch(:,2);
    y_events_batch = events_batch(:,3);
    pol_events_batch = 2 * (events_batch(:,4) - 0.5);
    
    % Get index into 1-D array of pixel grid
    idx_to_mat = x_events_batch*sensor_height + y_events_batch + 1;

    % Some assigned time to all events in the batch
    t_ev_mean = (t_events_batch(1) + t_events_batch(end)) * 0.5;
    
    if ( any(t_events_batch(:) < t0) )
        iEv_last = iEv;
        
        % Set rotations for each event. Simpler: all events are
        % assigned the same rotation
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
            % Update the sae
            event_map(idx_to_mat(ii)).sae = t_events_batch(ii);
            
            % Set rotations for each event
            event_map(idx_to_mat(ii)).rotation = Rot_prev;
        end
        
        continue;
    end
        
%     iEv
%     %%% <!-- DEBUG
%     figure(100);
%     keyboard
%     imshow(map,[]);
%     num_events_batch_plot = 600;
%     events_batch_plot = events(iEv + (0:num_events_batch_plot-1),:);
%     t_events_batch_plot = events_batch_plot(:,1);
%     x_events_batch_plot = events_batch_plot(:,2);
%     y_events_batch_plot = events_batch_plot(:,3);
%     t_ev_mean = (t_events_batch_plot(1) + t_events_batch_plot(end)) * 0.5;
%     % Using ground truth rotation
%     Rot = rotationAt(time_ctrl, rotmats_ctrl, t_ev_mean, f_r2a, f_a2r);
%     % Get bearing vector of the event
%     idx_to_mat_plot = x_events_batch_plot*sensor_height + y_events_batch_plot + 1;
%     one_vec_plot = ones(num_events_batch_plot,1);
%     bearing_vec = [undist_pix_calibrated(idx_to_mat_plot,:), one_vec_plot]'; % 3xN
%     % Get map point corresponding to current event
%     rotated_vec = Rot0' * Rot * bearing_vec;
%     pm = project_EquirectangularProjection(rotated_vec, pano_width, pano_height);
%     hold on, plot(pm(1,:),pm(2,:),'g.')
%     daspect([1 1 1]), axis tight;
%     
%     % Using currently estimated rotation:
%     Rot = expm( Cross2Matrix(rotvec_cur) );
%     % Get map point corresponding to current event
%     rotated_vec = Rot0' * Rot * bearing_vec;
%     pm = project_EquirectangularProjection(rotated_vec, pano_width, pano_height);
%     hold on, plot(pm(1,:),pm(2,:),'r.')
%     daspect([1 1 1]), axis tight;
%     %%% DEBUG -->

    
    % State and error covariance are updated on every event but are copied
    % to the trajectory array only every so many events or time (to avoid
    % storing a rotation per event)
    if (iEv_last + num_events_batch == iEv)
        % Debug: provide values to test
        rotmat_cur = reshape(traj.rotmat(:,end),3,3);
        covar_cur = reshape(traj.covar(:,end),3,3);
        disp('Last event before debugging:'); iEv
    end
    
    % EKF
    % 1. Prediction / Propagation
    t_cur_state = t_ev_mean;
    t_last_update = traj.time(end);
    delta_t_state = t_cur_state - t_last_update; % time elapsed since last measurement update
    rotmat_pred = rotmat_cur; % update with motion model? not needed
    % disp(['satisfy rot mat constraint?: ', num2str(norm(rotmat_pred'*rotmat_pred-eye(3),'fro'))]);
    % The covariance of the process noise should depend on the time elapsed
    % since the last measurement update: the longer, the larger the covariance
    % should be.
    covar_process_noise = scale_factor * delta_t_state * eye(3);
    covar_pred = covar_cur + covar_process_noise;
    
    % 2. Correction / Update
    % Compute innovation and its covariance:
    %t_prev_batch = [event_map(idx_to_mat).sae]';

%     %%%<!-- DEBUG
%     % 1. Get bearing vector of the event pixel
%     idx_to_mat_deb = x_events_batch * sensor_height + y_events_batch + 1;
%     one_vec =  ones(1,size(event,1));
%     bearing_vec = [undist_pix_calibrated(idx_to_mat_deb,:), one_vec]'; % 3xN
%     % 2a. Get current rotation and point on the map
%     Rot = expm( Cross2Matrix(rotvec_pred) );
%     rotated_vec = Rot * bearing_vec;
%     p_cur  = project_EquirectangularProjection(rotated_vec, pano_width, pano_height);
%     figure(100); hold on, plot(p_cur(1,:),p_cur(2,:),'c*')
%     % 2b. Get previous rotation and point on the map
%     % get previous rotation
%     idx_0 = find( traj.time <= t_prev_event, 1, 'last');
%     Rot_prev0 = expm( Cross2Matrix(traj.rotvec(:,idx_0)) );
%     if (idx_0 == numel(traj.time))
%         Rot_prev = Rot_prev0;
%     else
%         % Linear interpolation of rotation
%         idx_1 = idx_0 + 1;
%         Rot_prev1 = expm( Cross2Matrix(traj.rotvec(:,idx_1)) );
%         t_ctrl = traj.time([idx_0,idx_1]);
%         rotmats_ctrl = cat(3, Rot_prev0, Rot_prev1);
%         Rot_prev = rotationAt(t_ctrl, rotmats_ctrl, t_prev_event);
%     end
%     % get corresponding point on the map
%     rotated_vec_prev = Rot_prev * bearing_vec;
%     p_prev = project_EquirectangularProjection(rotated_vec_prev, pano_width, pano_height);
%     figure(100); hold on, plot(p_prev(1,:),p_prev(2,:),'co')
%     M = interp2(1:pano_width, 1:pano_height, map, [p_cur(1); p_prev(1)], [p_cur(2); p_prev(2)]);
%     contrast = M(1) - M(2);
%     %%% -->
    
%     f_contrast = @(x) compute_contrast(event, t_prev_batch, x, traj, map, ...
%         sensor_height, undist_pix_calibrated);

    % Discard nan values due to uninitialized rotation at previous time
    idx_to_mat_notnan = idx_to_mat;
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
        idx_to_mat_notnan(mask_uninitialized) = [];
        pol_events_batch_notnan(mask_uninitialized) = [];
    end
    % Compute innovation and Kalman gain
    
%     % Numerical derivative
%     f_contrast_num = @(x) compute_contrast_SO3(idx_to_mat_notnan, event_map, x, map, undist_pix_calibrated);
%     contrast_pred = f_contrast_num(rotmat_pred);
%     Jacobians = fdjac_SO3(rotmat_pred, f_contrast_num);
%     contrast_pred_num = contrast_pred;
%     Jacobians_num = Jacobians;
     
    % With analytical derivative (contrat wrt state)
    f_contrast = @(x) compute_contrast_SO3(idx_to_mat_notnan, event_map, x, map, undist_pix_calibrated, grad_map);
    [contrast_pred, Jacobians] = f_contrast(rotmat_pred);
    
%     norm(contrast_pred_num - contrast_pred) / norm(contrast_pred_num)
%     norm(Jacobians_num - Jacobians,'fro') / norm(Jacobians_num,'fro')
    
    innov = ( C_th - pol_events_batch_notnan .* contrast_pred );
    grad = (pol_events_batch_notnan * [1 1 1]) .* Jacobians;
    
    
    % Discard nan. Points that may fall out of the image 
    % Ideally one would pad the image before extrapolating pixel values
    % Just discard the event; there are many more. Padding is more expensive
    idx_nan = isnan(innov) | isnan(sum(grad,2));
    if any(idx_nan)
        innov(idx_nan) = [];
        grad(idx_nan,:) = [];
    end
    dqdstate = grad;
    %S_covar_innovation = var_meas_noise * eye(numel(innov)) + dqdstate * covar_pred * dqdstate';
    %S_covar_innovation = var_meas_noise * speye(numel(innov)) + (dqdstate * covar_pred) * dqdstate';
    
%     % Faster sum?
%     S_covar_innovation = dqdstate * covar_pred * dqdstate';
    num_innov = numel(innov);
%     S_covar_innovation(1:(num_innov+1):end) = S_covar_innovation(1:(num_innov+1):end) + var_meas_noise;

    %Kalman_gain = covar_pred * dqdstate' / S_covar_innovation;
    
    % Using the matrix inversion lemma we invert only a 3x3 matrix, 
    % and inverting the s*Id matrix is fast, too.
    precision_meas_noise = 1/var_meas_noise; % precision is inverse of variance
    
%     A_inv = precision_meas_noise * speye(numel(innov));
%     inv_S_covar_innovation = A_inv - (precision_meas_noise*precision_meas_noise) * ...
%         dqdstate * inv(inv(covar_pred) + precision_meas_noise*(dqdstate'*dqdstate)) * dqdstate';

    % faster
    inv_S_covar_innovation = - (precision_meas_noise*precision_meas_noise) * ...
        dqdstate * inv(inv(covar_pred) + precision_meas_noise*(dqdstate'*dqdstate)) * dqdstate';
    inv_S_covar_innovation(1:(num_innov+1):end) = inv_S_covar_innovation(1:(num_innov+1):end) + precision_meas_noise;
    
    Kalman_gain = covar_pred * dqdstate' * inv_S_covar_innovation;
    
    % Update rotation vector and covariance
    increment_state = Kalman_gain * innov;
    rotmat_cur = expm_v_SO3(increment_state)*rotmat_pred;
    if any(isnan(rotmat_cur))
        break;
    end
    
    % Joseph form requires to compute S_covar_innovation
    %covar_cur = covar_pred - Kalman_gain * S_covar_innovation * Kalman_gain';
    
    covar_cur = covar_pred - Kalman_gain * dqdstate * covar_pred;
    
    % Update last rotation and time of event (SAE)
    Rot = rotmat_cur;
    for ii = 1:num_events_batch
        event_map(idx_to_mat(ii)).sae = t_events_batch(ii);
        event_map(idx_to_mat(ii)).rotation = Rot;
    end
    
    % Store current state (rotation vector) in state history array (trajectory)
    if (delta_t_state > 1e-3)
        disp(['Storing control pose: Event # ' num2str(iEv)]);
        traj.time = [traj.time, t_cur_state];
        traj.rotmat = [traj.rotmat, rotmat_cur(:)];
        traj.covar = [traj.covar, covar_cur(:)];
    end
    
    if ( mod(iBatch, num_batches_display) == 0 )
        disp(['Update display: Event # ' num2str(iEv)]); % display current event number
        

        %%% <!-- DEBUG
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

legend({'\theta_1','\theta_2','\theta_3','true \theta_1','true \theta_2','true \theta_3'},'Location','southwest')
