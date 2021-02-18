function [contrast, Jac] = compute_contrast(idx_to_mat, event_map, rotvec_pred, map, ...
    undist_pix_calibrated, grad_map)
% compute_contrast Measurement function for EKF correction step
%
% Compute the contrast of a set of events using current rotation and past
% rotations. Contrast is computed by difference of map intensities.
%
% Input:
% -idx_to_mat(N,1): event locations using indices into 1-D array of undistorted 2D coordinates
% -event_map: DVS pixel map of times and rotations of previous events
% -rotvec_pred(3,1): predicted state (rotation, before measurement update)
% -map: mosaic image, containing brightness information
% -undist_pix_calibrated(Npix,2): look-up-table of undistorted 2D coordinates
% -grad_map: derivatives of the map, to compute analytical derivatives
%
% Output:
% -contrast(N,1): predicted contrast of current event(s) according to two
%            rotations and the map
% -Jac(N,3): each row is the gradient of the contrast with respect to the
%            current state (rotation vector)

% Minimal check
assert(numel(rotvec_pred)==3);

num_events_batch = numel(idx_to_mat);
[pano_height,pano_width] = size(map);

% 1. Get bearing vector of the event pixel
one_vec = ones(num_events_batch,1);
bearing_vec = [undist_pix_calibrated(idx_to_mat,:), one_vec]'; % 3xN

% 2a. Get map point corresponding to current event
Rot = expm_v_SO3(rotvec_pred);
rotated_bvec = Rot * bearing_vec;
if (nargout > 1)
    [pm,dpm_d3D] = project_EquirectangularProjection(rotated_bvec, pano_width, pano_height);
else
    pm = project_EquirectangularProjection(rotated_bvec, pano_width, pano_height);
end

% 2b. Get map point corresponding to previous event at same pixel
% Rotated vectors: manual matrix-vector multiplication is faster than alternatives
rbv = [event_map(idx_to_mat).rotation] .* ...
    (ones(3,1)*reshape(bearing_vec(:,1:num_events_batch),1,[]));
rotated_vec_prev = rbv(:,1:3:end) + rbv(:,2:3:end) + rbv(:,3:3:end);
pm_prev = project_EquirectangularProjection(rotated_vec_prev, pano_width, pano_height);

% 3. Subtract the intensities at both map points
% Bilinear interpolation of intensities.
% Using one function call instead of two is more efficient
% Faster code than interp2.
% https://www.somesolvedproblems.com/2017/12/how-do-i-do-fast-bilinear-interpolation.html
F = griddedInterpolant({1:pano_height,1:pano_width}, map);
M_ = F([pm(2,:),pm_prev(2,:)], [pm(1,:),pm_prev(1,:)]);
contrast = M_(1:num_events_batch) - M_(1+num_events_batch:end); % M_cur - M_prev
contrast = contrast(:);

% Compute Jacobians
if (nargout > 1) && (nargin > 5)
    % Derivative of the map wrt the 2D point pm
    % Faster code than interp2.
    Fx = griddedInterpolant({1:pano_height,1:pano_width}, grad_map.x);
    M_deriv_x = Fx(pm(2,:), pm(1,:));
    Fy = griddedInterpolant({1:pano_height,1:pano_width}, grad_map.y);
    M_deriv_y = Fy(pm(2,:), pm(1,:));
    
    % Derivative of the 2D point pm wrt rotated_bvec
    Jac = (M_deriv_x' * [1 1 1]) .* (squeeze(dpm_d3D(1,:,:))') + ...
        (M_deriv_y' * [1 1 1]) .* (squeeze(dpm_d3D(2,:,:))');
    % This is a source of inaccuracies: derivatives of the map by finite
    % differences differ from derivatives using bilinear interpolation
    
    % Derivative of rotated_bvec wrt rotvec_pred (the state)
    % Gallego and Yezzi, JMIV 2014. Formula for rotvec in [0,pi]
    v = rotvec_pred;
    matrix_factor = (v*v' + (Rot'-eye(3))*v2skew(v)) / (norm(v).^2);
    Jac = cross(rotated_bvec',Jac,2)*matrix_factor'; % Chain rule
end
