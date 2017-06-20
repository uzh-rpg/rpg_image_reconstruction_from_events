% Given an intrinsic calibration, compute a look-up-table (i.e., map) of
% undistorted pixel locations.

clc

%% Yaml DAVIS
% Copied from the yaml file
% image_width = 240;
% image_height = 180;
% Kint = [201.55659084412122, 0.0, 140.73977188104507; 0.0, 202.75456850041195, 106.65922789946397; 0.0, 0.0, 1.0];
% dist_coeffs = [-0.3734645841479531, 0.15098967465602398, 0.0006110153203668791, -0.0026759416877948096, 0];

% DVS synthetic experiment
image_width = 128;
image_height = 128;
Kint = [91.4014729896821, 0.0, 64.0; 0.0, 91.4014729896821, 64.0; 0.0, 0.0, 1.0];
dist_coeffs = [0,0,0,0,0];

% Matlab's camera structure
cameraParams = cameraParameters('IntrinsicMatrix',Kint.', ...
    'RadialDistortion', dist_coeffs([1,2,5]),...
    'TangentialDistortion', dist_coeffs(3:4))


%% Integer coordinates (pixels)
x = 0:image_width-1;
y = 0:image_height-1;
[xx,yy] = meshgrid(x,y);
pixels_grid = [xx(:) yy(:)];
num_pixels = image_width * image_height;


%% Create look-up-table (i.e., map) of undistorted pixel locations
num_chunks = 100;
undist_pix = zeros(num_pixels,2);
num_pix_subset = floor( num_pixels / num_chunks );
for i = 1:num_chunks
    i  % verbose
    idx = (1 +(i-1)*num_pix_subset) : (i*num_pix_subset);
    
    % Call Matlab's function from Computer Vision toolbox
    undist_pix(idx,:) = undistortPoints(pixels_grid(idx,:),cameraParams);
end


%% Remove the intrinsics
undist_pix_calibrated = zeros(size(undist_pix));
for i = 1:size(undist_pix,1)
    calib_pt = inv(Kint) * [undist_pix(i,:), 1].';
    undist_pix_calibrated(i,:) = calib_pt(1:2);
end
    

%% Visualize undistorted pixels
image_corners = [0,0,image_width-1,image_width-1;
    0,image_height-1,image_height-1,0];

% Too dense
figure, 
plot(undist_pix(:,1),undist_pix(:,2),'.')
grid on
% plot image limits rectangle
hold on
plot(image_corners(1,[1:4,1]),image_corners(2,[1:4,1]),'r')


%% Visualize undistorted pixels (not so dense)
undist_pix_mat_x = reshape(undist_pix(:,1),image_height,image_width);
undist_pix_mat_y = reshape(undist_pix(:,2),image_height,image_width);
figure,
skip_count = 5;
xu = undist_pix_mat_x(1:skip_count:end,1:skip_count:end);
yu = undist_pix_mat_y(1:skip_count:end,1:skip_count:end);
plot(xu(:),yu(:),'.')
grid on
% plot image limits rectangle
hold on
plot(image_corners(1,[1:4,1]),image_corners(2,[1:4,1]),'r')


%% Save undistortion map to disk

% save DAVIS_32_undistorted_pixels undist_pix pixels_grid image_width image_height undist_pix_mat_x undist_pix_mat_y Kint dist_coeffs undist_pix_calibrated

% save DVS_synth_undistorted_pixels undist_pix pixels_grid image_width image_height undist_pix_mat_x undist_pix_mat_y Kint dist_coeffs undist_pix_calibrated
