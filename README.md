# Image Reconstruction from an Event Camera

This repository contains code for brightness image reconstruction from a rotating event camera. 
For simplicity, we assume that the orientation of the camera is given, e.g., it is provided by a pose-tracking algorithm or by ground truth camera poses. 
The algorithm uses a per-pixel Extended Kalman Filter (EKF) approach to estimate the brightness image or gradient map that caused the events.

## Disclaimer and License

This code has been tested with MATLAB R2017a on Ubuntu 16.04.
This is research code, expect that it changes often and any fitness for a particular purpose is disclaimed.
The source code is released under a GNU General Public License (GPL).


## Instructions

Please run the file `matlab/test_image_reconstruction.m`.

This script reads a file of events and a file of camera rotations and produces a panoramic image with the reconstructed brightness that caused the events. See the example provided.

When running the script, a figure will emerge showing the evolution of the reconstructed brightness map as events are being processed:

![screenshot_image_reconstr_small](https://user-images.githubusercontent.com/8024432/27339188-6b586840-55d7-11e7-8999-fa8dba6399c1.png)

The reconstructed brightness image obtained after processing all events is the following: 
![reconstructed_image_log_scale](https://github.com/uzh-rpg/rpg_image_reconstruction_from_events/blob/master/data/synth1/sample_output/reconstructed_image_log_scale.jpg)

The above reconstructed map is displayed in logarithmic scale since the event camera known as Dynamic Vision Sensor (DVS) operates on logarithmic brightness.

### Some details on the EKF approach

**Two possible measurement functions** are provided for the EKF correction step: 
  - the *event rate* (the reciprocal of the time between two consecutive events within the same pixel), which gives an explicit EKF. See references [1] and [2] below.
  - the *brightness contrast* (the quantity thresholded by the event camera to generate events), which gives an implicit EKF. See reference [3] below.
  
A sample of the output produced by the algorithm can be found in [this folder](https://github.com/uzh-rpg/rpg_image_reconstruction_from_events/tree/master/data/synth1/sample_output).

  - EKF output: **estimated gradient map**, which is "integrated" using Poisson image reconstruction to yield the above **brightness map**. 
  The following image shows the magnitude and direction of the gradient map, combined: color represents direction, whereas saturation represents magnitude.
  ![mosaicing_grad_map_hsv](https://github.com/uzh-rpg/rpg_image_reconstruction_from_events/blob/master/data/synth1/sample_output/mosaicing_grad_map_hsv.png)

  - EKF output: trace of the **error covariance**. Points with smaller covariance (in red color) represent map points with a more confident estimation, due to a larger number of measurements.
  ![covariance_trace_colored_log_scale](https://github.com/uzh-rpg/rpg_image_reconstruction_from_events/blob/master/data/synth1/sample_output/covariance_trace_colored_log_scale.jpg)


## Publications

If you use this code in an academic context, please cite the following references:

  1. H. Kim, A. Handa, R. Benosman, S.-H. Ieng, A.J. Davison, 
  [*Simultaneous Mosaicing and Tracking with an Event Camera*](http://www.bmva.org/bmvc/2014/papers/paper066/).
  British Machine Vision Conference, 2014.

  2. H. Rebecq, T. Horstschaefer, G. Gallego, D. Scaramuzza, 
  [*EVO: A Geometric Approach to Event-based 6-DOF Parallel Tracking and Mapping in Real-time*](http://rpg.ifi.uzh.ch/docs/RAL16_EVO.pdf). 
  IEEE Robotics and Automation Letters (RA-L), Vol. 2, Issue 2, pp. 593-600, Apr. 2017.

  3. G. Gallego, C. Forster, E. Mueggler, D. Scaramuzza, 
  [*Event-based Camera Pose Tracking using a Generative Event Model*](https://arxiv.org/pdf/1510.01972v1).
  arXiv:1510.01972, 2015.
