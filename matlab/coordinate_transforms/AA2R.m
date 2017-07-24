function R = AA2R(AA)
% function R = AA2R(AA)
% Convert Axis Angle (AA) to Rotation Matrix (R)
% 
%  Input
%     AA: 1x4 axis angle rotation
% 
% Ouput:
%     R: a 3x3 rotation matrix
%     
% Written by Garrick Orchard July 2017
% Based on:
% https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis.E2.80.93angle


u_x = [ 0,    -AA(3), AA(2);...
        AA(3),     0, -AA(1);...
       -AA(2), AA(1), 0];
   

R = cos(AA(4))*eye(3) + sin(AA(4))*u_x + (1-cos(AA(4)))*AA(1:3)'*AA(1:3);