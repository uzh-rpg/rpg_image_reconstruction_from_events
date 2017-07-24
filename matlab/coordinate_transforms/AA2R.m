function R = AA2R(AA)
% function R = AA2R(AA)
% Convert Axis Angle (AA) to Rotation Matrix (R)
% 
%  Input
%     AA: 1x4 axis angle rotation.
% 
% Ouput:
%     R: a 3x3 rotation matrix
%     
% Written by Garrick Orchard July 2017
% Based on:
% http://mathworld.wolfram.com/RodriguesRotationFormula.html

AA(1:3) = AA(1:3)./norm(AA(1:3)); %ensure a unit length vector for the axis

omega = [ 0,    -AA(3), AA(2);...
        AA(3),     0, -AA(1);...
       -AA(2), AA(1), 0];
   
theta = AA(4);

R = eye(3) + omega*sin(theta) + omega*omega*(1-cos(theta));