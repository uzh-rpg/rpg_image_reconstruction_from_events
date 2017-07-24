function AA = R2AA(R)
% function AA = R2AA(R)
% Convert Rotation Matrix (R) to Axis Angle (AA)
% 
% Input
%     R: a 3x3 rotation matrix
% 
% Ouput:
%     AA: a 1x4 axis angle. The axis is not normalized
%     
% Written by Garrick Orchard July 2017
% Based on:
% https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis.E2.80.93angle

AA = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
AA = AA./norm(AA); %ensure a unit length axis direction
AA(4) = acos((trace(R)-1)/2);
