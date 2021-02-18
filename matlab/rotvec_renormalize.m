function rotvec = rotvec_renormalize(rotvec)
% rotvec_renormalize Normalize rotation vector to the ball of radius pi ?

ax = rotvec;
angle = norm(ax);
if (angle > 1e-8)
    ax = ax / angle;
    axang = rotm2axang(axang2rotm([ax; angle]'));
    rotvec = axang(1:3)' * axang(4);
end
