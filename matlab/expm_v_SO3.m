function R = expm_v_SO3(v)
assert(numel(v)==3);
% Based on Rodrigues formula (a 3x3 matrix)
angl = norm(v); % angle
a = v / angl;   % unit axis
c = cos(angl);
s = sin(angl);
R = c*eye(3) + ((1-c)*a)*a' + s*v2skew(a);
