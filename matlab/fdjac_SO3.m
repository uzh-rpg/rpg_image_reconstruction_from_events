% fdjac_SO3  Computes forward-difference approximation to Jacobian.
%
% On input,
%   - x is the rotation matrix at which the Jacobian is to be evaluated
%   - g is a user-supplied routine that returns a vector function at x
%   - varargin: optional parameters passed directly to the function g
%           g(x,varargin{:});
%
% On output
%   - df is the Jacobian matrix 
%      (as many rows as dependent variables, returned by the function g,
%       and with 3 columns: dim of elements in Lie Algebra so(3))


function df = fdjac_SO3(x_rot_op,g,varargin)

assert(all(size(x_rot_op)==3)); % 3x3, rotation mat

fx = feval(g,x_rot_op,varargin{:});

df = zeros(numel(fx(:)),3);
Id = eye(3);
h = 1e-8;
for j=1:3 % DOFS in a rotation
    x_rot = expm_v_SO3(h*Id(:,j)) * x_rot_op; % perturbation in Lie Theory
    fh = feval(g,x_rot,varargin{:});
    df(:,j) = (fh-fx)/h;
end
