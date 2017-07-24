
function [rec] = poisson_solver_function_neumann(gx,gy)

% Least squares solution
% Poisson Reconstruction Using Neumann boundary conditions
% Input gx and gy
% gradients at boundary are assumed to be zero
% Output : reconstruction
% Author: Amit Agrawal, 2004
% 
% Obtained online from (July 2017):
% http://www.cs.cmu.edu/~ILIM/projects/IM/aagrawal/software.html
% 
% Permission to use, copy and modify this software and its documentation without fee for educational, research and non-profit purposes, is hereby granted, provided that the above copyright notice and the following three paragraphs appear in all copies.
% 
% To request permission to incorporate this software into commercial products contact: Vice President of Marketing and Business Development; Mitsubishi Electric Research Laboratories (MERL), 201 Broadway, Cambridge, MA 02139
% 
% In no event shall MERL be liable to any party for direct, indirect, special, incidental, or consequential damages, including lost profits, arising out of the use of this software and its documentation, even if MERL has been advised of the possibility of such damages.
% 
% MERL specifically disclaims any warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. the software provided hereunder is on an "as is" basis, and MERL has no obligations to provide maintenance, support, updates, enhancements or modifications.
% 

% disp('=======================================================')
% disp('Solving Poisson Equation using Neumann');


[H,W] = size(gx);

gx(:,end) = 0;
gy(end,:) = 0;

% pad
gx = padarray(gx,[1 1],0,'both');
gy = padarray(gy,[1 1],0,'both');


gxx = zeros(size(gx)); 
gyy = zeros(size(gx)); 
f = zeros(size(gx)); 

j = 1:H+1;
k = 1:W+1;

% Laplacian
gyy(j+1,k) = gy(j+1,k) - gy(j,k);
gxx(j,k+1) = gx(j,k+1) - gx(j,k);
f = gxx + gyy;
f = f(2:end-1,2:end-1);
clear j k gxx gyy gyyd gxxd gx gy


%compute cosine transform
fcos = dct2(f);
clear f

%compute eigen values
[x,y] = meshgrid(0:W-1,0:H-1);
denom = (2*cos(pi*x/(W))-2) + (2*cos(pi*y/(H)) - 2);
clear x y

%divide. 1st element of denom will be zero and 1st element of fcos and
%after division should also be zero; so divided rest of the elements
fcos(2:end) = fcos(2:end)./denom(2:end);
clear denom

% compute inverse dct2
rec = idct2(fcos);
