function z = frankotchellappa(dzdx,dzdy)
% 
% Frankt-Chellappa Algrotihm
% Input gx and gy
% Output : reconstruction
% Author: Amit Agrawal, 2005
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
% disp('Solving Using Frankot Chellappa Algorithm');

[rows,cols] = size(dzdx);

[wx, wy] = meshgrid(-pi/2:pi/(cols-1):pi/2,-pi/2:pi/(rows-1):pi/2);

% Quadrant shift to put zero frequency at the appropriate edge
wx = ifftshift(wx); wy = ifftshift(wy);

DZDX = fft2(dzdx);   % Fourier transforms of gradients
DZDY = fft2(dzdy);

% Integrate in the frequency domain by phase shifting by pi/2 and
% weighting the Fourier coefficients by their frequencies in x and y and
% then dividing by the squared frequency.  eps is added to the
% denominator to avoid division by 0.
j = sqrt(-1);

% dd = wx.^2 + wy.^2;
Z = (-j*wx.*DZDX -j*wy.*DZDY)./(wx.^2 + wy.^2 + eps); 


z = real(ifft2(Z));  % Reconstruction
z = z - min(z(:));
z = z/2;
