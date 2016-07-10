function wimg = warp_a(img, p, dst)
% WARP_A - Affine warp the image
%   WIMG = WARP_A(IMG, P, DST)
%   Warp image IMG to WIMG. DST are the destination points, i.e. the corners
%   of the template image. P are the affine warp parameters that project
%   DST into IMG.
%
%   P = [p1, p3, p5
%        p2, p4, p6];

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: warp_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<3 error('Not enough input arguments'); end

% Convert affine warp parameters into 3 x 3 warp matrix
% NB affine parameterised as [1 + p1, p3, p5; p2, 1 + p4, p6]
M = [p; 0 0 1];
M(1,1) = M(1,1) + 1;
M(2,2) = M(2,2) + 1;

% Use bilinear filtering to warp image back to template
wimg = quadtobox(img, dst, M, 'bilinear');
