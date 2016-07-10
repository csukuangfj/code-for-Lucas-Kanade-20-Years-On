function wimg = warp_h(img, p, dst)
% WARP_H - Projective warp the image
%   WIMG = WARP_H(IMG, P, DST)
%   Warp image IMG to WIMG. DST are the destination points, i.e. the corners
%   of the template image. P are the affine warp parameters that project
%   DST into IMG.
%
%   P = [p1, p4, p7        (i.e. transpose of what is in LK20-1)
%        p2, p5, p8   
%        p3, p6, 1]; 

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: warp_h.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<3 error('Not enough input arguments'); end

M = p;
M(1,1) = M(1,1) + 1;
M(2,2) = M(2,2) + 1;

% Use bilinear filtering to warp back to template
wimg = quadtobox_h(img, dst, M, 'bilinear');
