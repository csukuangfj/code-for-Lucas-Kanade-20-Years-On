function wimg = quadtobox_homo(img, dst, M, f_type)
% QUADTOBOX_HOMO - Warp contents of a quadrilateral to a rectangle
%   W = QUADTOBOX_HOMO(IMG, DST, M, F_TYPE)
%
%   Destination corners, DST = [x1, x2, x3, x4; y1, y2, y3, y4];
%   These are the (assumed rectangular) corner points of the template
%   image and must be integers.
%
%   The projection matrix, M transforms DST into the image IMG 
%   [u' v' w] (image) = M * [x y 1] (DST - template);
%   u = u' / w; v = v' / w;
%
%   Matlab style indices, i.e. start from 1.

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: quadtobox_h.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% Check args
if nargin<4 f_type = 'bilinear'; end
if nargin<3 error('Invalid args'); end

% Dimensions of destination image - integers, assume rectangle
minv = min(dst');
maxv = max(dst');

% Get all points in destination to sample
[xg yg] = meshgrid(1:maxv(1), 1:maxv(2));
xy = [reshape(xg, prod(size(xg)), 1)'; reshape(yg, prod(size(yg)), 1)'];
xy = [xy; ones(1,size(xy,2))];

% Transform into source
uv = M * xy;

% Divide for homography
uv = uv ./ repmat(uv(3,:),3,1);

% Remove homogeneous
uv = uv(1:2,:)';

% Sample
xi = reshape(uv(:,1),maxv(1),maxv(2));
yi = reshape(uv(:,2),maxv(1),maxv(2));
wimg = interp2(img, xi, yi, f_type);

% Check for NaN background pixels - replace them with a background of 0
idx = find(isnan(wimg));
if ~isempty(idx)
	wimg(idx) = 0;
end
