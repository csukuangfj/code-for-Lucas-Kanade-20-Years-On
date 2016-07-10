% init_h.m
% Common initialisation things for all homography algorithms

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: init_h.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% Need to process image data, must be real
if ~isa(img, 'double')
	img = double(img);
end

% Projective warp
N_p = 8;
if prod(size(p_init)) ~= 9
	error('Number of warp parameters incorrect');
end

% Initial warp parameters
warp_p = p_init;

% Template size
h = size(tmplt, 1);
w = size(tmplt, 2);

% Template verticies, rectangular [minX minY; minX maxY; maxX maxY; maxX minY]
tmplt_pts = [1 1; 1 h; w h; w 1]';

% Verbose display of fitting?
if verbose
	verb_info = verb_init_h(img, tmplt, tmplt_pts, warp_p);
end
