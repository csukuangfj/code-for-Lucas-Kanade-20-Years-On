function verb_plot_a(verb, warp_p, tmplt_pts, error_img)
% VERB_PLOT_A - Verbose fitting plot
%   VERB_PLOT_A(V, WARP_P, TMPLT_PTS, ERROR_IMG)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: verb_plot_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<2 error('Not enough input arguments'); end

% Scaled error image
set(verb.ih_error, 'Cdata', (error_img + 256) / 2);

% Current parameters
M = [warp_p; 0 0 1];
M(1,1) = M(1,1) + 1;
M(2,2) = M(2,2) + 1;
warp_pts =  M * [tmplt_pts; ones(1, size(tmplt_pts,2))];
set(verb.lh, 'Xdata', [warp_pts(1,:) warp_pts(1,1)], 'Ydata', [warp_pts(2,:) warp_pts(2,1)]);
drawnow
