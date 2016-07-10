function verb = verb_init_h(img, tmplt, tmplt_pts, warp_p)
% VERB_INIT_H - Initialise verbose plot
%   V = VERB_INIT_H(IMG, TMPLT, TMPLT_PTS, WARP_P)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: verb_init_h.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<4 error('Not enough input arguments'); end

% Init figure
clf;
set(gcf,'DoubleBuffer','on');
colormap(gray(256));

% Input image
subplot(2,2,1);
image(img);
axis('image');
title('Image');

M = warp_p;
M(1,1) = M(1,1) + 1;
M(2,2) = M(2,2) + 1;
warp_pts = M * [tmplt_pts; ones(1, size(tmplt_pts,2))];
warp_pts = warp_pts ./ repmat(warp_pts(3,:), 3, 1);

hold on
verb.lh = plot([warp_pts(1,:) warp_pts(1,1)], [warp_pts(2,:) warp_pts(2,1)], 'g-');
hold off

% Template image
subplot(2,2,2);
image(tmplt);
axis('image');
title('Template');

% Error image
subplot(2,2,4);
verb.ih_error = image(zeros(size(tmplt)));
axis('image');
title('Error');
drawnow;
