function dW_dp = jacobian_h(nx, ny, warp_p);
% JACOBIAN_H - Compute Jacobian for projective warp
%   DW_DP = JACOBIAN_H(WIDTH, HEIGHT, P)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: jacobian_h.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% Easy bits
jac_x = kron([0:nx - 1],ones(ny, 1));
jac_y = kron([0:ny - 1]',ones(1, nx));
jac_zero = zeros(ny, nx);
jac_one = ones(ny, nx);

% Complicated bits are just homography of all image coordinates
xy = [repmat([1:ny]',nx,1) kron([1:nx]', ones(ny,1))];
xy = [xy ones(length(xy),1)]';
M = warp_p;
M(1,1) = M(1,1) + 1;
M(2,2) = M(2,2) + 1;
uv = M * xy;
uvc = uv ./ repmat(uv(3,:),3,1);

u_x = reshape(uvc(1,:),nx,ny)';
u_y = reshape(uvc(2,:),nx,ny)';
v = reshape(uv(3,:),nx,ny)';

% Divide each jacobian image by v
iv = 1 ./ v;
jac_x = iv .* jac_x;
jac_y = iv .* jac_y;
jac_one = iv .* jac_one;

dW_dp = [jac_x, jac_zero, -jac_x .* u_x, jac_y, jac_zero, -jac_y .* u_x, jac_one, jac_zero;
         jac_zero, jac_x, -jac_x .* u_y, jac_zero, jac_y, -jac_y .* u_y, jac_zero, jac_one];
