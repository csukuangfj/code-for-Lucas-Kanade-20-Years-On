function dW_dp = jacobian_a(nx, ny);
% JACOBIAN_A - Compute Jacobian for affine warp
%   DW_DP = JACOBIAN_A(WIDTH, HEIGHT)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: jacobian_a.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% See LK20-1 Equation (8)
jac_x = kron([0:nx - 1],ones(ny, 1));
jac_y = kron([0:ny - 1]',ones(1, nx));
jac_zero = zeros(ny, nx);
jac_one = ones(ny, nx);

dW_dp = [jac_x, jac_zero, jac_y, jac_zero, jac_one, jac_zero;
         jac_zero, jac_x, jac_zero, jac_y, jac_zero, jac_one];
