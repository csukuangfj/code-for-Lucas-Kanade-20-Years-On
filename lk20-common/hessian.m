function H = hessian(VI_dW_dp, N_p, w)
% HESSIAN - Compute Hessian
%   H = HESSIAN(VI_DW_DP, N_P, W)

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: hessian.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

if nargin<3 error('Not enough input arguments'); end

H = zeros(N_p, N_p);
for i=1:N_p
	h1 = VI_dW_dp(:,((i-1)*w)+1:((i-1)*w)+w);
	for j=1:N_p
		h2 = VI_dW_dp(:,((j-1)*w)+1:((j-1)*w)+w);
		H(j, i) = sum(sum((h1 .* h2)));
	end
end
