function fit = homo_ic(img, tmplt, p_init, n_iters, verbose)
% HOMO_IC - Homography image alignment using inverse-compositional algorithm
%   FIT = HOMO_IC(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
%   Align the template image TMPLT to an example image IMG using a
%   projective warp initialised using P_INIT. Iterate for N_ITERS iterations.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p4, p7     paper_equivalent = [p1, p3, p5
%             p2, p5, p8                         p2, p4, p6
%             p3, p6, 1];                        p7, p8, 1];
%
%   This assumes greyscale images and rectangular templates.
%
%   c.f. Baker-Matthews

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: homo_ic.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

if nargin<5 verbose = 0; end
if nargin<4 error('Not enough input arguments'); end

% Common initialisation
init_h;

% Pre-computable things ---------------------------------------------------

% 3) Evaluate gradient of T
[nabla_Tx nabla_Ty] = gradient(tmplt);

% 4) Evaluate Jacobian at (x;0)
p_0 = zeros(size(warp_p));  p_0(3,3) = 1;
dW_dp = jacobian_h(w, h, p_0);

% 5) Compute steepest descent images, VT_dW_dp
VT_dW_dp = sd_images(dW_dp, nabla_Tx, nabla_Ty, N_p, h, w);
	
% 6) Compute Hessian and inverse
H = hessian(VT_dW_dp, N_p, w);
H_inv = inv(H);


% Baker-Matthews, Inverse Compositional Algorithm -------------------------

for f=1:n_iters
	% 1) Compute warped image with current parameters
	IWxp = warp_h(img, warp_p, tmplt_pts);

	% 2) Compute error image - NB reversed
	error_img = IWxp - tmplt;
	
	% -- Save current fit parameters --
	fit(f).warp_p = warp_p;
	fit(f).rms_error = sqrt(mean(error_img(:) .^2));
	
	% -- Show fitting? --
	if verbose
		disp(['Inverse-Compositional [',num2str(f-1),']: RMS = ',num2str(fit(f).rms_error)]);
		verb_plot_h(verb_info, warp_p, tmplt_pts, error_img);
	end
	
	% -- Really iteration 1 is the zeroth, ignore final computation --
	if (f == n_iters) break; end

	% 7) Compute steepest descent parameter updates
	sd_delta_p = sd_update(VT_dW_dp, error_img, N_p, w);

	% 8) Compute gradient descent parameter updates
	delta_p = H_inv * sd_delta_p;
	
	% 9) Update warp parmaters
	warp_p = update_step(warp_p, delta_p);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function warp_p = update_step(warp_p, delta_p)
% Compute and apply the update

% Update warp
delta_M = [delta_p; 1];
delta_M = reshape(delta_M, 3, 3);
delta_M(1,1) = delta_M(1,1) + 1;
delta_M(2,2) = delta_M(2,2) + 1;

% Invert compositional warp
delta_M = inv(delta_M);

% Current warp
warp_M = warp_p;
warp_M(1,1) = warp_M(1,1) + 1;
warp_M(2,2) = warp_M(2,2) + 1;

% Compose
comp_M = warp_M * delta_M;

% Correct for parameterisation
comp_M(1,1) = comp_M(1,1) - 1;
comp_M(2,2) = comp_M(2,2) - 1;

warp_p = comp_M;
