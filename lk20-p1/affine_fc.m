function fit = affine_fc(img, tmplt, p_init, n_iters, verbose, step_size)
% AFFINE_FC - Affine image alignment using forwards-compositional algorithm
%   FIT = AFFINE_FC(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
%   Align the template image TMPLT to an example image IMG using an
%   affine warp initialised using P_INIT. Iterate for N_ITERS iterations.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p3, p5
%             p2, p4, p6];
%
%   This assumes greyscale images and rectangular templates.
%
%   c.f. Shum-Szeliski

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: affine_fc.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

if nargin<5 verbose = 0; end
if nargin<4 error('Not enough input arguments'); end

% Common initialisation
init_a;

% Pre-computable things ---------------------------------------------------

% 4) Evaluate Jacobian - constant for affine warps
dW_dp = jacobian_a(w, h);


% Shum-Szeliski, Forwards Compositional Algorithm -------------------------

for f=1:n_iters
	% 1) Compute warped image with current parameters
	IWxp = warp_a(img, warp_p, tmplt_pts);

	% 2) Compute error image
	error_img = tmplt - IWxp;
	
	% -- Save current fit parameters --
	fit(f).warp_p = warp_p;
	fit(f).rms_error = sqrt(mean(error_img(:) .^2));
	
	% -- Show fitting? --
	if verbose
		disp(['Forwards-Compositional [',num2str(f-1),']: RMS = ',num2str(fit(f).rms_error)]);
		verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
	end
	
	% -- Really iteration 1 is the zeroth, ignore final computation --
	if (f == n_iters) break; end

	% 3) Evaluate gradient at I(W(x;p))
	[nabla_Ix nabla_Iy] = gradient(IWxp);
	
	% 5) Compute steepest descent images, VI_dW_dp
	VI_dW_dp = sd_images(dW_dp, nabla_Ix, nabla_Iy, N_p, h, w);
	
	% 6) Compute Hessian and inverse
	H = hessian(VI_dW_dp, N_p, w);
	H_inv = inv(H);
	
	% 7) Compute steepest descent parameter updates
	sd_delta_p = sd_update(VI_dW_dp, error_img, N_p, w);

	% 8) Compute gradient descent parameter updates
	delta_p = H_inv * sd_delta_p;
	
	% 9) Update warp parmaters - compose W(x;p) and W(x;delta_p)
	warp_p = update_step(warp_p, delta_p);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function warp_p = update_step(warp_p, delta_p)
% Compute and apply the compositional update
% Note: parameterisation is 1+p1, 1+p4

delta_p = reshape(delta_p, 2, 3);
	
% Convert affine parameters into usual Matrix form
delta_M = [delta_p; 0 0 1];	
delta_M(1,1) = delta_M(1,1) + 1;
delta_M(2,2) = delta_M(2,2) + 1;

% Current warp
warp_M = [warp_p; 0 0 1];	
warp_M(1,1) = warp_M(1,1) + 1;
warp_M(2,2) = warp_M(2,2) + 1;

% Compose
comp_M = warp_M * delta_M;	
warp_p = comp_M(1:2,:);
warp_p(1,1) = warp_p(1,1) - 1;
warp_p(2,2) = warp_p(2,2) - 1;
