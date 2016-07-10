function fit = affine_ia(img, tmplt, p_init, n_iters, verbose, step_size)
% AFFINE_IA - Affine image alignment using inverse-additive algorithm
%   FIT = AFFINE_HB(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
%   Align the template image TMPLT to an example image IMG using an
%   affine warp initialised using P_INIT. Iterate for N_ITERS iterations.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p3, p5
%             p2, p4, p6];
%
%   This assumes greyscale images and rectangular templates.
%
%   c.f. Hager-Belhumeur

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: affine_ia.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

if nargin<5 verbose = 0; end
if nargin<4 error('Not enough input arguments'); end

% Common initialisation
init_a;

% Pre-computable things ---------------------------------------------------

% 3) Evaluate gradient of T
[nabla_Tx, nabla_Ty] = gradient(tmplt);

% 4) Evaluate Gamma_x
Gamma_x = jacobian_a(w, h);

% 5) Compute modified steepest descent images
VT_Gamma_x = sd_images(Gamma_x, nabla_Tx, nabla_Ty, N_p, h, w);

% 6) Compute modified Hessian and inverse
H_star = hessian(VT_Gamma_x, N_p, w);
H_star_inv = inv(H_star);


% Hager-Belhumeur, Inverse Additive Algorithm -----------------------------

for f=1:n_iters
	% 1) Compute warped image with current parameters
	IWxp = warp_a(img, warp_p, tmplt_pts);

	% 2) Compute error image - NB reversed
	error_img = IWxp - tmplt;
	
	% -- Save current fit parameters --
	fit(f).warp_p = warp_p;
	fit(f).rms_error = sqrt(mean(error_img(:) .^2));
	
	% -- Show fitting? --
	if verbose
		disp(['Inverse-Additive [',num2str(f-1),']: RMS = ',num2str(fit(f).rms_error)]);
		verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
	end
	
	% -- Really iteration 1 is the zeroth, ignore final computation --
	if (f == n_iters) break; end
	
	% 7) Compute steepest descent parameter updates
	sd_delta_p = sd_update(VT_Gamma_x, error_img, N_p, w);

	% 8) Compute gradient descent parameter updates
	delta_p_star = H_star_inv * sd_delta_p;

	% 9) Update warp parmaters
	warp_p = update_step(warp_p, delta_p_star);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function warp_p = update_step(warp_p, delta_p_star)
% Compute and apply the inverse additive update
% Note: parameterisation is 1+p1, 1+p4

% 9) Compute Sigma_p_inv
ss = warp_p(:,1:2);
ss(1,1) = ss(1,1) + 1;
ss(2,2) = ss(2,2) + 1;
Sigma_p_inv = kron(diag(ones(3,1)), ss);

% and update warp parmaters
delta_p = Sigma_p_inv * delta_p_star;
delta_p = reshape(delta_p, 2, 3);
warp_p = warp_p - delta_p;
