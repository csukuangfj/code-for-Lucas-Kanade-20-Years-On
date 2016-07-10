function fit = affine_ic_lm(img, tmplt, p_init, n_iters, verbose, step_size)
% AFFINE_IC_LM - Affine image alignment using inverse-compositional
% Levenberg-Marquardt algoritym
%   FIT = AFFINE_IC_LM(IMG, TMPLT, P_INIT, N_ITERS, VERBOSE)
%   Align the template image TMPLT to an example image IMG using an
%   affine warp initialised using P_INIT. Iterate for N_ITERS iterations.
%   To display the fit graphically set VERBOSE non-zero.
%
%   p_init = [p1, p3, p5
%             p2, p4, p6];
%
%   This assumes greyscale images and rectangular templates.
%
%   c.f. Baker-Matthews

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: affine_ic_lm.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

if nargin<5 verbose = 0; end
if nargin<4 error('Not enough input arguments'); end

% Common initialisation
init_a;

% Pre-computable things ---------------------------------------------------

% 0) Levenburg Marquardt
delta = 0.001;

% 1) Compute warped image with current parameters
IWxp = warp_a(img, warp_p, tmplt_pts);

% 2) Compute error image - NB reversed
error_img = IWxp - tmplt;
e = sqrt(mean(error_img(:) .^2));

% 3) Evaluate gradient of T
[nabla_Tx nabla_Ty] = gradient(tmplt);

% 4) Evaluate Jacobian - constant for affine warps
dW_dp = jacobian_a(w, h);

% 5) Compute steepest descent images, VT_dW_dp
VT_dW_dp = sd_images(dW_dp, nabla_Tx, nabla_Ty, N_p, h, w);
	
% 6) Precompute part of Hessian
H = hessian(VT_dW_dp, N_p, w);

% LM Extra Hessian bit
LM = zeros(N_p, N_p);
for i=1:N_p
	h1 = VT_dW_dp(:,((i-1)*w)+1:((i-1)*w)+w);
	LM(i, i) = sum(sum(sum((h1 .* h1))));
end

% Zeroth iteration
f = 1;
fit(f).warp_p = warp_p;
fit(f).rms_error = e;

% -- Show fitting? --
if verbose
	disp(['Inverse-Compositional LM [',num2str(f-1),']: RMS = ',num2str(fit(f).rms_error)]);
	verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
end

% Baker-Matthews, Inverse Compositional Algorithm -------------------------

for f=2:n_iters
	% 7) Compute steepest descent parameter updates
	sd_delta_p = sd_update(VT_dW_dp, error_img, N_p, w);
	
	% 8) Compute LM Hessian
	H_lm = H + (delta * LM);
	H_inv = inv(H_lm);

	% 8) Compute gradient descent parameter updates
	delta_p = H_inv * sd_delta_p;
	
	% 9) Update warp parmaters
	warp_p_lm = update_step(warp_p, delta_p);
	
	% 1) Compute warped image with current parameters
	IWxp = warp_a(img, warp_p_lm, tmplt_pts);

	% 2) Compute error image - NB reversed
	error_img_lm = IWxp - tmplt;
	e_lm = sqrt(mean(error_img_lm(:) .^2));
	
	% 10) Levenburg-Marquardt
	if e < e_lm
		% Bad step, do nothing
		delta = delta * 10;
	else
		% Good step, update
		delta = delta / 10;
		error_img = error_img_lm;
		e = e_lm;
		warp_p = warp_p_lm;
	end
	
	% -- Save current fit parameters --
	fit(f).warp_p = warp_p;
	fit(f).rms_error = e;
	
	% -- Show fitting? --
	if verbose
		disp(['Inverse-Compositional LM [',num2str(f-1),']: RMS = ',num2str(fit(f).rms_error)]);
		verb_plot_a(verb_info, warp_p, tmplt_pts, error_img);
	end
	
	% -- Really iteration 1 is the zeroth, ignore final computation --
	if (f == n_iters) break; end	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function warp_p = update_step(warp_p, delta_p)
% Compute and apply the update

delta_p = reshape(delta_p, 2, 3);
	
% Convert affine notation into usual Matrix form - NB transposed
delta_M = [delta_p; 0 0 1];	
delta_M(1,1) = delta_M(1,1) + 1;
delta_M(2,2) = delta_M(2,2) + 1;

% Invert compositional warp
delta_M = inv(delta_M);

% Current warp
warp_M = [warp_p; 0 0 1];	
warp_M(1,1) = warp_M(1,1) + 1;
warp_M(2,2) = warp_M(2,2) + 1;

% Compose
comp_M = warp_M * delta_M;	
warp_p = comp_M(1:2,:);
warp_p(1,1) = warp_p(1,1) - 1;
warp_p(2,2) = warp_p(2,2) - 1;
