function results = test_homo(tdata, pt_offsets, alg_list, n_iters, n_tests, n_freq_tests, spatial_sigma, image_pixel_sigma, tmplt_pixel_sigma, max_spatial_error, verbose)
% TEST_HOMO - Test homography algorithms
%
% See also: run_homo
%
% tdata has two fields:
%   tdata.img       unit8 greyscale image
%   tdata.tmplt     [x_start y_start x_end y_end] template rectangle corners
%                   matlab coordinates (1, 1) is top-left, start < end.
%
% pt_offsets is a N x 8 matrix of (x1, y1, x2, y2, x3, y3, x4, y4) deltas for 
% each corner. Zero mean, unit variance, 
% e.g. pt_offsets = randn(N, 8);
%
% alg_list is a cellstr list of algorithms to run, e.g.:
% alg_list = {'homo_fa' 'homo_fc' 'homo_ic'};
%
% The "template" is created by cutting out a distorted version of tdata.tmplt.
% The "target" is the image defined by tdata.tmplt.

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: test_homo.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

% Check args
if nargin<6 error('Invalid args'); end

% Target corner points (i.e. correct answer)
target_pts = [tdata.tmplt(1), tdata.tmplt(2);
              tdata.tmplt(1), tdata.tmplt(4);
			  tdata.tmplt(3), tdata.tmplt(4);
			  tdata.tmplt(3), tdata.tmplt(2)]';

% Template image dimensions
template_nx = tdata.tmplt(3) - tdata.tmplt(1) + 1;
template_ny = tdata.tmplt(4) - tdata.tmplt(2) + 1;

% Template test image corner points				   
template_pts = [1, 1; 
			    1, template_ny;
			    template_nx, template_ny;
			    template_nx, 1]';

% Initialise at unperturbed translation?
p_init = zeros(2,3);
p_init(1, 3) = tdata.tmplt(1) - 1;
p_init(2, 3) = tdata.tmplt(2) - 1;
p_init(3, 3) = 1;

% Translate by 0.5 pixels to avoid identity warp. Warping introduces a little
% smoothing and this avoids the case where the first iteration for a forwards
% algorithm is on the "unsmoothed" unwarped image
p_init(1, 3) = p_init(1, 3) + 0.5;
p_init(2, 3) = p_init(2, 3) + 0.5;

% Scale point offsets to have required sigma
pt_offsets = pt_offsets * spatial_sigma;

% Need image to be doubles
tdata.img = double(tdata.img);

% Space for results
results = [];

% Test counters
go = 1; offset_idx = 1; all_alg_converged = 0;

% Convergence counters in field 1
for l=1:length(alg_list)
	results = setfield(results, {1}, alg_list{l}, {1}, 'n_diverge', 0);
	results = setfield(results, {1}, alg_list{l}, {1}, 'n_converge', 0);
end

% Index of successfull convergence tests
results.all_converged_idx = [];

% Might not be doing any convergence tests at all
if n_tests > 0
	cnvg_testing = 1;
else
	cnvg_testing = 0;
	cnvg_result = [];
end

% Run
while go
	if cnvg_testing
		disp(['Convergence Test: ', num2str(all_alg_converged + 1), ' (of total ', num2str(offset_idx), ')']);
	else
		disp(['Divergence Test: ', num2str(offset_idx)]);
	end 
	
	% Test points: apply current point offset to target points
	test_pts = target_pts + reshape(pt_offsets(offset_idx,:), 2, 4);
		
	% Solve for homography to create test template image
	uv = template_pts';
	xy = test_pts';
	A = [uv, ones(4,1), zeros(4,3), -uv .* repmat(xy(:,1),1,2);
	     zeros(4,3), uv, ones(4,1), -uv .* repmat(xy(:,2),1,2)];
	M = A \ xy(:);
	M = [M; 1];  M = reshape(M, 3, 3);  M = M';
	
	% Get source and destination boxes
	dst = template_pts;
	src = M * [template_pts; ones(1,4)];
	src = src ./ repmat(src(3,:),3,1);
	
	% Warp image to get test template image
	tmplt_img = quadtobox_h(tdata.img, dst, M, 'bilinear');
	
	% Add noise to template
	if tmplt_pixel_sigma > 0
		tmplt_img = tmplt_img + (randn(size(tmplt_img)) * tmplt_pixel_sigma);
	end
	
	% Add noise to image
	if image_pixel_sigma > 0
		noisy_img = tdata.img + (randn(size(tdata.img)) * image_pixel_sigma);
	else
		noisy_img = tdata.img;
	end
	
	% Initial error in points. This is not quite sqrt(mean(pt_offset(offset_idx,:) .^ 2)) due to p_init
	rms_pt_init = ComputePointError(test_pts, template_pts, p_init);
	
	% Run each algorithm
	for l=1:length(alg_list)
		string = ['tic; fit = ', alg_list{l}, '(noisy_img, tmplt_img, p_init, n_iters, verbose); t = toc;'];
		eval(string);
		
		% Evaluate point spatial error for each iteration for convergence tests
		if cnvg_testing
			rms_pt_err = zeros(length(fit), 1);
			for f=1:length(fit)
				rms_pt_error(f) = ComputePointError(test_pts, template_pts, fit(f).warp_p);
			end
					
		% Only need final spatial rms for divergence tests. Don't need fitting results
		else
			rms_pt_error = ComputePointError(test_pts, template_pts, fit(end).warp_p);
			fit = [];
		end	
		
		% Save spatial errors
		results = setfield(results, {1}, alg_list{l}, {offset_idx}, 'rms_pt_error', rms_pt_error);
			
		% Save full fitting results
		results = setfield(results, {1}, alg_list{l}, {offset_idx}, 'fit', fit);
			
		% Save fitting time
		results = setfield(results, {1}, alg_list{l}, {offset_idx}, 'time', t);
	end
	
	% Evaluate final spatial errors for all algorithms
	disp('----------------------------------------------------');
	disp(['Initial spatial rms = ',num2str(rms_pt_init)]);
	
	all_cnvg_check = 1;
	
	for l=1:length(alg_list)
		final_rms_pt_error = getfield(results, {1}, alg_list{l}, {offset_idx}, 'rms_pt_error');
		final_rms_pt_error = final_rms_pt_error(end);

		string = [alg_list{l}, ': final spatial rms = ', num2str(final_rms_pt_error)];
		
		% Convergence test
		if final_rms_pt_error > max_spatial_error
			string = [string, '   **DIVERGED**'];
			
			% Update divergence counter in first field
			n_diverge = getfield(results, {1}, alg_list{l}, {1}, 'n_diverge');
			n_diverge = n_diverge + 1;
			results = setfield(results, {1}, alg_list{l}, {1}, 'n_diverge', n_diverge);
			
			% One or more algorithms diverged
			all_cnvg_check = 0;
			
		else
			% Update convergence counter in first field
			n_converge = getfield(results, {1}, alg_list{l}, {1}, 'n_converge');
			n_converge = n_converge + 1;
			results = setfield(results, {1}, alg_list{l}, {1}, 'n_converge', n_converge);
		end
		
		disp(string);
	end
	disp('----------------------------------------------------');
	
	% Everything converged?
	if all_cnvg_check
		all_alg_converged = all_alg_converged + 1;
		
		% Update index of fully converged test results
		results.all_converged_idx = [results.all_converged_idx; offset_idx];
	end
	
	% Finished convergence testing?
	if all_alg_converged >= n_tests
		cnvg_testing = 0;
		
		% Extra tests for divergence?
		if n_freq_tests > n_tests
			if offset_idx >= n_freq_tests
				go = 0;
			end
			
		% Or just stop
		else
			go = 0;
		end
	end
	
	% Always increment index into point offsets
	offset_idx = offset_idx + 1;
	if offset_idx > size(pt_offsets, 1)
		disp('Ran out of tests!');
		return;
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rms_pt_error = ComputePointError(test_pts, template_pts, warp_p)
% Compute point rms error

% Warp for this iteration
M = warp_p;
M(1,1) = M(1,1) + 1; 
M(2,2) = M(2,2) + 1;

% Homography points
iteration_pts = M * [template_pts; ones(1, 4)];
iteration_pts = iteration_pts ./ repmat(iteration_pts(3,:),3,1);

% Error in points
diff_pts = test_pts - iteration_pts(1:2,:);
diff_pts = diff_pts(:);
rms_pt_error = sqrt(mean(diff_pts .^ 2));
