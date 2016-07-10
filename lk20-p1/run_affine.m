function run_affine(data_name)
% Run an affine perturbation test
%
% e.g. run_affine('takeo');
%
% You should edit this file to define experiment parameters!

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: run_affine.m,v 1.1.1.1 2003/08/20 03:07:35 iainm Exp $

if nargin<1 error('Not enough input arguments'); end

% List of algorithms to run
alg_list = {'affine_fa' 'affine_fc' 'affine_ia' 'affine_ic'};

% Test parameters
verbose = 1;					% Show fitting?
n_iters = 15;					% Number of gradient descent iterations
n_tests = 10;					% Number of convergence tests
n_freq_tests = 100;				% Number of frequency of convergence tests
max_spatial_error = 1;			% Max location error for deciding convergence
opt_step_size = 0;				% For diagonal methods, compute best step size?	

img_pix_sig = [0 4 8 16 32];	% Noise power sigmas for image
tmp_pix_sig = [0 4 8 16 32];	% Noise power sigmas for template

all_spc_sig = [1:0.5:10];		% All spatial sigmas
cnv_spc_sig = [2 6 10];			% Those that should run convergence tests

% Should not need to modify anything below --------------------------------

% tdata - the image and initial template
load(['data/', data_name])

% pt_offset - precomputed random point offsets
load('data/affine_pt_offset');

% Check output directory exists
if ~exist('results/affine', 'dir')
	if ~mkdir('results/affine')
		error('Unable to create results directory');
	end
end

% Run tests
for i=1:length(img_pix_sig)
	image_pixel_sigma = img_pix_sig(i);
	
	for t=1:length(tmp_pix_sig)
		tmplt_pixel_sigma = tmp_pix_sig(t);
		
		for s=1:length(all_spc_sig)
			spatial_sigma = all_spc_sig(s);
			
			% Convergence and frequence of convergence test
			if ~isempty(find(cnv_spc_sig == spatial_sigma))
				disp(['CONV - Spatial: ', num2str(spatial_sigma), ...
					 ', Image: ', num2str(image_pixel_sigma), ...
					 ', Template: ', num2str(tmplt_pixel_sigma)]);
			
				results = test_affine(tdata, pt_offset, alg_list, n_iters, n_tests, n_freq_tests, spatial_sigma, image_pixel_sigma, tmplt_pixel_sigma, max_spatial_error, verbose, opt_step_size);
			
			% Just frequency of convergence
			else
				disp(['DIVG - Spatial: ', num2str(spatial_sigma), ...
					 ', Image: ', num2str(image_pixel_sigma), ...
					 ', Template: ', num2str(tmplt_pixel_sigma)]);
			
				results = test_affine(tdata, pt_offset, alg_list, n_iters, 0, n_freq_tests, spatial_sigma, image_pixel_sigma, tmplt_pixel_sigma, max_spatial_error, verbose, opt_step_size);
			end
			
			% Save results
			disp('********************************************');
			string = ['save results/affine/', data_name, '_I', num2str(image_pixel_sigma),...
				'_T', num2str(tmplt_pixel_sigma), ...
				'_S', num2str(spatial_sigma, '%.1f'), ...
				'.mat  results max_spatial_error opt_step_size'];
			disp(string);
			eval(string);
			disp('********************************************');
		end
	end
end
