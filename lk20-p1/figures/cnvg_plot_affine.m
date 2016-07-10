function cnvg_plot_affine(results)
% CNVG_PLOT_AFFINE(RESULTS)
% Plot convergence test data in struct RESULTS

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: cnvg_plot_affine.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% Load LK20 line styles
lk20style;
set(gca,'FontSize',font_size);

% Make assumptions about data!
lstyle = {ls.fa, ls.fc, ls.ia, ls.ic};
lwidth = [lw.fa, lw.fc, lw.ia, lw.ic];
msize = [ms.fa, ms.fc, ms.ia, ms.ic];
legend_string = {'Forwards Additive', 'Forwards Compositional' ,'Inverse Additive', 'Inverse Compositional'};
fnames = fieldnames(results);
cnvg_idx = results.all_converged_idx;
alg_idx = 1:length(fnames);
idx = strmatch('all_converged_idx', fnames);
alg_idx(idx) = [];

% May have done a divergence test too, in which case extra results
if length(cnvg_idx) > 100
	cnvg_idx = cnvg_idx(1:100);
end
if length(cnvg_idx) < 100
	disp('WARNING: Less than 100 converged');
end

rms_pt_error = getfield(results, {1}, fnames{1}, {cnvg_idx(1)}, 'rms_pt_error');
if isempty(rms_pt_error)
	error('No rms_pt_error data');
else
	n_samples = length(rms_pt_error);
end

% For each algorithm
hold on;
for l=1:length(alg_idx)
	alg_name = fnames{l};
	disp(alg_name);
	
	mean_rms = zeros(1, n_samples);

	% For each test
	for t=1:length(cnvg_idx)
		rms_pt_error = getfield(results, {1}, alg_name, {cnvg_idx(t)}, 'rms_pt_error');
		mean_rms = mean_rms + rms_pt_error;
	end
	
	mean_rms = mean_rms / length(cnvg_idx);
	
	plot(mean_rms, lstyle{l}, 'linewidth', lwidth(l), 'markersize', msize(l));
	
%	legend_string = strvcat(legend_string, alg_name);
end

%title([warp,': point sigma = 2:2:10, pixel sigma = ',num2str(pix_sigma)]);
xlabel('Iteration');
ylabel('RMS Point Error');
axis([0.75 15.25 -0.25 10]);
legend(legend_string);
drawnow;

hold off;
