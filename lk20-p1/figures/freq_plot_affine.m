function freq_plot_affine(dname, i_sigma, t_sigma)
% FREQ_PLOT_AFFINE(DATA_NAME, IMAGE_SIGMA, TEMPLATE_SIGMA)
% Plot freqence of convergence test results
%
% e.g. to plot results on takeo data for no image noise and no template noise:
%      freq_plot_affine('takeo', 0, 0);

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: freq_plot_affine.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% Load LK20 line styles
lk20style;
set(gca,'FontSize',font_size);

% Make assumptions about alogorithms run! Easier than trying to compute...
lstyle = {ls.fa, ls.fc, ls.ia, ls.ic};
lwidth = [lw.fa, lw.fc, lw.ia, lw.ic];
msize = [ms.fa, ms.fc, ms.ia, ms.ic];
legend_string = {'Forwards Additive', 'Forwards Compositional' ,'Inverse Additive', 'Inverse Compositional'};

% Load first example
string = ['load ',dname,'_I',num2str(i_sigma),'_T',num2str(t_sigma),'_S',num2str(1,'%.1f'),'.mat'];	
disp(string);
eval(string);

% Get fields
fnames = fieldnames(results);
alg_idx = 1:length(fnames);
idx = strmatch('all_converged_idx', fnames);
alg_idx(idx) = [];

% This needs to match experiment parameters in run_affine
sigma_idx = [1:0.5:10];

% For each spatial sigma
for s=1:length(sigma_idx)
	s_sig = sigma_idx(s);

	% Load data
	disp(['Sigma: ',num2str(s_sig)])
	string = ['load ',dname,'_I',num2str(i_sigma),'_T',num2str(t_sigma),'_S',num2str(s_sig,'%.1f'),'.mat'];	
	disp(string);
	eval(string);

	% For each algorithm
	for l=1:length(alg_idx)
		alg_name = fnames{l};
		disp(alg_name);
		
		n_div = getfield(results, {1}, alg_name, {1}, 'n_diverge');
		n_cnv = getfield(results, {1}, alg_name, {1}, 'n_converge');

		p_cnvg(l, s) = (n_cnv / (n_cnv + n_div)) * 100;
	end			
end

% Plot results
hold on;
for l=1:length(alg_idx)
	plot(p_cnvg(l, :), lstyle{l}, 'linewidth', lwidth(l), 'markersize', msize(l));
end

xlabel('Point Sigma');
ylabel('% Converged');
legend(legend_string, 3);
s = 1:10;
set(gca,'XTick', [1:2:20]);
set(gca,'XTickLabel',s);
axis([0.5 19.5 -2 102]);
drawnow;

hold off;
