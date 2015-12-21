% set up problem
p = 100;                     % dimension of matrix
n = 1000;                    % number of samples
samples = 1.2*rand(n, p);    % simulated random samples
S = cov(samples);            % covariance matrix

% ietrate through different number of k
all_num_nonzero = [10; 100; 1000];
for ii=1:length(all_num_nonzero)
    
    % set up parameters
    rho = 1.0;                             % factor appended for penalty
    r = 1.01;                              % rate of increase for rho
    num_nonzero = all_num_nonzero(ii);     % number of non-zero entries
    max_iter = 300;                        % maximum number of iterations
    
    % estimate sparse precision matrix
    [theta, obj_val] = spme(S, num_nonzero, rho, r, max_iter);

    % create objective value plots
    figure('visible', 'off');
    plot(obj_val, 'b-');
    xlabel('iteration', 'fontsize', 20);
    ylabel('objective value', 'fontsize', 20);
    set(gca,'FontSize',20)
    print(sprintf('prob_28_obj_k_%d', num_nonzero),'-depsc','-r0');
    
    % create sparsity pattern plot
    figure('visible', 'off');
    spy(theta);
    set(gca,'FontSize',20)
    print(sprintf('prob_28_sp_k_%d', num_nonzero),'-depsc','-r0');

end