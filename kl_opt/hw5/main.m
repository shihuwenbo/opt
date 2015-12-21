% simulate data
n = 100;                    % number of simulations
m = 100;                    % number of samples
p = 10;                     % number of dimensions
sigma_1 = 2.0;              % simulated sigma_1
sigma_2 = 2.0;              % simulated sigma_2
sim_rec = zeros(n, 5);      % record simulation result

% simulate data and estimate sigma_1 and sigma_2
for i=1:n
    
    % simulate y
    beta = rand(p, 1);       % beta vector for simulation
    A  = rand(m, p);         % measurement matrix
    V1 = eye(m);             % variance component 1
    V2 = diag(1:m)/100;      % variance component 2 (time)
    y  =  mvnrnd(A*beta, sigma_1*V1+sigma_2*V2)';
    
    % esimate beta, sigma_1, and sigma_2
    [beta_est, s1_est, s2_est, obj_val] = lmm(y, A, V1, V2);
    if(i == 1)
        figure('visible', 'off'); plot(obj_val, 'b-');
        xlim([1 size(obj_val,1)]); xlabel('iteration','fontsize',16);
        ylabel('objective value','fontsize',16); set(gca,'FontSize',12)
        print('prob_16_obj','-depsc','-r0');
    end
    
    % record result
    sim_rec(i,:) = [sigma_1 s1_est sigma_2 s2_est norm(beta-beta_est)];
end

% print summary
fprintf('s1: %.3f\n', sigma_1);
fprintf('estimated s1: %.3f %.3f\n', mean(sim_rec(:,2)), ...
    sqrt(var(sim_rec(:,2))/n));

fprintf('\n');

fprintf('s2: %.2f\n', sigma_2);
fprintf('estimated s2: %.3f %.3f\n', mean(sim_rec(:,4)), ...
    sqrt(var(sim_rec(:,4))/n));

fprintf('\n');

fprintf('s1/(s1+s2): %.3f\n', sigma_1/(sigma_1+sigma_2));
fprintf('estimated s1/(s1+s2): %.3f %.3f\n', ...
    mean(sim_rec(:,2)./(sim_rec(:,2)+sim_rec(:,4))),...
    sqrt(var(sim_rec(:,2)./(sim_rec(:,2)+sim_rec(:,4)))/n));