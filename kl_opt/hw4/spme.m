function [theta, obj_val] = spme(S, num_nonzero, rho, r, max_iter)
%SPPME Sparse precision matrix estimation
    %   Uses proximal distance algorithm
    
    % initialization
    p = size(S,1);
    theta = eye(p);                       % initialization of solution
    obj_val = zeros(max_iter+1, 1);       % store objective values

    % minimize surrogate for each rho
    for i=1:max_iter

        % compute objective value
        obj_val(i) = -log(det(theta))+trace(S*theta);

        % project theta
        theta = project(theta, num_nonzero);

        % update theta
        [U, D] = eig(S - rho*theta);
        e_vector = zeros(p,1);
        for j=1:p
            e_vector(j) = (-D(j,j)+sqrt(4*rho+D(j,j)^2))/(2*rho);
        end
        theta = U'*diag(e_vector)*U;

        % update rho
        rho = rho * r;
    end
    
    % project theta onto T^p_k
    theta = project(theta, num_nonzero);
    obj_val(max_iter+1) = -log(det(theta))+trace(S*theta);
    
end

