function [beta, s1, s2, obj_val] = lmm(y, A, V1, V2)
    % initialization
    eps = 10^-4; max_iter = 1000;
    
    s1 = 1.0; s2 = 1.0;
    omega = s1*V1+s2*V2; omega_inv = pinv(omega);
    beta = pinv(A'*omega_inv*A)*A'*omega_inv*y;
    
    obj_val = zeros(max_iter, 1);
    
    % iterate until convergence
    for i=1:max_iter
        
        % compute objval
        obj_val(i) = -0.5*log(det(omega))...
            -0.5*(y-A*beta)'*omega_inv*(y-A*beta);
        
        % compute new sigma
        u = omega_inv*(y-A*beta);
        s1_next = s1*sqrt(u'*V1*u/trace(omega_inv*V1));
        s2_next = s2*sqrt(u'*V2*u/trace(omega_inv*V2));
        
        % compute new omega
        omega = s1_next*V1+s2_next*V2;
        omega_inv = pinv(omega);
        
        % compute new beta
        beta_next = pinv(A'*omega_inv*A)*A'*omega_inv*y;
    
        % check convergence
        dist = norm([s1 s2 beta']-[s1_next s2_next beta']);
        if(dist < eps), break; end
        
        % update parameters
        s1 = s1_next; s2 = s2_next; beta = beta_next;
    end

    obj_val = obj_val(1:i);
end

