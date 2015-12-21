function theta_proj = project(theta, num_nonzero)
%PROJECT projects theta onto the set T_k^p
%   follows Example 5.5.3

    p = size(theta, 1);
    theta_proj = theta;
    idx = 1;
    abov_diag = zeros((p*p-p)/2, 3);
    for j=1:p
        for k=j+1:p
            abov_diag(idx,1) = abs(theta(j,k));
            abov_diag(idx,2) = j;
            abov_diag(idx,3) = k;
            idx = idx + 1;
        end
    end
    abov_diag = sortrows(abov_diag, -1);
    for j=(num_nonzero+1):length(abov_diag)
        theta_proj(abov_diag(j,2), abov_diag(j,3)) = 0.0;
    end
    
    % for below diagonal elements
    idx = 1;
    below_diag = zeros((p*p-p)/2, 3);
    for j=1:p
        for k=j+1:p
            below_diag(idx,1) = abs(theta(k,j));
            below_diag(idx,2) = k;
            below_diag(idx,3) = j;
            idx = idx + 1;
        end
    end
    below_diag = sortrows(below_diag, -1);
    for j=(num_nonzero+1):length(below_diag)
        theta_proj(below_diag(j,2), below_diag(j,3)) = 0.0;
    end

end