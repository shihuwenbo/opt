% initialization
y = 2.0; c = 2.0; x = -5.0; z = 0.0; r = 0.9;
eps = 10^-8; max_iter = 100;

% plot function
xmin = -6; xmax = 6; npts = 1000;
pts = linspace(xmin, xmax, npts);
fx = -c*z*pts + c*log(1+exp(pts)) + 0.5*(-pts+y).^2;
figure; plot(pts, fx, 'b-'); xlim([xmin xmax]); hold on;

% find prox using newton's method
for i=1:max_iter
    
    % compute first and second derivative
    derv1 = -c*z + c*exp(x)/(1+exp(x)) + x - y;
    derv2 = c*exp(x)/(1+exp(x))^2 + 1;

    % backtracking line search
    t = 1; next_x = x - t*derv1/derv2;
    cur_obj = -c*z*x + c*log(1+exp(x)) + 0.5*(y-x)^2;
    while(cur_obj < -c*z*next_x + c*log(1+exp(next_x)) + 0.5*(y-next_x)^2)
        t = t*r; next_x = x - t*derv1/derv2;
    end
    
    % plot step
    if(i <= 3), text(x, cur_obj+1.0, sprintf('%d', i),'fontsize', 25); hold on; end
    
    % check stop condition
    new_x = x - t*derv1/derv2;
    new_obj = -c*z*new_x + c*log(1+exp(new_x)) + 0.5*(y-new_x)^2;
    if(cur_obj - new_obj < eps || abs(new_x-x) < eps), break; end
    
    % update x
    x = new_x;
    
end

% save figure
xlabel('x', 'fontsize', 25); ylabel('f(x)', 'fontsize', 25);
title(sprintf('z=%.2f, c=%.2f, y=%.2f, prox_{cf}(y)=%.4f', z,c,y,x),...
    'fontsize', 25);
set(gca,'FontSize',25)
print(sprintf('prob_8_y_%d_c_%d_z_%d', y, c, z),'-depsc','-r0');