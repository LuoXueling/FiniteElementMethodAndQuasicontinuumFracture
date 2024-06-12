function [x] = solve_linear_system(K,R,method)
% Solve linear system of equations through direct or distributed method.
if strcmp(method,'direct')
    x=K\R;
elseif strcmp(method,'distributed')
    % ref: https://www.mathworks.com/help/parallel-computing/Use-Distributed-Arrays-to-Solve-Systems-of-Linear-Equations-with-Direct-Methods.html
    K_dist=distributed(K);
    R_dist=distributed(R);
    x=K_dist\R_dist;
    x=gather(x);
elseif strcmp(method,'pcg')
    % ref: https://www.mathworks.com/help/parallel-computing/Use-Distributed-Arrays-to-Solve-Systems-of-Linear-Equations-with-Iterative-Methods.html
    maxiter=floor(sqrt(length(K)));
    [~,x] = evalc('pcg(K,R,[],maxiter)');
elseif strcmp(method,'pcg-dist')
    % ref: https://www.mathworks.com/help/parallel-computing/Use-Distributed-Arrays-to-Solve-Systems-of-Linear-Equations-with-Iterative-Methods.html
    K_dist=distributed(K);
    R_dist=distributed(R);
    maxiter=floor(sqrt(length(K)));
    [~,x] = evalc('pcg(K_dist,R_dist,[],maxiter)');
    x=gather(x);
elseif strcmp(method,'gpu')
    K_gpu=gpuArray(K);
    R_gpu=gpuArray(full(R));
    x=K_gpu\R_gpu;
    x=gather(x);
end

