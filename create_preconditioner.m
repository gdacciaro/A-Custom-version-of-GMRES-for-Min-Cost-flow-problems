function [S, P, creation_time_S] = create_preconditioner(D, E, full)
    tic;
        D = diag(D);
        tic;
            S = -E * (D\E');
        creation_time_S = toc; 
        fprintf("S matrix created in %f seconds.\n",creation_time_S);

    if full
        dim = size(D, 1) + size(E, 1);
        P = zeros(dim, dim);
        P(1:size(D, 1), 1:size(D, 1)) = D;
        P(size(D, 1)+1:end, size(E, 2)+1:end) = -S;

    total_precond_time = toc;
        fprintf("P matrix created in %f seconds.\n",total_precond_time);
    else
    total_precond_time = toc;
        P = NaN;
    end
end
