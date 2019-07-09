function A = l2_norm(f)
    % Takes the L2 norm of a column vector
        A = sqrt(dot_prod(f,f));
end