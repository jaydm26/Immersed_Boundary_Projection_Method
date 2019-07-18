function A = l2_norm(f)
    %L2_NORM Takes the L2 norm of a column vector f.
    %
    % A = l2_norm(f)
    
    A = sqrt(dot_prod(f,f));
end