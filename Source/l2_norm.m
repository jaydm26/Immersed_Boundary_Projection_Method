function A = l2_norm(f)
    %L2_NORM Takes the L2 norm of a column vector f.
    %
    % A = l2_norm(f)
    %
    % Created by Jay Mehta (18 July 2019)
    
    A = sqrt(dot_prod(f,f));
end