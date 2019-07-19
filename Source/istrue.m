function flag = istrue(logic)
    %ISTRUE Checks if logic is true. If logic is true, flag = 1. If logic 
    % is false, flag = 0.
    %
    % flag = istrue(logic)
    %
    % Created by Jay Mehta (18 July 2019)
    
    if logic
        flag = 1;
    else
        flag = 0;
    end
end