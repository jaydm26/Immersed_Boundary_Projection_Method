function A = MatrixA_Generator(Type)
    
    global body_map
    
    switch Type
        case "vel"

            k = length(body_map(:,1));
            A = zeros(2*k,2*k);

            for i = 1:2*k
                x = zeros(2*k,1);
                x(i) = 1;
                X = afun(x);
                A(:,i) = X;
            end
        case "temp"
            k = length(body_map(:,1));
            A = zeros(k,k);
            
            for i = 1:k
                x = zeros(k,1);
                x(i) = 1;
                X = afun_temp(x);
                A(:,i) = X;
            end
    end
end