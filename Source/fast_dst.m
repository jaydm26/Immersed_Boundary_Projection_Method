function y = fast_dst(x,dim)
    
    % Perform a DST for x in the required dimension
    
    [Nx,Ny] = size(x);
    y = zeros(Nx,Ny);
    
    switch dim
        case 1
            for i = 1:Nx
                tmp = -imag(fft([0;transpose(x(i,:));transpose(zeros(1,Ny+1))]));
                y(i,:) = tmp(2:Ny+1);
            end
        case 2
            for j = 1:Ny
                tmp = -imag(fft([0;x(:,j);zeros(Nx+1,1)]));
                y(:,j) = tmp(2:Nx+1);
            end
    end
    