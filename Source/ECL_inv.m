function [x,y] = ECL_inv(gamma)

    global g_hat body_map Nx Ny
    
    sf = NodeData(Nx,Ny);
    gamma_hat = fft2(gamma.x, 2*Nx+1, 2*Ny+1);
    sf_hat = gamma_hat .* g_hat;
    sf_hat = ifft2(sf_hat);
    sf.x = sf_hat(Nx+1:end,Ny+1:end);
    
    q = curl_2(sf);
    
    x = zeros(length(body_map(:,1)),1);
    y = zeros(length(body_map(:,2)),1);
    for k = 1:length(body_map(:,1))
        x(k,1) = E_op(body_map(k,1),body_map(k,2),q,1);
        y(k,1) = E_op(body_map(k,1),body_map(k,2),q,2);
    end
end