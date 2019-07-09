function out = H_op(x,y,f)
    
    global body_map dx
    
    xi = body_map(:,1);
    eta = body_map(:,2);
    
    out = zeros(length(f),1);
    
    for k = 1:length(f)
        if f(k) == 0 || (x-xi(k))/dx > 1.5 || (y-eta(k))/dx > 1.5
            out(k) = 0;
        else
            out(k) = f(k) * ddf_roma_2D(x-xi(k),y-eta(k));
        end
    end
    out = sum(out);
end