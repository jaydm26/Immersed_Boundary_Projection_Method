function out = H_op_simple(x,y,f)
    global body_map
    xi = body_map(:,1);
    eta = body_map(:,2);
    out = zeros(length(f),1);
    for k = 1:length(f)
        out(k) = f(k) * ddf_roma_2D(x-xi(k),y-eta(k));
    end
    out = sum(out);
end