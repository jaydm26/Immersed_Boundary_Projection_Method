global body_map

k = length(body_map(:,1));
A_gen = zeros(2*k,2*k);

for i = 1:2*k
    x = zeros(2*k,1);
    x(i) = 1;
    X = afun(x);
    A_gen(:,i) = X;
end

disp('Done')