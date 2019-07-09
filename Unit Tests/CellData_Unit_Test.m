%% CellData Unit Test

flag = 0;
Nx = 3;
Ny = 3;

C = CellData(Nx,Ny);

if C.data ~= "cell" 
    flag = 1;
end

for i = 1:Nx+2
    for j = 1:Ny+2
        if C.x(i,j) ~= 0
            flag = 1;
        end
    end
end

[a,b] = C.size;

if a ~= Nx || b ~= Ny
    flag = 1;
end

if flag == 1
    disp("CellData is not working properly. Please check the code, or download a copy of the code from the git repository");
else
    disp("CellData is not working properly. Please check the code, or download a copy of the code from the git repository");
end