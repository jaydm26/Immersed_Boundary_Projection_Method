function op = curl_2(ip)
    %CURL_2 Execute the curl operation depending on the data field of ip. 
    % When ip is an Edge field, it takes a rot curl. When ip is a Node 
    % field, it takes a curl.
    % 
    % op = CURL_2(ip)
    %
    % Variable lookup:
    %
    % ip: Input. Can be an Edge field or a Node field.
    %
    % op: output after taking curl of ip.
    %
    % Created by Jay Mehta (18 July 2019)
    
    switch ip.data 
        case 'edge' 
            if nargin ~= 1
                error("Too many inputs")
            end
            Nx = ip.size(1);
            Ny = ip.size(2);
            op = NodeData(Nx,Ny);
            for i = 1:Nx+1
                for j = 1:Ny+1
                    op.x(i,j) = ip.x(i,j) - ip.x(i,j+1) + ip.y(i+1,j) - ip.y(i,j);
                end
            end
        case 'node'
            if nargin == 1
                Nx = ip.size(1);
                Ny = ip.size(2);
                op = EdgeData(Nx,Ny);
                for i = 1:Nx+1
                    for j = 2:Ny+1
                        op.x(i,j) = ip.x(i,j) - ip.x(i,j-1);
                    end
                end
                for i = 2:Nx+1
                    for j = 1:Ny+1
                        op.y(i,j) = ip.x(i-1,j) - ip.x(i,j);
                    end
                end
%             elseif nargin == 2
%                 Nx = ip.size(1);
%                 Ny = ip.size(2);
%                 op = EdgeData(Nx,Ny);
%                 if dim == 1
%                     for i = 1:Nx+1
%                         for j = 2:Ny+1
%                             op.x(i,j) = ip.x(i,j) - ip.x(i,j-1);
%                         end
%                     end
%                 elseif dim == 2
%                     for i = 2:Nx+1
%                         for j = 1:Ny+1
%                             op.y(i,j) = ip.x(i,j) - ip.x(i-1,j);
%                         end
%                     end
%                 else
%                     error("dim cannot exceed 2.");
%                 end
            else
                error("Too many inputs.")
            end
        otherwise
            error("ip is not an Edge field or a Node field")
    end    
end