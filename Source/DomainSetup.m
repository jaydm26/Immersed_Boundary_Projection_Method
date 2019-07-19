function [X,Y] = DomainSetup(params,domain,type)
    %DOMAINSETUP Creates the domain for the data fields.
    %
    % [X,Y] = DomainSetup(params,domain,type)
    %
    % Variable lookup:
    %
    % params: flow parameters.
    % 
    % domain: domain parameters.
    %
    % type: type of domain- "c" or "cell" for cell
    %                       "xe" or "x-edge" for X Edge Field
    %                       "ye" or "y-edge" for Y Edge Field
    %                       "n" or "node" for node
    %
    % X: Full Domain for X-direction
    %
    % Y: Full Domain for Y-direction
    %
    % Created by Jay Mehta (18 July 2019)
    
    x0 = domain.x_range(1);
    x_max = domain.x_range(2);
    y0 = domain.y_range(1);
    y_max = domain.y_range(2);    
    
    switch lower(type)
        case {"node","nodal","n"}
            x_node_range = x0:params.dx:x_max;
            y_node_range = y0:params.dx:y_max;
            [X,Y] = meshgrid(x_node_range,y_node_range);
        case {"y-edge","ye"}
            x_edge_range_y = x0-0.5*params.dx:params.dx:x_max+0.5*params.dx;
            y_edge_range_y = y0:params.dx:y_max;
            [X,Y] = meshgrid(x_edge_range_y,y_edge_range_y);
        case {"x-edge","xe"}
            x_edge_range_x = x0:params.dx:x_max;
            y_edge_range_x = y0-0.5*params.dx:params.dx:y_max+0.5*params.dx;
            [X,Y] = meshgrid(x_edge_range_x,y_edge_range_x);
        case {"cell","c"}
            x_cell_range = x0-0.5*params.dx:params.dx:x_max+0.5*params.dx;
            y_cell_range = y0-0.5*params.dx:params.dx:y_max+0.5*params.dx;
            [X,Y] = meshgrid(x_cell_range,y_cell_range);
    end
end