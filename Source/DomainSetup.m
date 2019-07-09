function [X,Y] = DomainSetup(x_range,y_range,Nx,Ny,type)
    
    global dx dy
    x0 = x_range(1);
    x_max = x_range(2);
    y0 = y_range(1);
    y_max = y_range(2);
    dx = (x_max-x0)/(Nx);
    dy = (y_max-y0)/(Ny);    
    
    switch lower(type)
        case {"node","nodal","n"}
            x_node_range = x0:dx:x_max;
            y_node_range = y0:dy:y_max;
            [X,Y] = meshgrid(x_node_range,y_node_range);
        case {"y-edge","ye"}
            x_edge_range_y = x0-0.5*dx:dx:x_max+0.5*dx;
            y_edge_range_y = y0:dy:y_max;
            [X,Y] = meshgrid(x_edge_range_y,y_edge_range_y);
        case {"x-edge","xe"}
            x_edge_range_x = x0:dx:x_max;
            y_edge_range_x = y0-0.5*dy:dy:y_max+0.5*dy;
            [X,Y] = meshgrid(x_edge_range_x,y_edge_range_x);
        case {"cell","c"}
            x_cell_range = x0-0.5*dx:dx:x_max+0.5*dx;
            y_cell_range = y0-0.5*dy:dy:y_max+0.5*dy;
            [X,Y] = meshgrid(x_cell_range,y_cell_range);
    end
end