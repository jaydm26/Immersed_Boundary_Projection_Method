function [map,L,theta] = Line_Builder(xL,xR,yL,yR)
    
    % Builds a bunch of points lying on the line dx distance apart from
    % each other. Map gives the X and Y coordinates of the line. L gives
    % the length. Theta gives the angle of the line with the x-axis.
    
    global dx dy
    
    L = sqrt((xR-xL)^2 + (yR-yL)^2);
    theta = abs(atan((yR-yL)/(xR-xL)));
    
    if dx * cos(theta) >= eps && dy * sin(theta) >= eps
        Nxl = ceil(abs((xR-xL)/(dx*cos(theta))));
        dxl = (xR-xL)/Nxl;
        map(:,1) = xL:dxl:xR;
        Nyl = ceil(abs((yR-yL)/(dx*sin(theta))));
        dyl = (yR-yL)/Nyl;
        map(:,2) = yL:dyl:yR;
    elseif dx * cos(theta) <= eps && dy * sin(theta) >= eps
        Nyl = ceil(abs((yR-yL)/(dx*sin(theta))));
        dyl = (yR-yL)/Nyl;
        map(:,2) = yL:dyl:yR;
        map(:,1) = xL;
    elseif dx * cos(theta) >= eps && dy * sin(theta) <= eps
        Nxl = ceil(abs((xR-xL)/(dx*cos(theta))));
        dxl = (xR-xL)/Nxl;
        map(:,1) = xL:dxl:xR;
        map(:,2) = yL;
    else
        map(:,1) = xL;
        map(:,2) = yL;
    end
end