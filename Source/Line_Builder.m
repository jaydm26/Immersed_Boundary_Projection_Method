function [xi,eta,L,theta] = Line_Builder(params,xL,xR,yL,yR)
    %LINE_BUILDER Builds a set of points lying on the line dx distance apart from
    % each other.
    %
    % [xi,eta,L,theta] = Line_Builder(params,xL,xR,yL,yR)
    %
    % Variable lookup:
    %
    % xi: X-coordinate of the Lagrangian points.
    %
    % eta: Y-corrdinate of the Lagrangian points.
    %
    % L: length of the line.
    %
    % theta: angle of the line with respect to the X-axis.
    %
    % params: flow parameters
    %
    % xL,yL: Starting coordinates of the line.
    %
    % xR,yR: Ending coordinates of the line.
    %
    % Created by Jay Mehta (18 July 2019)
    
    dx = params.dx;
    dy = params.dx;
    
    L = sqrt((xR-xL)^2 + (yR-yL)^2);
    theta = abs(atan((yR-yL)/(xR-xL)));
    
    if dx * cos(theta) >= eps && dy * sin(theta) >= eps
        if mod(abs((xR-xL)),abs(dx*cos(theta)))/(abs(dx*cos(theta))) > 0.475
            Nxl = ceil(abs((xR-xL)/(dx*cos(theta))));
        else
            Nxl = floor(abs((xR-xL)/(dx*cos(theta))));
        end
        dxl = (xR-xL)/Nxl;
        xi = xL:dxl:xR;
        xi = xi';
        
        if mod(abs((yR-yL)),abs(dx*sin(theta)))/(abs(dx*sin(theta))) > 0.475 * dx
            Nyl = ceil(abs((yR-yL)/(dx*sin(theta))));
        else
            Nyl = floor(abs((yR-yL)/(dx*sin(theta))));
        end
        dyl = (yR-yL)/Nyl;
        eta = yL:dyl:yR;
        eta = eta';
        
    elseif dx * cos(theta) <= eps && dy * sin(theta) >= eps
        if mod(abs((yR-yL)),abs(dx*sin(theta)))/(abs(dx*sin(theta))) > 0.475 * dx
            Nyl = ceil(abs((yR-yL)/(dx*sin(theta))));
        else
            Nyl = floor(abs((yR-yL)/(dx*sin(theta))));
        end
        dyl = (yR-yL)/Nyl;
        eta = yL:dyl:yR;
        eta = eta';
        xi = xL * ones(length(eta),1);
    elseif dx * cos(theta) >= eps && dy * sin(theta) <= eps
        if mod(abs((xR-xL)),abs(dx*cos(theta)))/(abs(dx*cos(theta))) > 0.475
            Nxl = ceil(abs((xR-xL)/(dx*cos(theta))));
        else
            Nxl = floor(abs((xR-xL)/(dx*cos(theta))));
        end
        dxl = (xR-xL)/Nxl;
        xi = xL:dxl:xR;
        xi = xi';
        eta = yL * ones(length(xi),1);
    else
        xi = xL;
        eta = yL;
    end
end