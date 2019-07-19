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
    dy = params.dy;
    
    L = sqrt((xR-xL)^2 + (yR-yL)^2);
    theta = abs(atan((yR-yL)/(xR-xL)));
    
    if dx * cos(theta) >= eps && dy * sin(theta) >= eps
        Nxl = ceil(abs((xR-xL)/(dx*cos(theta))));
        dxl = (xR-xL)/Nxl;
        xi = xL:dxl:xR;
        Nyl = ceil(abs((yR-yL)/(dx*sin(theta))));
        dyl = (yR-yL)/Nyl;
        eta = yL:dyl:yR;
    elseif dx * cos(theta) <= eps && dy * sin(theta) >= eps
        Nyl = ceil(abs((yR-yL)/(dx*sin(theta))));
        dyl = (yR-yL)/Nyl;
        eta = yL:dyl:yR;
        xi = xL;
    elseif dx * cos(theta) >= eps && dy * sin(theta) <= eps
        Nxl = ceil(abs((xR-xL)/(dx*cos(theta))));
        dxl = (xR-xL)/Nxl;
        xi = xL:dxl:xR;
        eta = yL;
    else
        xi = xL;
        eta = yL;
    end
end