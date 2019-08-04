function video = video_generator(params,domain,filename,xi,eta,ip1,ip2)
    %VIDEO_GENERATOR Generates video of the input field.
    %
    % video = video_generator(params,domain,filename,xi,eta,ip1,ip2)
    %
    % Variable lookup:
    %
    % params: flow parameters.
    %
    % domain: domain parameters.
    %
    % filename: File name.
    %
    % xi, eta: X and Y coordinate of the Lagrangian points.
    %
    % ip1,ip2: input fields.
    %
    % Created by Jay Mehta (18 July 2019)
    
    X_n = domain.X_n;
    Y_n = domain.Y_n;
    x_range = domain.x_range;
    y_range = domain.y_range;
    char_L = params.char_L; 
    U = params.U;
    Nx = domain.Nx;
    Ny = domain.Ny;
    
    xi = [xi;xi(1)];
    eta = [eta;eta(1)];
    
    switch nargin
        case 6
            [t,~,~] = size(ip1);
            s = NodeData(Nx,Ny);
            video = VideoWriter(filename,"MPEG-4");
            video.FrameRate = 1/params.dt;
            open(video);
            for i = 1:t
                F = figure("visible","off");
                s.x(:,:) = ip1(i,:,:);
                contour(X_n./(2*char_L),Y_n./(2*char_L),(s.x'),0:0.5:10)
                hold on
                contour(X_n./(2*char_L),Y_n./(2*char_L),(s.x'),-10:0.5:0,'--')
                plot(xi./(2*char_L),eta./(2*char_L),"LineWidth",2);
                hold off
                pbaspect([1 (y_range(2)-y_range(1))/(x_range(2)-x_range(1)) 1])
                title(strcat('Streamlines for a Cylinder of Diameter = '...
                    ,num2str(2*char_L),' in a flow of uniform velocity of U = ',num2str(U)));
                xlabel("X/D")
                ylabel("Y/D")
                writeVideo(video,getframe(gcf));
                close(F);
            end
            close(video)
            
        case 7
            [t,~,~] = size(ip1);
            qx = NodeData(Nx,Ny);
            qy = NodeData(Nx,Ny);
            velocity = EdgeData(Nx,Ny);
            video = VideoWriter(filename,"MPEG-4");
            video.FrameRate = 1/params.dt;
            open(video);
            for i = 1:t
                velocity.x(:,:) = ip1(i,:,:);
                velocity.y(:,:) = ip2(i,:,:);
                qx = interpol(qx,velocity,1);
                qy = interpol(qy,velocity,2);
                F = figure("visible","off");
                quiver(X_n./(2*char_L),Y_n./(2*char_L),qx.x',qy.x');
                hold on
                plot(xi./(2*char_L),eta./(2*char_L),"LineWidth",2);
                hold off
                pbaspect([1 1 1])
                title(strcat('Velocity around a Cylinder of D = ',num2str(2*char_L),' in a flow of uniform velocity of U = ',num2str(U)));
                xlabel("X/D")
                ylabel("Y/D")
                axis tight
                writeVideo(video,getframe(gcf));
                close(F);
            end
            close(video);
            close all
    end
end