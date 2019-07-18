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
    
    X_n = domain.X_n;
    Y_n = domain.X_n;
    char_L = params.char_L; 
    U = params.U;
    Nx = ip1.size(1);
    Ny = ip1.size(2);
    
    xi = [xi;xi(1)];
    eta = [eta;eta(1)];
    
    switch nargin
        case 2
            [t,~,~] = size(ip1);
            s = NodeData(Nx,Ny);
            video = VideoWriter(filename,"MPEG-4");
            video.FrameRate = 1/params.dt;
            open(video);
            for i = 1:t
                F = figure("visible","off");
                s.x(:,:) = ip1(i,:,:);
                contour(X_n./(2*char_L),Y_n./(2*char_L),s.x');
                hold on
                plot(xi./(2*char_L),eta./(2*char_L),"LineWidth",2);
                hold off
                pbaspect([1 1 1]);
                title(strcat('Streamlines for a Flat Plate of D = ',num2str(2*char_L),' flow is of uniform velocity of U = ',num2str(U)));
                xlabel("X/D");
                ylabel("Y/D");
                writeVideo(video,getframe(gcf));
                close(F);
            end
            close(video)
            
        case 3
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