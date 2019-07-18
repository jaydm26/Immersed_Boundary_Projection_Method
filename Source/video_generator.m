function video = video_generator(filename,X,Y)
    
    global dt body_map X_n Y_n R U Nx Ny
    
    xi = body_map(:,1);
    xi = [xi;xi(1)];
    eta = body_map(:,2);
    eta = [eta;eta(1)];
    switch nargin
        case 2
            [t,~,~] = size(X);
            s = NodeData(Nx,Ny);
            video = VideoWriter(filename,"MPEG-4");
            video.FrameRate = 1/dt;
            open(video);
            for i = 1:t
                F = figure("visible","off");
                s.x(:,:) = X(i,:,:);
                contour(X_n./(2*R),Y_n./(2*R),s.x');
                hold on
                plot(xi./(2*R),eta./(2*R),"LineWidth",2);
                hold off
                pbaspect([1 1 1]);
                title(strcat('Streamlines for a Flat Plate of D = ',num2str(2*R),' flow is of uniform velocity of U = ',num2str(U)));
                xlabel("X/D");
                ylabel("Y/D");
                writeVideo(video,getframe(gcf));
                close(F);
            end
            close(video)
            
        case 3
            [t,~,~] = size(X);
            qx = NodeData(Nx,Ny);
            qy = NodeData(Nx,Ny);
            velocity = EdgeData(Nx,Ny);
            video = VideoWriter(filename,"MPEG-4");
            video.FrameRate = 1/dt;
            open(video);
            for i = 1:t
                velocity.x(:,:) = X(i,:,:);
                velocity.y(:,:) = Y(i,:,:);
                qx = interpol(qx,velocity,1);
                qy = interpol(qy,velocity,2);
                F = figure("visible","off");
                quiver(X_n./(2*R),Y_n./(2*R),qx.x',qy.x');
                hold on
                plot(xi./(2*R),eta./(2*R),"LineWidth",2);
                hold off
                pbaspect([1 1 1])
                title(strcat('Velocity around a Cylinder of D = ',num2str(2*R),' in a flow of uniform velocity of U = ',num2str(U)));
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