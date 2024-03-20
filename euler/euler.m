x_len = 1.0;
y_len = 1.0;

dx = 0.01;
target_dt = 0.1;
rho = 1;

% compute the optimal dt
solver = euler2d(x_len/dx,y_len/dx,dx,target_dt);
% solver.v(20:50,20:60)=2;
% solver.particle(20:50,20:60) = 10;
% solver.v(51:81,40:80)=-2;
% solver.particle(51:81,40:80) = -10;


for itr = 1:1000
    solver.u(30:70,2) = 1;
    solver.particle(30:50,2) = -10;
    solver.particle(51:70,2) = 10;

    dt = solver.get_safe_dt(target_dt);

    for i=1:y_len/dx
        for j = 1:x_len/dx
            if((i-50)*(i-50)+(j-30)*(j-30)<100)
                solver.obstacle(i,j) = 1;
                picture_data(i,j) = 50;
                solver.u(i,j) = 0;
                solver.v(i,j) = 0;
            end
            picture_data(i,j) = solver.u(i,j);
        end
    end
    imagesc (picture_data)    
    drawnow
    str = string(itr);
    filename = insertAfter(str,str,".png");
    set(gca,'XTick',[]) % Remove the ticks in the x axis
    set(gca,'YTick',[]) % Remove the ticks in the y axis
    set(gca,'Position',[0 0 1 1]) % Make the axes occupy the hole figure
    saveas(gcf,filename,'png')
    % project
    solver = solver.project(dt);
    % advection
    solver = solver.advection_particle(dt);
    solver = solver.advection(dt);
end

% imagesc (solver.particle)
% imagesc (solver.particle)
% solver = solver.advection_particle(dt);
% imagesc (solver.particle)
% solver = solver.advection_particle(dt);
% imagesc (solver.particle)