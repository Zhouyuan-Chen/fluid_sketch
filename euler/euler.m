x_len = 1.0;
y_len = 1.0;

dx = 0.01;
target_dt = 0.1;
rho = 1;

% compute the optimal dt
solver = euler2d(x_len/dx,y_len/dx,dx,target_dt);
solver.v(20:50,20:60)=2;
solver.particle(20:50,20:60) = 10;
solver.v(51:81,40:80)=-2;
solver.particle(51:81,40:80) = -10;
dt = solver.get_safe_dt(target_dt);

for itr = 1:50
    imagesc (solver.particle)    
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