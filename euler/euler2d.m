classdef euler2d
    %EULER2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        u,v,pressure,rho,size_x,size_y,dx,target_dt,particle,obstacle
    end
    
    methods
        function obj = euler2d(size_x_,size_y_,dx_,target_dt_)
            %EULER2D Construct an instance of this class
            %   Detailed explanation goes here
            obj.size_x = size_x_;
            obj.size_y = size_y_;
            obj.dx = dx_;
            obj.target_dt = target_dt_;
            obj.u = zeros(size_x_,size_y_);
            obj.v = zeros(size_x_,size_y_);
            obj.pressure = zeros(size_x_,size_y_);
            obj.particle = zeros(size_x_,size_y_);
            obj.obstacle = zeros(size_x_,size_y_);
            obj.rho = 1;
        end
        
        function dt = get_safe_dt(obj,target_dt)
            max_uv_norm = 0;
            for i=2:obj.size_x-1
                for j=2:obj.size_y-1
                    max_uv_norm = max(max_uv_norm,sqrt((obj.u(i,j)^2 + obj.v(i,j)^2)));
                end
            end
            u_max = max_uv_norm+sqrt(5*obj.dx*9.8);
            dt = 5*obj.dx/u_max;
            if dt > target_dt
                dt = target_dt;
            end
        end

        function nobj = project(obj,dt)
            for itr = 1:500
                for i=2:obj.size_y-1
                    for j=2:obj.size_x-1
                        obj.pressure(i,j) = -obj.rho/dt*obj.dx*(obj.u(i,j)-obj.u(i,j-1)+obj.v(i,j)-obj.v(i-1,j))...
                            + obj.pressure(i+1,j) + obj.pressure(i-1,j) + obj.pressure(i,j+1) + obj.pressure(i,j-1);
                        obj.pressure(i,j) = obj.pressure(i,j) *0.25;
                    end
                end
            end

            for i=2:obj.size_y-1
                for j=2:obj.size_x-1
                    obj.u(i,j) = obj.u(i,j) - dt/obj.rho*(obj.pressure(i,j+1)-obj.pressure(i,j))/obj.dx;
                    obj.v(i,j) = obj.v(i,j) - dt/obj.rho*(obj.pressure(i+1,j)-obj.pressure(i,j))/obj.dx;
                end
            end
            nobj = obj;
        end

        function nobj = advection(obj,dt)
            tmpu = zeros(obj.size_x,obj.size_y);
            tmpv = zeros(obj.size_x,obj.size_y);
            for i=2:obj.size_y-1
                for j=2:obj.size_x-1
                    ou = obj.u(i,j);
                    ov = obj.v(i,j);
                    px = (j+0.5)*obj.dx - 0.5*ou*dt;
                    py = (i+0.5)*obj.dx - 0.5*ov*dt;
                    pj = floor(px/obj.dx);
                    pj = max(pj,1);
                    pj = min(pj,obj.size_x);
                    pi = floor(py/obj.dx);
                    pi = max(pi,1);
                    pi = min(pi,obj.size_y);
                    ou = obj.u(pi,pj);
                    ov = obj.v(pi,pj);
                    px = px - 0.5*ou*dt;
                    py = py - 0.5*ov*dt;
                    pj = floor(px/obj.dx);
                    pj = max(pj,1);
                    pj = min(pj,obj.size_x);
                    pi = floor(py/obj.dx);
                    pi = max(pi,1);
                    pi = min(pi,obj.size_y);
                    tmpu(i,j) = obj.u(pi,pj);
                    tmpv(i,j) = obj.v(pi,pj);
                end
            end
            obj.u = tmpu;
            obj.v = tmpv;
            nobj = obj;
        end
        function nobj = advection_particle(obj,dt)
            tmp_particle = zeros(obj.size_x,obj.size_y);
            for i=2:obj.size_y-1
                for j=2:obj.size_x-1
                    ou = obj.u(i,j);
                    ov = obj.v(i,j);
                    px = (j+0.5)*obj.dx - 0.5*ou*dt;
                    py = (i+0.5)*obj.dx - 0.5*ov*dt;
                    pj = floor(px/obj.dx);
                    pj = max(pj,1);
                    pj = min(pj,obj.size_x);
                    pi = floor(py/obj.dx);
                    pi = max(pi,1);
                    pi = min(pi,obj.size_y);
                    ou = obj.u(pi,pj);
                    ov = obj.v(pi,pj);
                    px = px - 0.5*ou*dt;
                    py = py - 0.5*ov*dt;
                    pj = floor(px/obj.dx);
                    pj = max(pj,1);
                    pj = min(pj,obj.size_x);
                    pi = floor(py/obj.dx);
                    pi = max(pi,1);
                    pi = min(pi,obj.size_y);
                    tmp_particle(i,j) = obj.particle(pi,pj);
                end
            end
            obj.particle = tmp_particle;
            nobj = obj;
        end
    end
end

