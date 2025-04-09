function [x_ref, y_ref, th_ref, v_ref, w_ref] = Path(dt, traj)

    if traj == "circle"
        theta   = 0:dt:2*pi;
        x_ref   = sin(theta);
        y_ref   = cos(theta);
    elseif traj == "u"
        a = 0;
        x_ref = [];
        y_ref = [];
        theta = [];

        while(a < 6)
            x_ref = [x_ref; a*cos(pi/2)];
            y_ref = [y_ref; a*sin(pi/2)];
            theta = [theta; pi/2];
            a     = a+dt;
        end

        angle = pi - dt;

        while(angle > (pi/2 + dt))
            x_ref = [x_ref; 0.5 + 0.5 * cos(angle)];
            y_ref = [y_ref; 6 + 0.5 * sin(angle)];
            theta = [theta; angle];
            angle = angle - (dt/0.5);
        end

        a = 0;

        while(a < 8)
            x_ref = [x_ref; 0.5 + a*cos(0)];
            y_ref = [y_ref; 6.5 + a*sin(0)];
            theta = [theta; 0];
            a     = a + dt;
        end

        angle = pi/2 - dt;

        while(angle > dt)
            x_ref = [x_ref; 8.5 + 0.5 * cos(angle)];
            y_ref = [y_ref; 6 + 0.5 * sin(angle)];
            theta = [theta; angle];
            angle = angle - (dt/0.5);
        end

        a = 0;

        while(a < 6)
            x_ref = [x_ref; 9 + a*cos(-pi/2)];
            y_ref = [y_ref; (a-6)*sin(-pi/2)];
            theta = [theta; -pi/2];
            a     = a + dt;
        end

        x_ref = x_ref';
        y_ref = y_ref';
        theta = theta';

    elseif traj == "infinito"
        eta     = 1;
        alpha   = 10;
        theta   = 0:dt:2*pi*alpha*2;

        x_ref   = eta*sin(theta/alpha);
        y_ref   = eta*sin(theta/(2*alpha));
    end

    n = length(x_ref);

    dx_ref = zeros(1,n);
    dy_ref = zeros(1,n);
    ddx_ref = zeros(1,n);
    ddy_ref = zeros(1,n);

    for i=2:n-1
        dx_ref(i) = (x_ref(i+1) - x_ref(i-1)) / (2*dt);
        dy_ref(i) = (y_ref(i+1) - y_ref(i-1)) / (2*dt);
    end

    % Difference for the first and last points
    % Taylor expansion first order
    dx_ref(1) = (x_ref(2) - x_ref(1)) / dt;
    dy_ref(1) = (y_ref(2) - y_ref(1)) / dt;
    dx_ref(n) = (x_ref(n) - x_ref(n-1)) / dt;
    dy_ref(n) = (y_ref(n) - y_ref(n-1)) / dt;

    th_ref = zeros(1, n);

    th_ref(1) = atan2(dy_ref(1), dx_ref(1));
    for i=2:n
        th_ref(i) = atan2(dy_ref(i), dx_ref(i));
        if (th_ref(i) - th_ref(i-1) > pi)
            while(th_ref(i) - th_ref(i-1) > pi)
                th_ref(i) = th_ref(i) - 2*pi;
            end
        elseif (th_ref(i) - th_ref(i-1) < -pi)
            while(th_ref(i) - th_ref(i-1) < -pi)
                th_ref(i) = th_ref(i) + 2*pi;
            end
        end
    end

    for i=2:n-1
        ddx_ref(i) = (dx_ref(i+1) - dx_ref(i-1)) / (2*dt);
        ddy_ref(i) = (dy_ref(i+1) - dy_ref(i-1)) / (2*dt);
    end

    % Difference for the first and last points
    % Taylor expansion first order
    ddx_ref(1) = (dx_ref(2) - dx_ref(1)) / dt;
    ddy_ref(1) = (dy_ref(2) - dy_ref(1)) / dt;
    ddx_ref(n) = (dx_ref(n) - dx_ref(n-1)) / dt;
    ddy_ref(n) = (dy_ref(n) - dy_ref(n-1)) / dt;

    v_ref = zeros(1,n);
    w_ref = zeros(1,n);
 
    % disp(dy_ref)

    % Velocities
    for i=1:n
        v_ref(i) = sqrt(dx_ref(i)*dx_ref(i) + dy_ref(i)*dy_ref(i));
        % To avoid division by zero
        if v_ref(i) > 0.001
            w_ref(i) = (ddy_ref(i)*dx_ref(i) - ddx_ref(i)*dy_ref(i)) / (dx_ref(i)*dx_ref(i) + dy_ref(i)*dy_ref(i));
        else
            w_ref(i) = 0;
        end
    end

    % v_ref = 0.3*v_ref;
    % w_ref = 0.3*w_ref;

end