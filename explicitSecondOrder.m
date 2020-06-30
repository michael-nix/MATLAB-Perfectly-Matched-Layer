clear;   close all;   showpotential = true;

%% Initialize computational domain:
n = 256;   c = 1;   tmax = 2^12;

xmin = -50; xmax = 50;   x = linspace(xmin, xmax, n);   dx = x(2) - x(1);
[x, y] = meshgrid(x, x);

phi = exp(-(x.^2+y.^2)/2)/(2*pi);

%% Initialize perfectly matched layer & simulation variables:
[sigmax, sigmay] = setupPML(x, dx);

dt = 0.5 * dx / c;
s_xplusy = dt/2 * (sigmax + sigmay);
s_xtimesy = dt^2 * sigmax .* sigmay;

u_prev = zeros(size(x));
u_now = zeros(size(x));

vx = zeros(size(x));
vy = zeros(size(x));

s = c^2 * dt^2;

%% Run the simulation:
if showpotential
    fig = mesh(x, y, u_now);
    axis([xmin xmax xmin xmax -0.2 0.2]);
    xlabel('x');   ylabel('y');
    zlabel('Wave Amplituide');
    title('Generic Wave with PML at Boundary');
end

tavg = 0;
for t = 1:tmax
    tic;
    
    [dudx, dudy] = gradient(u_now, dx);
    vx = vx + dt*(dudx - sigmax.*vx);
    vy = vy + dt*(dudy - sigmay.*vy);
    
    u_next = (s *(4 * del2(u_now, dx) - divergence(x, y, sigmax.*vx, sigmay.*vy) + cos(dt*t) * phi) + ...
             -(s_xtimesy - 2).*u_now - (1 - s_xplusy).*u_prev)./(1 + s_xplusy);
    
    u_next(1,:) = 0; u_next(:,1) = 0; u_next(end,:) = 0; u_next(:,end) = 0;
    u_prev = u_now;
    u_now = u_next;
    
    tavg = tavg + toc;
    
    if showpotential && ~mod(t, 4)
        fig.ZData = u_now;
        drawnow;
    end
end

disp(['Average single-pass calculation time: ', num2str(tavg / tmax * 1e3), ' ms']);