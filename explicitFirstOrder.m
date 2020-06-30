clear;   close all;   showpotential = true;

%% Initialize computational domain:
n = 256;   c = 1;   tmax = 2^12;

xmin = -50; xmax = 50;   x = linspace(xmin, xmax, n);   dx = x(2) - x(1);
[x, y] = meshgrid(x, x);

phi = exp(-(x.^2+y.^2)/2)/(2*pi);

%% Initialize perfectly matched layer & simulation variables:
[sigmax, sigmay] = setupPML(x, dx);

dt = 0.25 * dx / c;
s_xplusy = 1/c^2*(sigmax + sigmay);
s_xtimesy = 1/c^2*(sigmax.*sigmay);

psi = zeros(size(x));
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
    vx = vx + dt * (dudx - vx.*sigmax);
    vy = vy + dt * (dudy - vy.*sigmay);
    
    [dvxdx, ~] = gradient(vx, dx);
    [~, dvydy] = gradient(vy, dx);
    psi = psi + dt * (sigmay.*dvxdx + sigmax.*dvydy - s_xtimesy.*u_now + cos(dt*t) * phi);
    
    u_now = u_now + dt * c^2 * (dvxdx + dvydy - s_xplusy.*u_now + psi);
    
    tavg = tavg + toc;
    
    if showpotential && ~mod(t, 4)
        fig.ZData = u_now;
        drawnow;
    end
end

disp(['Average single-pass calculation time: ', num2str(tavg / tmax * 1e3), ' ms']);