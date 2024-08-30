%%Phase diagram of simulations
Simulations, r = 1.0, α = 0.1 
clear all
clc
%Common parameters
alpha = 0.1;    % Damping coefficient
r = 1.0;        % Reflection coefficient
dt = 0.001;      % Time step

%%(a) Period 2T motion, β = 0.65, ω = 4.0. 
beta1 = 0.65;    % Forcing amplitude for Type I motion
omega1 = 4.0;    % Angular frequency for Type I motion
T1 = 2*pi/omega1;  % Period for Type I motion
tspan1 = 0.1:dt:2*T1; 
x1 = -1;
y1 = -0.63;
t1 = 0;
data1 = [x1 y1];
% Simulate Type I motion
for t = tspan1
    y1 = y1 + (dt * ((beta1*cos(omega1*t)) - (2*alpha*y1) + (x1)));
    x1 = x1 + (dt * y1);
    
    % Check for impact
    if abs(x1) >= 1
        y1 = -r * y1;   % Reflect the velocity
        x1 = sign(x1);  % Correct position to exactly 1 or -1
    end
    
    % Store data
    data1 = [data1; x1 y1];
end

% Plotting
figure()
plot(data1(:,1), data1(:,2), 'k', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('(a) Type I motion, \beta = 0.65, \omega = 4.0, period T');
axis equal;
xlim([-1 1]);
ylim([-1.5 1.5]);

%%(b) Unsymmetric period T motion, β = 2.1, ω = 1.8. 
beta2 = 2.3;    % Forcing amplitude for Type II motion
omega2 = 1.8;  % Angular frequency for Type II motion
T2 = 2*pi/omega2;  % Period for Type II motion
tspan2 = 0.988:dt:1.35*T2; % Simulate over 1 period for Type II
x2 = 0;
y2 = -1.3;
t2 = 0;
data2 = [x2 y2];

% Simulate Type II motion
for t = tspan2
    y2 = y2 + (dt * ((beta2*cos(omega2 * t)) - (2*alpha*y2) + (x2)));
    x2 = x2 + (dt * y2);
    
    % Check for impact
    if abs(x2) >= 1
        y2 = -r * y2;   % Reflect the velocity
        x2 = sign(x2);  % Correct position to exactly 1 or -1
    end
    
    % Store data
    data2 = [data2; x2 y2];
end

%Plotting
figure()
plot(data2(:,1), data2(:,2), 'k', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('(b) Type II motion, \beta = 2.3, \omega = 1.8, period T');
xlim([-1 1]);
ylim([-4 4]);

%%c) Apparently chaotic motion, β = 2.0, ω = 4.0. 
beta3 = 2.0;    % Forcing amplitude for Type II motion
omega3 = 4.0;  % Angular frequency for Type II motion

% Simulation settings
dt = 0.001;      % Time step
T3 = 2*pi/omega3;  % Period for Type II motion
tspan3 = 0:dt:50*T3; % Simulate over 1 period for Type II

% Initialize variables for chaotic motion
x3 = 0;
y3 = -1.3;
t3 = 0;
data3 = [x3 y3];

% Simulate Type II motion
for t3 = tspan3
    y3 = y3 + (dt * ((beta3*cos(omega3 * t3)) - (2*alpha*y3) + (x3)));
    x3 = x3 + (dt * y3);
    
    % Check for impact
    if abs(x3) >= 1
        y3 = -r * y3;   % Reflect the velocity
        x3 = sign(x3);  % Correct position to exactly 1 or -1
    end
    
    % Store data
    data3 = [data3; x3 y3];
end

% Plotting
figure()
plot(data3(:,1), data3(:,2), 'k', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('(c) Chaotic motion, \beta = 2.0, \omega = 4.0, period T');
xlim([-1 1]);
ylim([-4 4]);
