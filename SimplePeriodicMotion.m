clear all;
close all;
clc;

% Parameters
alpha = 0.1;    
r = 1.0;        
beta1 = 0.6;    % Forcing amplitude for Type I motion
omega1 = 4.0;   % Angular frequency for Type I motion
beta2 = 2.3;    % Forcing amplitude for Type II motion
omega2 = 1.95;  % Angular frequency for Type II motion

% Simulation settings
dt = 0.001;      
T1 = 2*pi/omega1;  % Period for Type I motion
T2 = 2*pi/omega2;  % Period for Type II motion
tspan1 = 0:dt:0.95*T1; 
tspan2 = 0.916:dt:1.3*T2; 

% Initialize variables for Type I motion
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
        y1 = -r * y1;
        x1 = sign(x1); 
    end

    data1 = [data1; x1 y1];
end

% Initialize variables for Type II motion
x2 = 0;
y2 = -1;
t2 = 0;
data2 = [x2 y2];

% Simulate Type II motion
for t = tspan2
    y2 = y2 + (dt * ((beta2*cos(omega2 * t)) - (2*alpha*y2) + (x2)));
    x2 = x2 + (dt * y2);
    
    % Check for impact
    if abs(x2) >= 1
        y2 = -r * y2;   
        x2 = sign(x2);  
    end

    data2 = [data2; x2 y2];
end

%%Type I motion
% Plotting
figure()
plot(data1(:,1), data1(:,2), 'k', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('(a) Type I motion, \beta = 0.6, \omega = 4.0, period T');
axis equal;
xlim([-1 1]);
ylim([-1.5 1.5]);

%%Type II motion
% Plotting
figure()
plot(data2(:,1), data2(:,2), 'k', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('(b) Type II motion, \beta = 2.3, \omega = 1.95, period T');
%axis equal;
xlim([-1 1]);
ylim([-4 4]);
