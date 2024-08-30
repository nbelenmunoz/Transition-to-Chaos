clear all
close all
clc

alpha = 0;   
beta = 0;    

x_range = linspace(-1, 1, 20); 
y_range = linspace(-1, 1, 20); 

[x, y] = meshgrid(x_range, y_range);

% Differential equations 
dx = y + alpha*x;
dy = -x + beta*y;

% Plot the phase portrait
figure;
quiver(x, y, dx, dy); 
hold on;
plot(x_range, zeros(size(x_range)), 'k--'); 
plot(zeros(size(y_range)), y_range, 'k--'); 
xlabel('x');
ylabel('dy/dt');
title('Unperturbed Phase Portrait');
axis tight;
grid on;
