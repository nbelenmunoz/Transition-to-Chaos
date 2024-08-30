
clear all
clc
% Parameters
beta = [0.400, 0.425, 0.450];  % Values of beta close to beta infinity
r = 1.0; 
q = 1.0;  
alpha = 0.1;  
omega = 4.0;  

% Time/phase vector
phi = linspace(0, 1.5, 100);

% Calculate yU(phi) and yS(phi) for different values of beta
for k = 1:length(beta)
    b = beta(k);

    s1 = -alpha + sqrt(alpha^2 + omega^2);
    s2 = -alpha - sqrt(alpha^2 + omega^2);
   
    c1 = 1 - (b/q) * cos(omega * phi); 
    c2 = 1 - (b/q) * cos(omega * phi); 

    yU = s1 * (c1 - (b/q) * omega * sin(omega * phi));
    yS = s2 * (c2 - (b/q) * omega * sin(omega * phi));
    
    %Plotting
    figure;
    hold on;
    plot(phi, yU, 'b', 'LineWidth', 1.5); % Unstable manifold
    plot(phi, yS, 'r', 'LineWidth', 1.5); % Stable manifold
    title(['Trajectories for \beta = ', num2str(b)]);
    xlabel('\Phi');
    ylabel('y');
    legend({'y^U', 'y^S'}, 'Location', 'northeast');
    grid on;
    hold off;
end
