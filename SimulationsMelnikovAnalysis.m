clear all
clc
% Parameters
beta = [0.400, 0.425, 0.450];  % Values of beta close to beta infinity
r = 1.0;   
alpha = 0.1;  
omega = 4.0;  
q = sqrt(((1 + (omega^2))^2) + ((2*alpha*omega)^2)); %1.0;
sol_tan = (2*alpha*omega) / (1 + (omega^2));
psi = atan(sol_tan);
cap_omega = sqrt((alpha^2) + 1);


% Time/phase vector
phi = linspace(0, 1.5, 10000);

% Calculate yU(phi) and yS(phi) for different values of beta
for k = 1:length(beta)
    b = beta(k);

    s1 = -alpha + cap_omega;
    s2 = -alpha - cap_omega;
   
    c1 = 1 - ((b/q) * cos((omega * phi) + psi)); 
    c2 = 1 - ((b/q) * cos((omega * phi) + psi)); 

    yU = (s1 * c1) - ((b/q) * omega * sin((omega * phi) + psi));
    yS = (s2 * c2) - ((b/q) * omega * sin((omega * phi) + psi));
    
    %Plotting
    figure;
    hold on;
    plot(phi, yU, 'b', 'LineWidth', 1.5); % Unstable manifold
    plot(phi, abs(yS), 'r', 'LineWidth', 1.5); % Stable manifold
    title(['Trajectories for \beta = ', num2str(b)]);
    xlabel('\Phi');
    ylabel('y');
    legend({'y^U', 'y^S'}, 'Location', 'northeast');
    ylim([0 1.5])
    grid on;
    hold off;
end
