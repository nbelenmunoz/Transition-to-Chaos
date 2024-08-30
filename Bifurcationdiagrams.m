%%Bifurcation diagrams in (β, ω) space with r=1.0 and α =0.2
% Clear previous data and close figures
clear all
clc

%Parameters
n = 1:5;
alfa = 0.2; 
r = 1; 

%Frequency range
dw = 0.0001;
w = 0.1:dw:5.0;

Hn = zeros(length(n), length(w));
Gn = zeros(length(n), length(w));
In = zeros(length(n), length(w));

Jn = zeros(length(n), length(w));
Kn = zeros(length(n), length(w));
Ln = zeros(length(n), length(w));

q = zeros(1, length(w)); 

% Loop over each subharmonic order
for i = 1:length(n)
    for j = 1:length(w)
        T = (2 * pi) / w(j);
        cap_w = sqrt(1 + (alfa^2));
        
        % Calculate q for each frequency
        q(j) = sqrt(((1 + (w(j)^2))^2) + ((2 * alfa * w(j))^2));
        s_om = sinh(n(i) * cap_w * T);
        s_om_2 = sinh(n(i) * cap_w * T/2);
        c_om = cosh(n(i) * cap_w * T);
        c_om_2 = cosh(n(i) * cap_w * T/2);
        s_al = sinh(n(i) * alfa * T);
        s_al_2 = sinh(n(i) * alfa * T/2);
        c_al = cosh(n(i) * alfa * T);
        c_al_2 = cosh(n(i) * alfa * T/2);
        
        % Compute Hn, Gn, and In
        Hn(i, j) = c_om - c_al; % H_n
        Gn(i, j) = ((1+r)/(2*cap_w)) * s_om; % G_n
        In(i, j) = (((1-r)/2) * Hn(i, j)) + (((1+r)/2) * (((alfa/cap_w)*(s_om)) - (s_al)));

        % Compute Jn, Ln, and Kn
        Jn(i, j) = c_om_2 + c_al_2;
        Ln(i, j) = ((1+r)/(2*cap_w)) * s_om_2;
        Kn(i, j) = (((1-r)/2) * Jn(i, j)) + (((1+r)/2) * (((alfa/cap_w)*(s_om_2)) + (s_al_2)));

    end
end

% Preallocate matrices for additional results
yn_SN = zeros(length(w), length(n));
yn_II_SN = zeros(length(w), length(n));
yn_PD = zeros(length(w), length(n));
yn_PF = zeros(length(w), length(n));
betaI_SN = zeros(length(w), length(n));
betaII_SN = zeros(length(w), length(n));
betaI_PD = zeros(length(w), length(n));
betaII_PF = zeros(length(w), length(n));
betaI_inf = zeros(length(w), length(n));
betaII_inf = zeros(length(w), length(n));

% Compute yn_SN, yn_PD, betaI_SN, betaI_PD, and betaI_inf
for ii = 1:length(n)
    for jj = 1:length(w)
        % Solving the quadratic equation for y_n (Saddle-node and Period Doubling)
        A = (In(ii, jj)^2) + (w(jj)^2 * Gn(ii, jj)^2);
        B = -2 * Gn(ii, jj) * Hn(ii, jj) * w(jj)^2;
        
        % Solve for betaI_SN first to use in C
        betaI_SN(jj, ii) = (q(jj) * abs(In(ii, jj))) / sqrt(A);
        
        % Calculate C using betaI_SN
        C = w(jj)^2 * Hn(ii, jj)^2 * (1 - (betaI_SN(jj, ii)^2 / q(jj)^2));
        
        % Solve the quadratic equation for yn_SN and yn_PD
        yn_SN(jj, ii) = (+s_om*(1+r)*(w(jj)^2)*Hn(ii, jj)) / ((cap_w*Hn(ii, jj)*(c_al+s_al+((r^2)*(c_al-s_al)))) + ((s_om*(1+r)*((Gn(ii, jj)*(1+(w(jj)^2)))+(2*alfa*In(ii, jj)))) - (2*r*cap_w*Hn(ii, jj)*c_om)));
        yn_PD(jj, ii) = (-s_om*(1+r)*(w(jj)^2)*Hn(ii, jj)) / ((cap_w*Hn(ii, jj)*(c_al+s_al+((r^2)*(c_al-s_al)))) - ((s_om*(1+r)*((Gn(ii, jj)*(1+(w(jj)^2)))+(2*alfa*In(ii, jj)))) - (2*r*cap_w*Hn(ii, jj)*c_om)));

        % Calculate betaI_PD
        betaI_PD(jj, ii) = q(jj) * ((1.0 + (((A * yn_PD(jj, ii)^2) + (B * yn_PD(jj, ii))) / (w(jj)^2 * Hn(ii,jj)^2)))^0.5); 
        
        % Calculate betaI_inf
        betaI_inf(jj, ii) = (q(jj) * abs(((1 - r) * cap_w) + ((1 + r) * alfa))) / sqrt((((1 - r) * cap_w) + ((1 + r) * alfa))^2 + (w(jj) * (1 + r))^2);
        
        % Solution of the Type II Motion
        Delta = -(((2*alfa*Kn(ii, jj))+((Ln(ii, jj))*(1 + (w(jj)^2))))/(Jn(ii, jj)));
        det = (r*exp(-alfa*n(ii)*T/2))^4;
        rho = (exp(-alfa*n(ii)*T/2))/cap_w;

        coef_A = ((s_om_2*((2*(r^2)*(cap_w^2))+((Delta^2)*((1+r)^2)))) + (4*c_om_2*s_om_2*cap_w*r*Delta*(1+r)) + ((c_om_2^2)*2*(r^2)*(cap_w^2)));
        coef_B = ((2*Delta*((1+r)^2)*(w(jj)^2)*s_om_2) + (4*c_om_2*s_om_2*cap_w*r*(1+r)*(w(jj)^2)));
        coef_C = (((1+r)^2)*(w(jj)^4)*(s_om_2^2));

        P = +coef_A*rho;
        Q = 1 + det + (coef_B*rho);
        R = +coef_C*rho;

        coef_S = (Q^2) - (4*P*R);
        if coef_S >= 0
            yn_II_SN(jj, ii) = (-Q + sqrt(coef_S)) / (2 * P);
            yn_PF(jj, ii) = (-Q - sqrt(coef_S)) / (2 * P);
        else
            yn_II_SN(jj, ii) = NaN;
            yn_PF(jj, ii) = NaN;
        end

        % Calculate betaII_SN and betaII_PF
        betaII_SN(jj, ii) = (q(jj) * abs(Kn(ii, jj)))/(sqrt(((Kn(ii, jj))^2) + (((w(jj))^2)*((Ln(ii, jj))^2))));
        betaII_PF(jj, ii) = q(jj) * ((1 + (((((Kn(ii, jj)^2) + ((w(jj)^2)*(Ln(ii, jj)^2)))*(yn_PF(jj, ii)^2)) - (2*Ln(ii, jj)*Jn(ii, jj)*(w(jj)^2)*(yn_PF(jj, ii)))) / ((w(jj)^2)*(Jn(ii, jj)^2)))))^0.5;
    end
end

betaII_inf = betaI_inf;

% Plot Type I Motion results
figure()
for kk = 1:3 
    plot(w, betaI_PD(:,kk), 'LineWidth', 1.5, 'DisplayName', ['$$\beta_{PD}$$, n=' num2str(n(kk))])
    hold on
end
for kk = 1:3 
    plot(w, betaI_SN(:,kk), 'LineWidth', 1.5, 'DisplayName', ['$$\beta_{SN}$$, n=' num2str(n(kk))])
end
plot(w, betaI_inf(:,1), 'DisplayName', '$$\beta_{\infty}$$', 'LineWidth', 2, 'Color', 'k')

title('Type I Motion')
legend('show', 'Interpreter', 'latex','Location','southwest')
xlabel('$$\omega$$', 'Interpreter', 'latex');
ylabel('$$\beta$$', 'Interpreter', 'latex');
grid on
ylim([0 1])
xlim([0 5])
hold off

% Plot Type II Motion results
figure()
for kk = 1:length(n)
    plot(w, betaII_SN(:,kk), 'LineWidth', 1.5, 'DisplayName', ['$$\beta_{SN}$$, n=' num2str(n(kk))])
    hold on
end
for kk = 1:length(n)
    plot(w, betaII_PF(:,kk), 'LineWidth', 1.5, 'DisplayName', ['$$\beta_{PF}$$, n=' num2str(n(kk))])
end
plot(w, betaII_inf(:,1), 'DisplayName', '$$\beta_{\infty}$$', 'LineWidth', 2, 'Color', 'k')

title('Type II Motion')
legend('show', 'Interpreter', 'latex','Location','northwest')
xlabel('$$\omega$$', 'Interpreter', 'latex');
ylabel('$$\beta$$', 'Interpreter', 'latex');
grid on
ylim([0 2])
xlim([0 5])
hold off
