clear all; close all; clc;

% Caricamento dei Dati
T = readtable('data.csv'); % le colonne devono essere: time,T_out,G,T_in 
Tout = T.T_out;
G = T.G;
Tin = T.T_in;
N = length(Tin);

% Modello Lineare: Tin = a + b*Tout + c*G
X = [ones(N,1), Tout, G];
theta = (X'*X)\(X'*Tin);
a = theta(1); b = theta(2); c = theta(3);
fprintf('Equazione Modello Lineare: Tin = %.3f + %.3f*Tout + %.6f*G\n', a,b,c);

% Modello Dinamico: Tin(t+1) = gamma*Tin(t) + beta1*Tout(t) + beta2*G(t)
Y = Tin(2:end);
Phi = [Tin(1:end-1), Tout(1:end-1), G(1:end-1)];
phi_theta = (Phi'*Phi)\(Phi'*Y);
gamma = phi_theta(1);
beta1 = phi_theta(2);
beta2 = phi_theta(3);
fprintf('Equazione Modello Dinamico: Tin(t+1)=%.3f*Tin(t) + %.3f*Tout(t) + %.6f*G(t)\n', gamma,beta1,beta2);

% Scarto Residuo ed Errore Medio
Yhat_static = X*theta;
rmse_static = sqrt(mean((Tin - Yhat_static).^2));
Yhat_dyn = [Tin(1); gamma*Tin(1:end-1) + beta1*Tout(1:end-1) + beta2*G(1:end-1)];
rmse_dyn = sqrt(mean((Tin - Yhat_dyn).^2));
fprintf('RMSE Modello Lineare: %.3f ; RMSE Modello Dinamico: %.3f\n', rmse_static, rmse_dyn);

% Fattori di Scala serra
A_f_ref = 61.654;    % [m^2] Area Irradiata di Riferimento (escluso il terreno)
A_ref = 79.564;      % [m^2] Area di dispersione totale di Riferimento (incluso il terreno)
V_ref = 45.9;        % [m^3] Volume interno di Riferimento
A_f = 154;            % [m^2] Area Irradiata della Serra (esluso il terreno)
A_tot = 101;         % [m^2] Area di dispersione totale della Serra (incluso il terreno)
V_tot = 61.2;        % [m^3] Volume interno della Serra
U_ref = 5.8;         % [W/m^2K] Trasmittanza di riferimento
U = 2.8;             % [W/m^2K] Trasmittanza della Serra

% Calcolo dei Rapporti di Scala
Av_ref = A_ref / V_ref;
Av = A_tot / V_tot;
scale_b = (U / U_ref) * (Av / Av_ref);                                    % per il Coefficiente di Tout
scale_c = A_f / A_f_ref;                                                  % per Coefficiente di Irradianza
b_scaled = b * scale_b;
c_scaled = c * scale_c;
beta1_scaled = beta1 * scale_b;
beta2_scaled = beta2 * scale_c;
a_scaled = mean(Tin) - b_scaled*mean(Tout) - c_scaled*mean(G);            % per Intercetta di Tin (Modello Lineare)
a_dyn_scaled = mean(Tin(2:end)) - gamma*mean(Tin(1:end-1)) ...            % per Intercetta di Tin (Modello Dinamico)
                                - beta1_scaled*mean(Tout(1:end-1)) ...
                                - beta2_scaled*mean(G(1:end-1));

fprintf('Modello Lineare Scalato:\n');
fprintf('Tin = %.3f + %.3f*Tout + %.6f*G\n', a_scaled, b_scaled, c_scaled);
fprintf('Modello Dinamico Scalato:\n');
fprintf('Tin(t+1) = %.3f + %.3f*Tin(t) + %.3f*Tout(t) + %.6f*G(t)\n', ...
        a_dyn_scaled, gamma, beta1_scaled, beta2_scaled);

fprintf('Fattore scala per Tout = %.3f\n', scale_b);
fprintf('Fattore scala per G    = %.3f\n', scale_c);

% --- Plot Temperatura Interna della Serra ---
figure;
plot(T.time, Tin, 'b-', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Temperatura interna T_{in} [°C]');
title('Andamento della Temperatura Interna di Riferimento');  % come verifica del codice
grid on;
