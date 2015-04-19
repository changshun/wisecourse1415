
% EXAMPL0.M:
% Solving the stochastic neoclassical growth model with the recursive law
% of motion by hand, and implement the impulse responses.

%
clear all;

disp('EXAMPLE 0: The stochastic neoclassical growth model');


disp('Hit any key when ready...');
pause;


% Setting parameters:

Z_bar     = 1;    % Normalization
rho       = .36;  % Capital share
delta     = 0.025; % Depreciation rate for capital
betta     = .99;  % Discount factor beta 
eta       = 1;  % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)
psi       = .95;  % autocorrelation of technology shock
sigma_eps = .712; % Standard deviation of technology shock.  Units: Percent.

% Calculating the steady state:

R_bar   = 1.0/betta;  % One percent real interest per quarter
K_bar   = ((rho*Z_bar)/(R_bar - 1 + delta))^(1.0/(1 - rho));
Y_bar   = Z_bar*K_bar^rho;
C_bar   = Y_bar - delta*K_bar;


% Calculating the intermediate parameters:
v_rk = -(1-rho)*(1-betta*(1-delta));
v_rz = 1-betta*(1-delta);

% Define parameters for solving eta_kk:
% v_kk^2 - gamma * v_kk + R_bar = 0;
gamma = -v_rk*C_bar/(K_bar*eta) + 1 + R_bar;

v_kk = gamma/2 - ((gamma/2)^2 - R_bar)^0.5;

v_ck = (R_bar - v_kk)*K_bar/C_bar;

v_kz = (v_rz*psi + eta*(1-psi)*Y_bar/C_bar)/(-v_rk + eta*v_ck + eta*(1-psi)*K_bar/C_bar);
v_cz = Y_bar/C_bar - v_kz*K_bar/C_bar;



%% EX1. Calculate impulse responses to a shock to technology.
TT = -1:1:20;
K = zeros(size(TT,2),1); Z = zeros(size(TT,2),1); C = Z; Y = Z; R = Z;
Z(2) = 1; % Shock to Z at t=1;

for t = TT(2):TT(end);
    i = t+2;
    Z(i)=psi^(t);
    K(i) = v_kk*K(i-1)+v_kz*Z(i);
    C(i) = v_ck*K(i-1)+v_cz*Z(i);
    Y(i) = rho*K(i-1) + Z(i);
    R(i) = v_rk*K(i-1) +v_rz*Z(i);
end;

 
figure
plot(TT,K,TT,Z,TT,C,TT,Y,TT,R,TT,TT*0,'linewidth',2);
legend('Capital','Technology','Consumption','Output','Return');
xlim([-1 20])
%ylim([-0.2 1])
xlabel('Time');
ylabel('Percentage deviation from steady state');

%% EX2. Calculate impulse responses to an initial deviation to capital.
TT = -1:1:20;
K = zeros(size(TT,2),1); Z = zeros(size(TT,2),1); C = Z; Y = Z; R = Z;
K(1) = 0.1;
for t = TT(2):TT(end);
    i = t+2;
    K(i) = v_kk*K(i-1)+v_kz*Z(i);
    C(i) = v_ck*K(i-1)+v_cz*Z(i);
    Y(i) = rho*K(i-1) + Z(i);
    R(i) = v_rk*K(i-1) +v_rz*Z(i);
end;
 
figure
plot(TT,K,TT,Z,TT,C,TT,Y,TT,R,TT,TT*0,'linewidth',2);
legend('Capital','Technology','Consumption','Output','Return');
xlim([-1 TT(end)])
xlabel('Time');
ylabel('Percent deviation from steady state');


%%
A = [v_kk v_kz*psi; 0 psi]; 
B = [0 v_kz; 0 1];
sigma = B*[0 0; 0 sigma_eps]*B';
VecGamma_0 =inv(eye(4) - kron(A, A))*reshape(sigma, 4,1);
Gamma_0 = reshape(VecGamma_0, 2,2);
Gamma_1 = A*Gamma_0;
<center>页面执行时间：15.625 毫秒</center>