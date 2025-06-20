clc;
clear all;
close all;

%% BSE 
% dV/dt = (sigma)^2/2*d^2V/dx^2 + (r-sigma^2/2)dV/dx -rV
% E =100;
% t = [0,T-0.01] T=1
% x =[ln2,ln3/2]
% B.C V(x (min), t )= max(Ee^x(min)-E,0)=100
% B.C V(x (max), t )= max(Ee^x(max)-E,0)=50
% I.C V(x, 0 ) = max(Ee^x-E,0);

%tetta scheme : M u(t) = (sigma^2/2)Ku+(r-sigma^2/2)Nu-rMu 
x1 = log(0.5);
x2 = log(1.5);
Nx = 53;
dx = (x2 - x1) / Nx;
x = linspace(x1, x2, Nx+1);  
h = diff(x);

sigma = 0.2;
r= 0.02;

t1 = 0;
t2 = 2;
Nt =20;
dt = (t2 - t1) / Nt;
t = linspace(t1, t2, Nt+1);

%% Matrices C K M for only inner points 
K = zeros(Nx+1, Nx+1);
M = zeros(Nx+1, Nx+1);
C = zeros(Nx+1, Nx+1);

for i = 2:Nx
    %
    K(i,i) = -(1/dx + 1/dx);
    K(i,i+1) = 1/h(i);
    K(i,i-1) = 1/h(i-1);

    % 
    M(i,i) = (h(i-1) + h(i)) / 3;
    M(i,i+1) = h(i) / 6;
    M(i,i-1) = h(i-1) / 6;

    % 
    C(i,i-1) = +1/(2*dx);
    C(i,i+1) = -1/(2*dx);
end

%% Initial condition
E =100;

x_inner = x(2:end-1);  % Exclude boundaries
f = max((E * exp(x_inner) - E), 0);
u=f';
u_all = zeros(Nx-1, Nt+1);
u_all(:, 1) = u;

%% Boundary condition vectors
BC1 = max(E*exp(min(x))-E,0);
BC2 = max(E*exp(max(x))-E,0);

K = K(2:Nx,2:Nx);
M = M(2:Nx,2:Nx);
C = C(2:Nx,2:Nx);

BCK = zeros(Nx-1,1);
BCM = zeros(Nx-1,1);
BCC = zeros(Nx-1,1);

BCK(1,:) = 1/dx*BC1*(sigma^2/2);
BCK(end,:) = 1/dx*BC2*(sigma^2/2);

BCM(1,:) = r*(dx/6)*BC1;
BCM(end,:) = r*(dx/6)*BC2;

BCC(1,:) = 1/(2)*BC1*(r-sigma^2/2);
BCC(end,:) = -(r-sigma^2/2)*1/(2)*BC2;


%tetta scheme : M u(t) = (sigma^2/2)Ku+(r-sigma^2/2)Nu-rMu 
%% Tetta scheme
theta = 0.5;
A = M - theta*dt*((sigma^2)/2)*K -theta*dt*(r-((sigma^2)/2))*C +theta*dt*r*M;
B = M + (1-theta)*dt*((sigma^2)/2)*K +(1-theta)*dt*(r-((sigma^2)/2))*C -(1-theta)*dt*r*M;
for n = 1:Nt
    bc = BCK-BCM+BCC;
    b = (B * u)+bc*dt;
    
    u_next = A \ b;
    
    
   
    
    
    u = u_next;
    
   
    u_all(:, n+1) = u;
end
% Построение графика
figure;
[X, T_grid] = meshgrid(x(2:Nx), t); % Сетка (x, t)
surf(X, T_grid, u_all');  % u_all транспонируем, чтобы строки — время
xlabel('x');
ylabel('t');
zlabel('V(x,t)');
shading interp;
colorbar;

S_full = E * exp(x);  % Включая граничные точки
U_full = [BC1 * ones(1, Nt+1); u_all; BC2 * ones(1, Nt+1)];
disp(f);

% Построение графика в координатах S и t
figure;
S = E * exp(x(2:Nx));             % Преобразуем x в реальные значения актива S
[S_grid, T_grid] = meshgrid(S, t); 
surf(S_grid, T_grid, u_all');     % u_all транспонируем, чтобы строки — время
xlabel('S');
ylabel('t');
zlabel('V(S,t)');
title('Option Price V(S,t) in (S, t) coordinates');
shading interp;
colorbar;


K = 100;          
     
T_max = 1;        



S = linspace(1\2*K, 3\2*K, Nx + 1);       
T = linspace(0.01, T_max, Nt + 1);

[S, T] = meshgrid(S, T);

d1 = (log(S ./ K) + (r + 0.5 * sigma^2) .* T) ./ (sigma .* sqrt(T));
d2 = d1 - sigma .* sqrt(T);

price = S .* normCDF(d1) - K * exp(-r .* T) .* normCDF(d2);

figure;
surf(S, T, price);
xlabel('Stock Price S');
ylabel('Time to Maturity T');
zlabel('Option Price');

% Пользовательская функция для вычисления CDF стандартного нормального распределения
function cdf = normCDF(x)
    % Аппроксимация с использованием функции ошибок (error function)
    cdf = 0.5 * (1 + erf(x / sqrt(2)));
end


S_grid_analytical = E * exp(x(2:Nx));     % Только внутренние узлы
T_grid_analytical = t;                   % Временная сетка

[SS, TT] = meshgrid(S_grid_analytical, T_grid_analytical);   % сетка (S, T)
d1 = (log(SS / E) + (r + 0.5 * sigma^2) .* (t2 - TT)) ./ (sigma * sqrt(t2 - TT));
d2 = d1 - sigma * sqrt(t2 - TT);

d1(TT == t2) = 0;
d2(TT == t2) = 0;

V_exact = SS .* normCDF(d1) - E * exp(-r * (t2 - TT)) .* normCDF(d2);


V_numeric = u_all';   % u_all: строки — время, столбцы — x
error_matrix = abs(V_numeric - V_exact);