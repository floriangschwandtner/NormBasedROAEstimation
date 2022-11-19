clear all

%% System definition
syms xk yk


%% parameters
g = 10; % gravitational coefficient
m = 0.15; % mass
l = 0.5; % length
mu = 0.05; % frictional coefficient
dt = 0.02; % sampling period

% System from yin et al.
% xk entspricht theta
% yk entspricht theta_dot
AG = [1,      dt;...
    g/l*dt, 1-mu/(m*l^2)*dt];
% describes how q enters the system
BG1 = [0;...
    -g/l*dt];
% describes how u enters the system
BG2 = [0;...
    dt/(m*l^2)];

fname = 'Wb_s32_tanh/';
load([fname 'W1.csv'])
load([fname 'W2.csv'])
load([fname 'W3.csv'])

q = xk-sin(xk);
% build u as a NN
% all b are 0
n1 = size(W1,1);
n2 = size(W2,1);
n3 = size(W3,1);
nphi = n1+n2;

b1 = zeros(n1,1);
b2 = zeros(n2,1);
b3 = zeros(n3,1);

% build u from NN
x = [xk; yk];
v1eq = W1*x + b1;
w1eq = tanh(v1eq);
v2eq = W2*w1eq + b2;
w2eq = tanh(v2eq);
u = W3*w2eq + b3;

% u = sign(u)*min(abs(u),0.7);
u = 0.7*tanh(u);
%u = max(-0.7, u);
%u = min(u ,0.7);
% closed loop representation
f_x = AG*x+BG1*q+BG2*u;
J = jacobian(f_x);
A = subs(J,[xk yk], [0.0 0.0]);

%P = eye(2);
P = [[1.0, 0.0];[0.0, 1.0]];
Q = dlyap(A', P); % A' becuase lyapunov eq. is defined differently in dlyap()
%Q = eye(2);

%% Visualize
% fun = @(x) sqrt(x*Q*transpose(x));
% 
% nonlcon = @(x) boundaryConditionsQNorm(x,transpose(f_x),Q);
% x0 = [10 10];
% A_constr = []; % No other constraints
% b_constr = [];
% Aeq = [];
% beq = [];
% lb = [];
% ub = [];
% 
% 
% options = optimoptions(@fmincon,'Display','iter');
% 
% x_opt = fmincon(fun,x0,A_constr,b_constr,Aeq,beq,lb,ub,nonlcon, options);
% 
% R =  sqrt(x_opt*Q*x_opt');
[R, x_opt] = NN_casadi_calcR(P);

n_grid = 400;

x_lim = 5.0;
y_lim = 25.0;
x_vec = linspace(-x_lim,x_lim, n_grid);
y_vec = linspace(-y_lim, y_lim, n_grid);

openfig("NN_invertedPendulum_ROA.fig");
hold on

% xfun = @(xf,yf) QNorm([xf;yf],Q);
% fxfun = @(xf,yf) QNorm(subs(f_x, [xk; yk], [xf; yf]),Q);
% 
% [X,Y]=ndgrid(x_vec,y_vec);
% 
% xfun(X,Y)

f_xM = matlabFunction(f_x);

Z = zeros(n_grid, n_grid);

for k = 1:n_grid
    for j = 1:n_grid
        if(QNorm([x_vec(k);y_vec(j)], Q)>=QNorm(f_xM(x_vec(k),y_vec(j)),Q))
            Z(k,j) = 1;
        end
    end
end

surf(x_vec, y_vec, Z,  'edgecolor','none','FaceAlpha',0.4)

plot(x_opt(1), x_opt(2), 'g*', 'LineWidth',4)



%             plot(x_vec(k),y_vec(j),'g*')
%         else
%             plot(x_vec(k),y_vec(j),'r*')

