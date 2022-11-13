clear

%% import and set up casadi
addpath('/Users/floriangschwandtner/repos/casadi-osx-matlabR2015a-v3.5.5')
import casadi.*

%opti = casadi.Opti();

%xk = opti.variable();
%yk = opti.variable();
xk = SX.sym('xk',1);
yk = SX.sym('yk',1);
x = [xk; yk];

%% Load NN

% parameters
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
v1eq = W1*x + b1;
w1eq = tanh(v1eq);
v2eq = W2*w1eq + b2;
w2eq = tanh(v2eq);
u = W3*w2eq + b3;

%u = sign(u)*min(abs(u),0.7);
%u = 0.7*tanh(u);
%u = max(-0.7, u);
%u = min(u ,0.7);
% closed loop representation
f_x = AG*[xk;yk]+BG1*q+BG2*u;

J = jacobian(f_x,x);
Jfun = Function('Jfun', {xk,yk}, {J});
A_res = Jfun(0.0,0.0);

%A = subs(J,[xk yk], [0.0 0.0]);

A = A_res.full();

%P = eye(2);
P = [[1.0, 0.0];[0.0, 1.0]];
Q = dlyap(A', P);


%% Optimization
% f = transpose(x)*Q*x;
% g = (transpose(f_x)*Q*f_x)-transpose(x)*Q*x;
% nlp = struct('x',vertcat(xk,yk), 'f',f, 'g',g);
% 
% S = nlpsol('S', 'ipopt', nlp);
% 
% r = S('x0',[0,0],'lbg',0,'ubg',0);
% x_opt = r.x;
% disp(x_opt)
% x_optD = zeros(2,1);
% x_optD(1) = to_double(x_opt(1));
% x_optD(2) = to_double(x_opt(2));
% 
% R = x_optD'*Q*x_optD;

[R, x_optD] = NN_casadi_calcR(P);

%%figure
openfig("NN_invertedPendulum_ROA.fig");

hold on

%% reload NN in symbolic form to evaluate
clear x xk yk q p f_x v1eq w1eq v2eq w2eq

syms xk yk
x = [xk; yk];
q = xk-sin(xk);
v1eq = W1*x + b1;
w1eq = tanh(v1eq);
v2eq = W2*w1eq + b2;
w2eq = tanh(v2eq);
u = W3*w2eq + b3;
u = 0.7*tanh(u);

f_x = AG*x+BG1*q+BG2*u;

%% Calc lyapunov function
P_set = [0];
eststrarr = strings(2*size(P_set));

for i=1:length(P_set)
    p = P_set(i);
    V_x = 0;

%     f_xk = [xk; yk];
% 
%     for k = 0:p
%         %f_xk = f_x.^k;
%         V_x=V_x+sqrt(transpose(f_xk)*Q*f_xk)
%         f_xk = subs(f_xk, [xk yk], [f_x(1) f_x(2)]);
%     end
    gridsize = 100;
    X = linspace(-5,5,gridsize);
    Y = linspace(-12,12,gridsize);


%     V_xMatlab = matlabFunction(V_x);
%     V_xnew = V_xMatlab(X,Y);
    
    V_xnew = evaluateSymQNorm(f_x,X,Y,Q,p);
    M_xp = evaluate_Mx_grid(f_x, X,Y,Q,p);

    %contour(X,Y,V_xnew,[1:6]);
    c=(p+1)*R^2;
    contour(X,Y,V_xnew,[c c], '-k', 'LineWidth', 2);
    contour(X,Y,M_xp,[R^2 R^2], '-b', 'LineWidth', 2);
    eststr = sprintf(['Vx Estimation p=',num2str(p),', c=',num2str(c)]);
    eststrarr(i) = eststr;
    eststr = sprintf(['Mp Estimation p=',num2str(p),', c=',num2str(c)]);
    eststrarr(i+1) = eststr;

end

legend(eststrarr)

