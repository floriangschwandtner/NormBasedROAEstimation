function [R, x_optD] = NN_casadi_calcR(P)
%NN_CASADI_CALCR Summary of this function goes here
%   Detailed explanation goes here

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

% A = zeros(2,2);
% A(1,1) = to_double(A_res(1,1));
% A(1,2) = to_double(A_res(1,2));
% A(2,1) = to_double(A_res(2,1));
% A(2,2) = to_double(A_res(2,2));
A = A_res.full();
%P = eye(2);
%P = [[1.0, 0.0];[0.0, 1.0]];
Q = dlyap(A', P);


%% Optimization
f = sqrt(transpose(x)*Q*x);
g = -sqrt((transpose(f_x)*Q*f_x))+sqrt(transpose(x)*Q*x);
nlp = struct('x',x, 'f',f, 'g',g);

S = nlpsol('S', 'ipopt', nlp);

r = S('x0',[10.0 -10.0],'lbg',0);
x_opt = r.x;
disp(x_opt)
x_optD = x_opt.full();

R = x_optD'*Q*x_optD;
end

