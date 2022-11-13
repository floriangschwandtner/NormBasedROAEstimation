function [R, x_optD] = VdP_casadi_calcR(P)
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

% closed loop representation
a=2;
f_x = [xk-yk;xk+(1-a)*yk+a*xk^2*yk]

%f_x(1) = xk-yk;
%f_x(2) = xk+(1-a)*yk+a*xk^2*yk;

%f_x = transpose(f_x);

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
f = -sqrt(transpose(x)*Q*x);
g = -sqrt((transpose(f_x)*Q*f_x))+sqrt(transpose(x)*Q*x);
nlp = struct('x',x, 'f',f, 'g',g);

S = nlpsol('S', 'ipopt', nlp);

r = S('x0',[0.0 0.0],'ubg',0,'lbx',0);
x_opt = r.x;
disp(x_opt)
x_optD = x_opt.full();

R = x_optD'*Q*x_optD;
end

