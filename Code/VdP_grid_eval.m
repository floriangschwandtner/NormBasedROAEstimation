clear all

%% System definition
syms xk yk


%% parameters
a=2;

f_x(1) = xk-yk;
f_x(2) = xk+(1-a)*yk+a*xk^2*yk;

f_x = transpose(f_x);

J = jacobian(f_x);
A = subs(J,[xk yk], [0.0 0.0]);

%P = eye(2);
P = [[1.0, 0.0];[0.0, 1.0]];
Q = dlyap(A', P); % A' becuase lyapunov eq. is defined differently in dlyap()
%Q = eye(2);

%% R finden;
%fun = @(x) sqrt(x(1)^2+x(2)^2);
fun = @(x) sqrt(x*Q*transpose(x));
nonlcon = @(x) boundaryConditionsQNorm(x,transpose(f_x),Q);
%nonlcon = @boundaryConditionsVdP;
x0 = [10 10];
A_constr = []; % No other constraints
b_constr = [];
Aeq_constr = [];
beq_constr = [];
lb = [];
ub = [];
tic
x_opt = fmincon(fun,x0,A_constr,b_constr,Aeq_constr,beq_constr,lb,ub,nonlcon);
toc

R = sqrt(x_opt*Q*x_opt');


%% Visualize

n_grid = 800;

x_lim = 2.0;
y_lim = 1.5;
x_vec = linspace(-x_lim,x_lim, n_grid);
y_vec = linspace(-y_lim, y_lim, n_grid);

openfig("ROA_discrVdP_noFunHandles.fig");
hold on
% xfun = @(xf,yf) QNorm([xf;yf],Q);
% fxfun = @(xf,yf) QNorm(subs(f_x, [xk; yk], [xf; yf]),Q);
% 
% [X,Y]=ndgrid(x_vec,y_vec);
% 
% xfun(X,Y)

f_xM = matlabFunction(f_x);

Z = zeros(n_grid, n_grid);
Z_val = zeros(n_grid,n_grid);

for k = 1:n_grid
    for j = 1:n_grid
        if(QNorm([x_vec(k);y_vec(j)], Q)>=QNorm(f_xM(x_vec(k),y_vec(j)),Q))
            Z(k,j) = 1;
            %Z_val(k,j) = QNorm([x_vec(k);y_vec(j)], Q)-QNorm(f_xM(x_vec(k),y_vec(j)),Q);
        end
        Z_val(k,j) = QNorm([x_vec(k);y_vec(j)], Q);

        %         res = f_xM(x_vec(k),y_vec(j));
%         Z(k,j)= sqrt(x_vec(k)^2+y_vec(j)^2)-sqrt(res(1)^2+res(2)^2);
%         if(sqrt(x_vec(k)^2+y_vec(j)^2)>sqrt(res(1)^2+res(2)^2))
%             Z(k,j) = 1;
%         end
    end
end

surf(x_vec, y_vec, Z,  'edgecolor','none','FaceAlpha',0.4)
plot(x_opt(1), x_opt(2), 'g*', 'LineWidth',4)

figure
%surf(x_vec, y_vec, Z_val, 'edgecolor','none');
%zlim([-1, inf])
contour(x_vec, y_vec, Z_val,[1.0 1.0])

%             plot(x_vec(k),y_vec(j),'g*')
%         else
%             plot(x_vec(k),y_vec(j),'r*')

