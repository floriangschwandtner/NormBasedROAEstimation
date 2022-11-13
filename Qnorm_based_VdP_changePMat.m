%% Program to determine lyapunov function and ROA of discrete system
clear all

%% System definition
syms xk yk

a=2;

f_x(1) = xk-yk;
f_x(2) = xk+(1-a)*yk+a*xk^2*yk;

%f_x(1) = yk; 
%f_x(2) = 2*(1-xk^2)*yk-xk;

%% Determine Q Matrix
J = jacobian(f_x);
A = subs(J, [xk;yk], [0.0;0.0]);
%g_x = transpose(f_x)-A*[xk; yk];

% Calculate Q
%P = [[0.10,1.0];[1.0,0.1]];
% figure
% 
% [TOUT,Yout] = ode45(@myvdp,[0 20],[-0.2 0.2]); %simulate van pol cycle
% plot(Yout(:,1),Yout(:,2),'r')
openfig("ROA_discrVdP_noFunHandles.fig");
hold on
%% Calc lyapunov function
%P_set = {eye(2);[[1.0, 1.0];[1.0, 0.0]]}; %;[[1.0, 3.0];[1.0, 0.0]]};
%P_set = {eye(2);10.0*eye(2)}; %;[[1.0, 3.0];[1.0, 0.0]]};
%P_set = {eye(2);[[10.0, 0.0];[0.0, 1.0]];[[0.1, 0.0];[0.0, 1.0]]}; 
%P_string = {'I';'[[10.0, 0.0];[0.0, 1.0]]';'[[0.1, 0.0];[0.0, 1.0]]'}; %;'[[1.0, 3.0];[1.0, 0.0]]'};
P_set = {eye(2)}; 
P_string = {'I'};
% transpose to make 2x1A'
f_xopt = f_x;
f_x = transpose(f_x);
eststrarr = strings(size(P_set));
for i=1:length(P_set)
    
    P = P_set{i};
    Q = dlyap(A', P)
    
    %% R finden;
    %fun = @(x) sqrt(x(1)^2+x(2)^2);
    fun = @(x) sqrt(x*Q*transpose(x));
    nonlcon = @(x) boundaryConditionsQNorm(x,f_xopt,Q);
    x0 = [10 10];
    A_constr = []; % No other constraints
    b_constr = [];
    Aeq_constr = [];
    beq_constr = [];
    lb = [];
    ub = [];
    options = optimoptions(@fmincon,'Display','iter');
    tic
    x = fmincon(fun,x0,A_constr,b_constr,Aeq_constr,beq_constr,lb,ub,nonlcon, options)
    toc
    
    R = sqrt(x*Q*x');

    p = 0;
    V_x = 0;
    
    gridsize = 100;

    X = linspace(-2,2,gridsize);
    Y = linspace(-2,2,gridsize);
    
    %V_xMatlab = matlabFunction(V_x);
    %V_xnew = V_xMatlab(X,Y');

    V_xtest = evaluateSymQNorm(f_x,X,Y,Q,p);

    %contour(X,Y,V_xnew,[1:6]);
    c=(p+1)*R^2;
    %contour(X,Y,V_xnew,[c c], '-k', 'LineWidth',3);
    %eststr = sprintf(['QNorm ROA Estimation p=',num2str(p),', c=',num2str(c)]);
    %eststrarr(end+1) = eststr;

    h = contour(X,Y,V_xtest,[c c], '-', 'LineWidth',2);
    %legend(h, 'QNorm ROA Estimation');
    eststr = sprintf(['QNorm ROA Estimation P=', P_string{i}]);
    eststrarr(i) = eststr;
end

axis([-1.5 1.5 -2.0 2.0])

xlabel("x_k")
ylabel("y_k")
title("Q-Norm based approach");

legend(eststrarr)
