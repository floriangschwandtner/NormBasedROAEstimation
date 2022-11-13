%% Program to determine lyapunov function and ROA of discrete system
clear all

%% System definition
syms xk yk

system = "vdp";
a=2;

if system=="vdp"
    f_x(1) = xk-yk;
    f_x(2) = xk+(1-a)*yk+a*xk^2*yk;
elseif system=="ellipse"
    f_x(1) = 1/2*xk*(1+xk^2+2*yk^2);
    f_x(2) = 1/2*yk*(1+xk^2+2*yk^2);
    A = [0.5, 0; 0, 0.5];
    g_x(1) = 1/2*xk*(xk^2+2*yk^2);
    g_x(2) = 1/2*yk*(xk^2+2*yk^2);
end

%% Visualize
%%findBallR(f_x,[xk yk]);

if system=="vdp"
    
    fun = @(x) sqrt(x(1)^2+x(2)^2);
    nonlcon = @boundaryConditionsVdP;
    x0 = [10 10];
    A = []; % No other constraints
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    tic
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
    toc
    R = norm(x);

elseif system=="ellipse"
    
    fun = @(x) sqrt(x(1)^2+x(2)^2);
    nonlcon = @boundaryConditionsEll;
    x0 = [10 10];
    A = []; % No other constraints
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
    
    R = norm(x);
end

%R=0.365;
%contour(X,Y,V_xnew,(p+1)*R^2);

%figure

if system=="vdp"
    %[TOUT,Yout] = ode45(@myvdp,[0 20],[-1.5 0.5]); %simulate van pol cycle
    %plot(Yout(:,1),Yout(:,2),'r')
    openfig("ROA_discrVdP_noFunHandles.fig")
elseif system=="ellipse"
    t = linspace(0,2*pi) ;
    a_ell = 1 ; 
    b_ell = 1/sqrt(2) ;
    x_ell = a_ell*cos(t);
    y_ell = b_ell*sin(t);
    plot(x_ell,y_ell,'r')
end

hold on
%% Calc lyapunov function

P_set = [0, 4];
eststrarr = strings(length(P_set));

for i=1:length(P_set)
    p = P_set(i);
    V_x = 0;
    
    f_xk = [xk yk];
    
    for k = 0:p
        %f_xk = f_x.^k;
        V_x=V_x+sqrt(f_xk(1)^2+f_xk(2)^2)^2;
        f_xk = subs(f_xk, [xk yk], [f_x(1) f_x(2)]);
    end

    X = -2:0.01:2;
    Y = transpose(X);
    V_xMatlab = matlabFunction(V_x);
    V_xnew = V_xMatlab(X,Y);

    %contour(X,Y,V_xnew,[1:6]);
    c=(p+1)*R^2;
    contour(X,Y,V_xnew,[c c], '-', 'LineWidth',2);
    eststr = sprintf(['ROA Estimation p=',num2str(p),', c=',num2str(c)]);
    eststrarr(i+1) = eststr;
end

if system=="vdp"
    axis([-1.5 1.5 -2.0 2.0])
elseif system=="ellipse"
    axis([-1 1 -1 1])
end


xlabel("x_k")
ylabel("y_k")
title("Euclidean Norm based approach");

legend(eststrarr)
%figure
%surf(X,Y,V_xnew)


function [c,ceq] = boundaryConditionsEll(x)
    ceq = -sqrt((1/2*x(1)*(1+x(1)^2+2*x(2)^2))^2+(1/2*x(2)*(1+x(1)^2+2*x(2)^2))^2)+sqrt(x(1)^2+x(2)^2);
    c = [];
end


function dydt = myvdp(t,y)
    %mycdp  Evaluate the van der Pol ODEs for a = 2
    dydt = [y(2); 2*(1-y(1)^2)*y(2)-y(1)];
end
