%% Program to determine lyapunov function and ROA of discrete system
clear all

%% System definition
syms xk yk

system = "NN";

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
elseif system=="NN"
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
    J_00 = subs(J,[xk yk], [0.0 0.0]);
    J_abs = abs(J_00);
end



%% Visualize

%%findBallR(f_x,[xk yk]);

if system=="vdp"
    
    fun = @(x) sqrt(x(1)^2+x(2)^2);
    nonlcon = @boundaryConditionsVdPQNorm;
    x0 = [10 10];
    A = []; % No other constraints
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
    
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

elseif system=="NN"
    fun = @(x) sqrt(x(1)^2+x(2)^2);
    nonlcon = @boundaryConditionsNN;
    x0 = [10 10];
    A = []; % No other constraints
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    options = optimoptions(@fmincon,'Display','iter');

    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
    
    R = norm(x);

end

%R=0.365;
%contour(X,Y,V_xnew,(p+1)*R^2);

figure

if system=="vdp"
    [TOUT,Yout] = ode45(@myvdp,[0 20],[-1.5 0.5]); %simulate van pol cycle
    plot(Yout(:,1),Yout(:,2),'r')
elseif system=="ellipse"
    t = linspace(0,2*pi) ;
    a_ell = 1 ; 
    b_ell = 1/sqrt(2) ;
    x_ell = a_ell*cos(t);
    y_ell = b_ell*sin(t);
    plot(x_ell,y_ell,'r')
elseif system =="NN"
    load("NN_invertedPendulum_ROA.fig");
end

hold on
%% Calc lyapunov function
eststrarr = ["ROA"];
P_set = [0 1 8];

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
    contour(X,Y,V_xnew,[c c], '-b');
    eststr = sprintf(['ROA Estimation p=',num2str(p),', c=',num2str(c)]);
    eststrarr(i+1) = eststr;
end

if system=="vdp"
    axis([-3 3 -3 3])

elseif system=="ellipse"
    axis([-1 1 -1 1])
end


xlabel("x_k")
ylabel("y_k")

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
