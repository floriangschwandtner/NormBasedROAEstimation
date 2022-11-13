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
%P = double(eye(2));
P = [[1.0, 0.0];[0.0, 1.0]];
Q = dlyap(A', P);

%% R finden;
%fun = @(x) sqrt(x(1)^2+x(2)^2);
[R, x_opt] = VdP_casadi_calcR(P);

% figure
% 
% [TOUT,Yout] = ode45(@myvdp,[0 20],[-0.2 0.2]); %simulate van pol cycle
% plot(Yout(:,1),Yout(:,2),'r')
openfig("ROA_discrVdP_noFunHandles.fig");
hold on
%% Calc lyapunov function
P_set = [4];
eststrarr = strings(size(P_set));

% transpose to make 2x1
f_x = transpose(f_x);

for i=1:length(P_set)
    p = P_set(i);
    V_x = 0;
    
    gridsize = 200;

    X = linspace(-2,2,gridsize);
    Y = linspace(-2,2,gridsize);



%     f_xk = [xk; yk];
%     
%     for k = 0:p
%         %f_xk = f_x.^k;
%         V_x=V_x+transpose(f_xk)*Q*f_xk;
%         f_xk = subs(f_xk, [xk yk], [f_x(1) f_x(2)]);
%     end
    
    %V_xMatlab = matlabFunction(V_x);
    %V_xnew = V_xMatlab(X,Y');

    V_xtest = evaluateSymQNorm(f_x,X,Y,Q,p);

    %contour(X,Y,V_xnew,[1:6]);
    c=(p+1)*R^2;
    %contour(X,Y,V_xnew,[c c], '-k', 'LineWidth',3);
    %eststr = sprintf(['QNorm ROA Estimation p=',num2str(p),', c=',num2str(c)]);
    %eststrarr(end+1) = eststr;

    contour(X,Y,V_xtest,[c c],'-', 'LineWidth',2);
 
    eststr = sprintf(['QNorm ROA Estimation p=',num2str(p),', c=',num2str(c)]);
    eststrarr(i) = eststr;
end

axis([-1.5 1.5 -2.0 2.0])

xlabel("x_k")
ylabel("y_k")
title("Q-Norm based approach");

legend(eststrarr)

function dydt = myvdp(t,y)
%mycdp  Evaluate the van der Pol ODEs for a = 2
dydt = [y(2); 2*(1-y(1)^2)*y(2)-y(1)];
%dydt = [y(1)-y(2); y(1)-y(2)+2*y(1)^2*y(2)];

end
