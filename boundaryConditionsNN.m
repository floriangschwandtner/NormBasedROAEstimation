function [c,ceq] = boundaryConditionsNN(x, Q)
    a=2;
    
    P = [2 3];
    ceq = zeros(length(P));
    %p_f=3;

    syms xk yk;
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

    for j = 1:length(P) 
        f_xk = [xk yk];

        for k = 1:P(j)
            f_xk = subs(f_xk, [xk yk], [f_x(1) f_x(2)])
        end
        f_p = matlabFunction(f_xk);

        ceq(j) = -norm(f_p(x(1), x(2)))+sqrt(x(1)^2+x(2)^2);
    end
    
    c = [];
end
