function [c,ceq] = boundaryConditionsVdP(x)
    a=2;
    
    P = [2 3];
    ceq = zeros(length(P));
    %p_f=3;

    syms xk yk;
    f_x(1) = xk-yk;
    f_x(2) = xk+(1-a)*yk+a*xk^2*yk;

    for j = 1:length(P) 
        f_xk = [xk yk];

        for k = 1:P(j)
            f_xk = subs(f_xk, [xk yk], [f_x(1) f_x(2)]);
        end
        f_p = matlabFunction(f_xk);

        ceq(j) = -norm(f_p(x(1), x(2)))+sqrt(x(1)^2+x(2)^2);
    end
    
    c = [];
end
