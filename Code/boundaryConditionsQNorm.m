function [c,ceq] = boundaryConditionsQNorm(x, f_x, Q)
    syms xk yk;
    ceq = -sqrt(f_x*Q*transpose(f_x))+sqrt(x*Q*transpose(x));
    ceq = subs(ceq, [xk yk], [x(1) x(2)]);
    ceq = double(ceq);
    c = [];
end
