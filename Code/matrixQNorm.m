function [QNormA] = matrixQNorm(A,Q)
    %MATRIXQNORM calculates the Q Norm of A w.r.t Q
    %   A is nxn matrix
    %   Q is p.d. nxn matrix
    
    fun = @(x) -sqrt(transpose(A*transpose(x))*Q*(A*transpose(x)));
    
    nonlcon = @(x) matrixQNormBC(x, Q);
    x0 = [10 10];
    A_constr = []; % No other constraints
    b_constr = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    
    
    options = optimoptions(@fmincon,'Display','iter');
    
    x = fmincon(fun,x0,A_constr,b_constr,Aeq,beq,lb,ub,nonlcon, options);
    
    QNormA = QNorm(A*transpose(x),Q);
end

function [c, ceq] = matrixQNormBC(x,Q)
    ceq = QNorm(transpose(x),Q)-1;
    c = [];
end