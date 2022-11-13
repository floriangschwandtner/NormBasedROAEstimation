function [V_xp] = evaluateSymQNorm(f_x,X,Y,Q,p)
%EVALUATESYMQNORM evaluates a symbolic function f_x p times and returns the
%lyapunv function made by the q norm
%   f_x symbolic system function
%   X row Vector with values in X Range
%   Y row Vector with values in Y Range
%   Q Matrix for Q Norm. Simplifies to euclidean norm if Q=I
%   p degree of Lyapunov function

%V_xp = zeros(length(X),length(Y));

% convert to matlab function to speed up evaluation
f_xM = matlabFunction(f_x);

syms xk yk;
x = [xk; yk];
xQ = transpose(x)*Q*x;
xQM = matlabFunction(xQ);

%V_xp = xQM(X,Y');
dimx = length(X);
dimy = length(Y);

V_xp = zeros(dimx, dimy);
[Xg, Yg] = meshgrid(X,Y);

for k = 0:p
    % for each gridpoint calc xQ and next step
    for m=1:dimx
        for j=1:dimy
            V_xp(m,j)=V_xp(m,j)+xQM(Xg(m,j),Yg(m,j));
            fxmyj = f_xM(Xg(m,j),Yg(m,j));
            Xg(m,j) = fxmyj(1);
            Yg(m,j) = fxmyj(2);
        end
    end
end

end