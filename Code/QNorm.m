function qnorm = QNorm(x,Q)
%QNORM calc Q Norm of nx1 vector with nxn Q matrix.
%   Detailed explanation goes here
    qnorm = sqrt(x'*Q*x);
end

