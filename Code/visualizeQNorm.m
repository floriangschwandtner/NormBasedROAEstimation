%% Tool to visualize QNorm


%Q = [3.000   -2.0000; -2.0000    3.0000];
Q = [93.3165606900217	3.70733512793459; 3.70733512793459	1.55233174728444]
%Q = eye(2);
figure
hold on
x=[0;0];
for phi=0.0:0.1:2*pi
    x = [cos(phi);sin(phi)];
    qnorm = QNorm(x,Q);
    x_plot = qnorm*x;
    plot(x_plot(1),x_plot(2), '.k');
end