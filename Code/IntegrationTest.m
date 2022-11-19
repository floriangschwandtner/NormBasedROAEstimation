%% Test 

clear all
syms xk yk;
a=2.0;
f_x(1) = xk-yk;
f_x(2) = xk+(1-a)*yk+a*xk^2*yk;

figure
hold on

x_start=-1.2;
x_end=1.2;
y_start=-2.0;
y_end=2.0;

for k=x_start:0.05:x_end
    for j=y_start:0.1:y_end
        k_vec = 1:1:20;
        x_vec = zeros(length(k_vec));
        y_vec = zeros(length(k_vec));
        
        x_vec(1) = k;
        y_vec(1) = j;
        
        for i=1:length(k_vec)-1
            x_vec(i+1) = subs(f_x(1), [xk yk], [x_vec(i) y_vec(i)]);
            y_vec(i+1) = subs(f_x(2), [xk yk], [x_vec(i) y_vec(i)]);
        end
        
        if(abs(x_vec(i)^2+y_vec(i)^2)<10^3)
            plot(x_vec(1), y_vec(1), '.r', 'HandleVisibility','off')
        else
            %plot(x_vec(1), y_vec(1), '.k')
        end
    end
end

axis([x_start x_end y_start y_end])