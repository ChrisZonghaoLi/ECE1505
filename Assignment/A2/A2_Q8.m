close all
clear all

p_mean = [1.1 1.35 1.25 1.05];
% sigma = [0.2 -0.2 -0.12 0.02;
%          -0.2 1.4 0.02 0;
%          -0.12 0.02 1 -0.4;
%          0.02 0 -0.4 0.02];

sigma = [0.2 -0.2 -0.12 0.02;
         -0.2 1.4 0.02 0;
         -0.12 0.02 1 -0.4;
         0.02 0 -0.4 0.02];

H = sigma;
f = [0;0;0;0]; 
A = [-1 0 0 0;      % -x1 <= 0
     0 -1 0 0;      % -x2 <= 0
     0 0 -1 0;      % -x3 <= 0
     0 0 0 -1;      % -x4 <= 0
     -1.*p_mean];   % trans(-p_mean)x <= -r_min
Aeq = [1 1 1 1];
beq = 1;
 
for i = 1:1:15
    b(:,i) = [0;0;0;0;-0.1*i]; % -r_min = 0.1*i
    [x(:,i), fval(:,i), exitflag(:,i), output, lambda] = ...
    quadprog(H, f, A, b(:,i), Aeq, beq);
end

% x, fval, exitflag
exp_return = p_mean * x;

figure(1)
p1 = plot(exp_return, fval);
hold on
p2 = plot([1.3 1.3], [-0.5 1], 'r');
set(gca,'linewidth',2)
set(p1, 'linewidth', 3)
set(p2, 'linewidth', 3)
xlabel('expected return')
ylabel('risk of profolio')
legend('expected return vs risk of profolio')
grid
hold off

figure(2)
bar([x(:,1); x(:,2); x(:,3); x(:,4); x(:,5); x(:,6); x(:,7); x(:,8); x(:,9);
    x(:,10); x(:,11); x(:,12); x(:,13)])