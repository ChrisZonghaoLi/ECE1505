close all
clear all

% p = optimvar('p',10);
% y = optimvar('y',10);
% %objec = transpose(p)*p;
% objec = (p(1))+(p(2))+(p(3))+(p(4))+(p(5))+(p(6))+(p(7))+(p(8))+(p(9))+(p(10));
% prob = optimproblem('Objective',objec);
% prob.Constraints.cons1 = sum(y) == 1;
% prob.Constraints.cons2 = p*y == 0;
% 
% problem = prob2struct(prob);
% [x,fval] = quadprog(problem)

H = ones(10,10);
%%%%constraint for a, uncomment this%%%%%%
%Aeq = [1 1 1 1 1 1 1 1 1 1;9.5 8.5 7.5 6.5 5.5 4.5 3.5 2.5 1.5 0.5];
%beq = [0;1];
%%%%constraint for b, uncomment this%%%%%%
Aeq = [1 1 1 1 1 1 1 1 1 1;9.5 8.5 7.5 6.5 5.5 4.5 3.5 2.5 1.5 0.5;...
    4.5 3.5 2.5 1.5 0.5 0 0 0 0 0];         % constraint for b                               
beq = [0;1;0];                              % constraint for b

[x,fval,exitflag,output,lambda] = ...
   quadprog(H,[],[],[],Aeq,beq);
% x = f, force

v = zeros(10,1);        % velocity
for t = 1:1:10
    v(t,:) = sum(x(1:t));
end

a = [0;x];              % add a zero to the front, acceleration
s = zeros(10,1);        % displacement 
for t = 1:1:10
    s(t,:) = sum(a(1:t)) + 0.5 * x(t,:);
end

p = zeros(10,1);        % position
for t = 1:1:10
    p(t,:) = sum(s(1:t));
end

figure(1)
p1 = plot(x);
set(gca,'linewidth',2)
set(p1, 'linewidth',3)
xlabel('time')
ylabel('force')
legend('force')
grid on

figure(2)
p1 = plot(v);
set(gca,'linewidth',2)
set(p1, 'linewidth',3)
xlabel('time')
ylabel('velocity')
legend('velocity')
grid on

figure(3)
p1 = plot(p);
set(gca,'linewidth',2)
set(p1, 'linewidth',3)
xlabel('time')
ylabel('position')
legend('position')
grid on

