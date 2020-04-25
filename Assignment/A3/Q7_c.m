close all
clear all

S = [9.5 8.5 7.5 6.5 5.5 4.5 3.5 2.5 1.5 0.5];
t = 10; % 10 time instances
y = zeros(t,1);

obj = [zeros(1,t),1];
A = [eye(t,t),-ones(t,1); -eye(t,t), -ones(t,1)];
b = [y;y];
Aeq = [ones(1,t),0; S, 0];
beq = [0; 1];

% formulate the LP problem in MATLAB
[opt_sol, fval] = linprog(obj,A,b,Aeq,beq );
optima = opt_sol(1:10);
x = optima;

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

