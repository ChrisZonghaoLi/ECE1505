close all
clear all

digitsOld = digits(7);

L = 13;
N = 8;
c = ones(L,1); % capacity
s = [1.2;0.6;0.6;zeros(8-3-1,1)];
A = [
     1 1 1 0 0 0 0 0 0 0 0 0 0;
    -1 0 0 1 0 1 0 0 0 0 0 0 0;
     0 0 -1 0 1 0 0 1 0 0 0 0 0;
     0 -1 0 -1 -1 0 1 0 0 0 0 0 0;
     0 0 0 0 0 0 -1 0 1 1 0 1 0;
     0 0 0 0 0 -1 0 0 -1 0 1 0 0;
     0 0 0 0 0 0 0 -1 0 -1 0 0 1;
     0 0 0 0 0 0 0 0 0 0 -1 -1 -1;
    ];

A_plus = A(1:end-1,:);

B = [
    0 0 0 1 0 1 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 1 0 1 0;
    ];
b = ones(3,1);

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 t
x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13]';    % the input are real anyway, so just to remove conj in front of each x
f0 = 0;
for i = 1:1:L
    f0 = f0 + x(i,:)./(1-x(i,:));
end

B_hat = B*x - b;

f1 = (-1/t) * log(-1^(size(B,1))*prod(B_hat));
f2 = (-1/t) * log(prod(x));
f3 = (-1/t) * log((-1)^L*prod(x-1));
f = f0 + f1 + f2 + f3;

m = size(b,1) + size(x,1) + size(c, 1);

H = hessian(f,[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13]);     % Hessian of f
J = jacobian(f,[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13])';   % Jacobian of f

alpha = 0.5*rand;      % alpha in the range of (0,0.5)
beta = rand;           % beta in the range of (0,1)

%% Infeasible start Newton method with log-barrier %%
% for setting the accuracy of log-barrier
t_barrier = 1; 
mu = 10;    % used for update t in each iteration

% intialize step size for Newton's method
t_newton = 1;

% initialize x
x = 0.3*ones(L,1);  % avoid 0.5 since there will be zero in the denominator

% initialize v
v = 0*ones(size(A_plus',2),1);

% define how many iterations one wants
iteration = 10;

% https://www.mathworks.com/help/symbolic/subs.html
% symbolic substitution for J and H
% initialize Jacobian and Hessian
Jacobian = subs(J,[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 t],[x' t_barrier]);
Hessian = subs(H,[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 t],[x' t_barrier]);

% repeat
t_int = zeros(iteration,1);     % t for log-barrier in each iteration
delay = zeros(iteration,1);     % delay in each iteration
residual = zeros(iteration,1);  % total residual in eaech iteration

for k=1:1:iteration
    
    % old dual residual
    r_dual1 = Jacobian + A_plus' * v;
    
    % old primal residual
    r_primal1 = A_plus * x - s;    
    
    % calculate del_x
    del_x = pinv(A_plus) * (s - A_plus*x);
    %del_x = (Hessian)\(-A_plus'*w-Jacobian);

    % calculate del_v
    w = pinv(A_plus') * (-Jacobian-Hessian*del_x);
    %w = -inv((A_plus/(Hessian)*A_plus'))*A_plus*inv(Hessian)*Jacobian;
    del_v = w - v;    
    
    % update Jacobian and Hessian
    Jacobian = subs(J,[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 t],[x'+t_newton*del_x' t_barrier]);
    Hessian = subs(H,[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 t],[x'+t_newton*del_x' t_barrier]);

    % new dual residual
    r_dual2 = vpa(Jacobian + A_plus' * (v + t_newton*del_v))
    
    % new primal residual
    r_primal2 = A_plus * (x + t_newton*del_x) - s
    
    while (norm([r_dual2, r_primal2], 2) > (1 - alpha * t_newton) * norm([r_dual1, r_primal1],2)) && (k>=2)
       t_newton = beta * t_newton;
    end
    
    % update x
    x = x + t_newton*del_x;
    
    % update v
    v = v + t_newton*del_v;
    
    % update t_barrier
    t_barrier = mu*t_barrier;
    t_int(k,:) = t_barrier;
    
    delay(k,:) = subs(f0,[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 t],[x' t_barrier]);
    
    % total residual
    residual(k,:) = norm([r_dual2, r_primal2], 2) + m/t_barrier;

end
    
    
