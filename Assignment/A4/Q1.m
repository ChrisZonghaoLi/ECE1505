close all
clear all

%% 11.23(b) %%
load('hw4data.mat');
%W = W5;     % small case
W = W10;   % medium case
%W = W50;   % large case
[n,n] = size(W);
cvx_begin sdp
    variable X(n,n) symmetric
    minimize ( trace(W*X) )
    subject to
        diag(X) == 1;
        X >= 0;
cvx_end

% SDP lower bound d_prime
opt_value = trace(W*X)

% V is the eigenvector, diag(D) is the eigenvalue
[V,D] = eig(X);    
% return the max eigenvalue and its index
[val,idx] = max(diag(D));     
% return the eigenvector that corresponds to the max eigenvalue
v = V(:,idx);                   

x_hat = sign(v);
% heuristic value
heuristic_b = x_hat.'*W*x_hat
% compute the difference between heuristic and SDP lower bound d*
difference = heuristic_b - opt_value;
% small case: -2.2013e-07
% medium case: 8.8061
% large case: 613.5235

%% 11.23(c) %%
% generate randomized sample x with zero mean and X as the covariance
K = 100;
mu = zeros(n,1);
sigma = X;
rng('default')  % For reproducibility
x = mvnrnd(mu,sigma,K); % K rows, n column, so covariance is n*n
x_hat = sign(x);
for i = 1:1:K
    heuristic_c(i,:) = x_hat(i,:)*W*x_hat(i,:).';
end
min(heuristic_c)

%% 11.23(da) %%
% initialize vector x
x = ones(n,1);
iteration = 20;     % iteration 
heuristic_da = zeros(iteration,1);
for i = 1:1:iteration
    heuristic_da(i,:) = x.'*W*x;
    for j = 1:1:n
        x(j,:) = -x(j,:);
        buffer(j,:) = x.'*W*x;  % buffer of W to see which one returns the minimum objective value
        x(j,:) = -x(j,:);     % flip it back before the next run
    end
    [val,idx] = min(buffer);  % return the minimum value from buffer value and the corrresponding index
    x(idx,:) = -x(idx,:);
end 
min(heuristic_da)

%% 11.23(db) %%
% initialize vector x with uniform distribution {-1,1}
K = 100;
x = -2*round(rand(K,n))+ones(K,n);
iteration = 20;     % iteration 
heuristic_db = zeros(iteration,1);
for k = 1:1:K
    for i = 1:1:iteration
        heuristic_db(i,k) = x(k,:)*W*x(k,:).';
        for j = 1:1:n
            x(k,j) = -x(k,j);
            buffer(j,:) = x(k,:)*W*x(k,:).';  % buffer of W to see which one returns the minimum objective value
            x(k,j) = -x(k,j);     % flip it back before the next run
        end
        [val,idx] = min(buffer);  % return the minimum value from buffer value and the corrresponding index
        x(k,idx) = -x(k,idx);
    end 
end
[val,idx] = min(min(heuristic_db))
x_best = x(idx,:)

%% 11.23(dc) %%
% initialize vector x
K = 100;
mu = zeros(n,1);
sigma = X;
rng('default')  % For reproducibility
x = mvnrnd(mu,sigma,K); % K rows, n column, so covariance is n*n
x = sign(x);
iteration = 20;     % iteration 
heuristic_dc = zeros(iteration,1);
for k = 1:1:K
    for i = 1:1:iteration
        heuristic_dc(i,k) = x(k,:)*W*x(k,:).';
        for j = 1:1:n
            x(k,j) = -x(k,j);
            buffer(j,:) = x(k,:)*W*x(k,:).';  % buffer of W to see which one returns the minimum objective value
            x(k,j) = -x(k,j);     % flip it back before the next run
        end
        [val,idx] = min(buffer);  % return the minimum value from buffer value and the corrresponding index
        x(k,idx) = -x(k,idx);
    end 
end
[val,idx] = min(min(heuristic_dc))
x_best = x(idx,:)
    

    
    
    
    
    
    
    