clear all
load('hw2data.mat','h')
load('hw2data.mat','y')

C = zeros(29,20);                                           % initialize circulant matrix
H = transpose([transpose(h) zeros(1,10)]);                  % impulse response with zero padding
B = zeros(size(H));                                         % buffer
m = 20;                                                     % dimension of g
n = 10;                                                     % dimension of h

%% Formulating Circulating Matrix
for i=-(n+(n-1)):1:-1
    B(-i+1:end) = H(1:end+i);
    C(i+m,:) = fliplr(transpose(B));
end

for q=m:1:(m+n-1)
    k = 1;                                                  % right shifting each row by one
    B = C(q-1,:);
    B(k+1:end)=B(1:end-k);
    C(q,:) = B;
end

%% Formulating the coefficient for Q.P., defined by matrix P
A = zeros(m,m,m+n-1);
for j=1:1:(m+n-1)                                           % in total there are 29 time samples are meaningful
    A(:,:,j) = transpose(C(j,:))*(C(j,:));                  % second-order term coefficicent for each time step
end

P = zeros(m:m);
for k=1:1:(m+n-1)  
    P = P + A(:,:,k);                                       % pull out the variable g, the summation of all second-order term coefficicent for each time step
end

%% Formulating Q.P., with variable D (delay), and variable g (equalizer)
H = P;
f = transpose(zeros(1,m)); 
objective = (zeros(1,m+n-1));
for D = 1:1:(m+n-1)
    Aeq = C(D,:);                                           %  the constraint at t = D
    beq = 1;
    [x(:,D),fval(:,D),exitflag(:,D),output,lambda] = ...
   quadprog(H,f,[],[],Aeq,beq);
	objective(:,D) = transpose(x(:,D))*P*(x(:,D)) - 1;      % objective function, {sum(cov(g*m))^2(t), t not equal to D}
end

%% Plot objective function versus D
D = 1:1:(m+n-1);
figure(1)
p1 = plot(D, objective);
set(gca,'linewidth',2)
set(p1, 'linewidth',3)
xlabel('delay D')
ylabel('sum of squared errors')
legend('delay vs sum of squared errors')
grid on

% D = 9 gives the optimal solution, check C(9,:)*x(:,9) to see if the
% constraint is satisfied

%% Apply deconvolution to the data
g = x(:,9);                                                 % best least-square equalizer
z = conv(y,g);
z = z(n-1:10000+(n-1)-1);                                   % 10000 samples, reconsrtucted

% normalize the data  
% for r=1:1:10000
%     z(r,:) = z(r,:)/abs(z(r,:));                                             
% end
figure(2)
hist(y)
xlabel('amplitude of y')
ylabel('count')

figure(3)
hist(z)
xlabel('amplitude of z')
ylabel('count')


        





