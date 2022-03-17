close all; clear all; clc

% Known Parameters
h = 0.01;
tol = 10e-4;
t0 = 0.;
tf = 10.;

% Calculate number of iterations
N = int16((tf-t0)/h);

% Initialize
tival = zeros(N);
X = zeros(N,3);

% Initial Conditions s(0)=1, c(0)=0, p(0)=0
X(1,1) = 1.; 
x1 = X(1,1);
x2 = X(1,2);
x3 = X(1,3);
[dummy, neq] = size(X);

% Run RKF45 method
[X,hused,tival] = RKF45enzk(X,tival,tf,h,neq,tol);


% Extract substrate, complex and product solutions
s = X(:,1);
c = X(:,2);
p = X(:,3);

% Plot
figure(1)
plot(tival,c)
hold on
plot(tival,s)
hold on
plot(tival,p)
xlabel('Time')
ylabel('Concentration')
legend('Complex','Substrate','Product')
title('Enzyme Kinetics with s(0)=1, c(0) = 1, p(0) = 1')


% Plot iteration over time;

[dummy2, n] = size(tival); % get exact number of iterations

% Construct an array from 1 to n, incrementing by 1 (number of iterations)
iteration= linspace(1,n,n); 

% Plot
figure(2)
plot(iteration,hused(1:n,1))
xlabel('Iteration Number')
ylabel('Step Size')
title('Step Size vs. Iteration for RKF45 Method')
