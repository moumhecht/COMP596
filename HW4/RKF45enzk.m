function [X,hused,tival] = RKF45enzk(X,tival,tf,h,neq,tol)
% Function used for implementation of Runge-Kutta-Fehlberg of 5th order
% Adapted from Mathews and Fink (2004). Numerical Methods using Matlab.
% Output:
% X = matrix of solutions
% hused= time step used at each subinterval

% Input:
% t0,tf = initial and final values for time interval
% X     = matrix of solutions
%         such that at tival(1) --> X(1,:) = Initial Conditions
% h     = Initial time step used for solution
% neq   = Number of equations for the Right Hand Side of the system (symbolic)
%         dX/dt = F(t,X). This number should be equal to the number of
%         input parameters of X in F(t,X)
% tol   = Tolerance used for adaptive step size

% Requirements:
% The user must write a MATLAB function called frhsrkf45 with the following
% definition:
% function  Fout = frhsrkf45(t,x1,x2,x3)  or similar (e.g. using varargin)
%     where: Fout: is a column vector of 3 components
%               t: time input
%              x1: first input variable 
%              x2: first input variable
%              x3: first input variable

% Parameters used for Runge-Kutta implementation
a2=1/4;b2=1/4;
a3=3/8;b3=3/32;c3=9/32;
a4=12/13;b4=1932/2197;c4=-7200/2197;d4=7296/2197;
a5=1;b5=439/216;c5=-8;d5=3680/513;e5=-845/4104;
a6=1/2;b6=-8/27;c6=2;d6=-3544/2565;e6=1859/4104;f6=-11/40;
r1=1/360;r3=-128/4275;r4=-2197/75240;r5=1/50;r6=2/55;

% Parameters used for 4th and 5th order Runge-Kutta implementation
n1=25/216;n3=1408/2565;n4=2197/4104;n5=-1/5;
% m1=16/135;m3=6656/12825;m4=28561/56430;m5=-9/50;m6=2/55;

% Limits for step adjustment
hmin = h/64;
hmax = 64*h;
blowup = 1e15;

% Initializing variables/vectors
k1 = zeros(neq,1);
k2 = zeros(neq,1);
k3 = zeros(neq,1);
k4 = zeros(neq,1);
k5 = zeros(neq,1);
k6 = zeros(neq,1);
hused = zeros(size(tival));
hused(1) = h;
ii = 2;
iint = 0; 

while ( (ii <= size(tival,2))&&(iint <= size(tival,2))&&(tival(ii-1) < tf) )
        
    iint = iint + 1;
    
    k1 = h*frhsrkf45( tival(ii-1) , X(ii-1,1) , X(ii-1,2) , X(ii-1,3) );
    if blowup < norm(abs(k1),2)
        fprintf('Parameters exceed limits \n');
        break
    end;
    
    k2 = h*frhsrkf45( tival(ii-1)+a2*h , ...
        X(ii-1,1)+k1(1)*b2 , X(ii-1,2)+k1(2)*b2 , X(ii-1,3)+k1(3)*b2);
    if blowup < norm(abs(k2),2) 
        fprintf('Parameters exceed limits \n');        
        break
    end;
    
    k3 = h*frhsrkf45( tival(ii-1)+a3*h , ...
        X(ii-1,1)+k1(1)*b3+k2(1)*c3 , ...
        X(ii-1,2)+k1(2)*b3+k2(2)*c3 , ...
        X(ii-1,3)+k1(3)*b3+k2(3)*c3);
    if blowup < norm(abs(k3),2) 
        fprintf('Parameters exceed limits \n');        
        break
    end;
    
    k4 = h*frhsrkf45( tival(ii-1)+a4*h , ...
        X(ii-1,1)+k1(1)*b4+k2(1)*c4+k3(1)*d4, ... 
        X(ii-1,2)+k1(2)*b4+k2(2)*c4+k3(2)*d4, ...
        X(ii-1,3)+k1(3)*b4+k2(3)*c4+k3(3)*d4);
    if blowup < norm(abs(k4),2) 
        fprintf('Parameters exceed limits \n');        
        break
    end;
    
    k5 = h*frhsrkf45( tival(ii-1)+a5*h , ...
        X(ii-1,1)+k1(1)*b5+k2(1)*c5+k3(1)*d5+k4(1)*e5, ... 
        X(ii-1,2)+k1(2)*b5+k2(2)*c5+k3(2)*d5+k4(2)*e5, ...
        X(ii-1,3)+k1(3)*b5+k2(3)*c5+k3(3)*d5+k4(3)*e5);
    if blowup < norm(abs(k5),2) 
        fprintf('Parameters exceed limits \n');        
        break
    end;
    
    k6 = h*frhsrkf45( tival(ii-1)+a6*h , ...
        X(ii-1,1)+k1(1)*b6+k2(1)*c6+k3(1)*d6+k4(1)*e6+k5(1)*f6, ... 
        X(ii-1,2)+k1(2)*b6+k2(2)*c6+k3(2)*d6+k4(2)*e6+k5(2)*f6);
    if blowup < norm(abs(k6),2) 
        fprintf('Parameters exceed limits \n');        
        break
    end;

    
    error = norm(r1*k1+r3*k3+r4*k4+r5*k5+r6*k6,Inf);
    
    xnew = X(ii-1,:)' + n1*k1 + n3*k3 + n4*k4 + n5*k5;
    
    if ((error<tol)||(h<2*hmin))
        
        X(ii,:) = xnew';
        hused(ii) = h;
        
        if ( (tival(ii-1)+h) > tf )
            tival(ii) = tf;
        else
            tival(ii) = tival(ii-1)+h;
        end
        ii =ii+1;
        
    end
    
    if (error == 0)
        s = 0;
    else
        s = 0.84*(tol*h/error)^(0.25);
    end
    
    if ((s < 0.75)&&(h > 2*hmin))
        h = h/2;
    end
    
    if ((s > 1.5)&&(2*h < hmax))
        h = 2*h;
    end
    
    if (blowup < norm(abs(X(ii-1,:)),2))       
        fprintf('Solution exceeds limits \n');
        break;
    end
    
end

if ( iint >= size(tival,2) )
    fprintf('Maximum number of iterations exceeded \n'); 
end
tival = tival(1:(ii-1));
X = X(1:(ii-1),:);
return
end