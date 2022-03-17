function Fout = frhsrkf45(t,x1,x2,x3)

epsilon = 0.6;
lambda = 0.5;
k3 = 1;

Fout(1,1) = -epsilon*x1 + epsilon*(x1 + k3 - lambda)*x2;
Fout(2,1) =  x1 - (x1 + k3)*x2;
Fout(3,1) =  lambda*x2;

end
