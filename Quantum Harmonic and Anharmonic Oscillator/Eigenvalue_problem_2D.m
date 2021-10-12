%%%   part (a) 
%Create a variable t as syms t;
syms t;
%expm(X) = V*diag(exp(diag(D)))/V
A=[0 1 0 0; -1 0 0.1 0; 0 0 0 1; 0.1 0 -1 0];
At = t*A;
%[V,D]= eigs(At,4)
Bt =expm(At);
%we have matrix R^bar=exp(AT)*R^bar(0)
R0 = [1;0;0;0];
Rbar=Bt*R0;
%t=0:0.0627:62.7;
X1=Rbar(1,1);
X2=Rbar(3,1);
S1=simplify(X1)
S2=simplify(X2)
%Simplify does not work here
% so the solution can be written directly as 
t=0:0.1:62.7;
R1= cos((10^(1/2)*t*3)/10)/2 +cos((110^(1/2)*t)/10)/2;
R2= cos((10^(1/2)*t*3)/10)/2 -cos((110^(1/2)*t)/10)/2;
plot(t,R1,t,R2)

%% part (b)







