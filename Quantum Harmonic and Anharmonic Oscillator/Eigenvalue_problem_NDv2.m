clc
clear ALL
close all
% define constants
m     = 20;
% size of matrix
N     = 2*m+1;
% Define x(0)
jar   = linspace(1,N,N);
% sigma can have two value 10 and 100
sigma = 100;
x0    = exp(-(jar-21).^2/sigma);
plot(jar,x0)
% define matrix
M     = diag(ones(N,1),0)-0.2*(diag(ones(2*m,1),1)+diag(ones(2*m,1),-1));
% Eigen value and Eigen Vectors where 'sm' arrange eigen values from low to high 
[V,D] = eigs(M,N,'sm');
% First 4 Eigen Vectors
v1    = V(1:N,1);
v2    = V(1:N,2);
v3    = V(1:N,3);
v4    = V(1:N,4);
% Coefficents of Eigen vectors for x(0)
c     = V'*x0';
% time has 10000 points
nt    = 10000;
% maximum time
tmax  = 1094;
% define a zero matrix
xar   = zeros(N,nt);

for it = 1:nt
    % Define time matrix, it goes from t=0 to t=1094
    t       = (it-1)*tmax/(nt-1);
    %define matrix tar for time
    tar(it) = t;
    x       = 0*x0;
    % vary diffent value of mu and calculate x(t)
    for mu = 1:N
        v   = V(1:N,mu);
        Om  = sqrt(D(mu,mu));
        % summ of all the eigen vectors evolving in time
        % time evolution of eigen vectors and summation of cofficent * time evolving eigen vactors 
        x   = x + c(mu)*cos(Om*t)*v;
    end
    % xar is all the vectors from time [0,1094]
    xar(1:N,it) = x(1:N);
end
% plot of irst four eigen vectors.
figure(1)
plot(jar,v1,jar,v2,jar,v3,jar,v4,jar,0*jar)
set(gca,'FontSize',15);
xlabel('\mu');
ylabel('x_\mu');
title('Eigenvectors');

%plot of x_\mu vs \mu for t=0,547 and 1094
figure(2)
subplot(3,1,1)
plot(jar,xar(:,1))
set(gca,'FontSize',15);
ylabel('x_\mu');
title('t = 0');
%
subplot(3,1,2)
plot(jar,xar(:,nt/2))
set(gca,'FontSize',15);
ylabel('x_\mu');
title('t = 547');
%
subplot(3,1,3)
plot(jar,xar(:,nt))
set(gca,'FontSize',15);
xlabel('\mu');
ylabel('x_\mu');
title('t = 1094');
% x_21 plot vs t for sigma = 100
figure(3)
plot(tar,xar(21,1:nt))
set(gca,'FontSize',15);
axis tight
xlabel('t');
ylabel('x_{21}(t)');
title('Center amplitude');
% c_j vs j
figure(4)
plot([1:20],abs(c(1:20)),'*','markersize',12)
set(gca,'FontSize',15);
axis tight;
xlabel('j');
ylabel('|c_j|');
title('Mode spectrum');
% For sigma =100 the gaussian function is broader than for sigma =10. 