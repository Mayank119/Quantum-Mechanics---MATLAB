% part (a) of Anharmoni Oscillator
clc
clear all 
close all
%size of matrix
% Hamiltonian with eigen values
hbar=1; 
m=1; 
%size of matrix
N=200; 
% Maximum limit of potential zmax = 4*z_o
zmax=10e-6;
%Define axis of the potential in 1D
z=linspace(-zmax,zmax,N);
% Define Delta z
dz=2*zmax/N;
% Term other than potential in the hamiltonian
cz=hbar^2/(2*m*dz^2);
% wave length
omega=1e12;
% In the problem zo = hbar/(m*omega) = 1e-6 meter
zo = sqrt(hbar/(m*omega));
%Define Psi(0)
Psi0=exp(-((z-2*zo).^2)/(2*zo*zo));
%Finding normlization constant
SN= sum(conj(Psi0).*Psi0)*dz;
N0=1/sqrt(SN);
Psi0=N0*exp(-((z-2*zo).^2)/(2*zo*zo));
%Plot of Psi0
figure(1);
plot(Psi0);
%Harmonic Oscillator potential
V=(m/2*omega^2)*z.^2; 
%figure(1);
%plot(z*1e6,V)
%xlabel('Distance z [micrometer]'); ylabel('Potential V [eV]')
%title('Anharmonic Oscillator Potential')

%anharmonic Oscillator harmiltonian
H=cz*(diag(2*ones(N,1))+diag(-1*ones(N-1,1),1)+diag(-1*ones(N-1,1),-1))+diag(V);

%eigen value and eigenvector
[A,B]= eigs(H,N,'sm');
A=A./sqrt(dz);
%probability amplitude plot
figure(2)
plot(z*1e6, (A(:,1:4)))
xlabel('Distance [nmicrometer]'); ylabel('|\Psi|')
title('Probability amplitude vs distance')

% Coefficents of Eigen vectors for c_j = Sum over all z (dz * eigen vector * Psi0) 
c     = dz*A'*Psi0';
% time has 10000 points
nt    = 1000;
% maximum time
tmax  = 8*pi/omega;
% define a zero matrix
xar   = zeros(N,nt);
for it = 1:nt
    % Define time matrix, it goes from t=0 to t=1094
    t       = (it-1)*tmax/(nt-1);
    %define matrix tar for time
    tar(it) = t;
    Psit       = 0*Psi0;
    % vary diffent value of mu and calculate x(t)
    for mu = 1:N
        v   = A(1:N,mu);
        E  = B(mu,mu);
        % summ of all the eigen vectors evolving in time
        % time evolution of eigen vectors and summation of cofficent * time evolving eigen vactors 
        Psit   = Psit + c(mu)*exp(i*E*t/hbar)*v;
    end
    % xar is all the vectors from time [0,1094]
    xar(1:N,it) = Psit(1:N);
end
% |Psit|^2
Mp2= conj(xar).*xar;
% Plot of |Psit|^2
figure(3)
imagesc(Mp2)
colorbar
set(gca,'FontSize',15);
axis tight
xlabel('time');
ylabel('\mu');
title('|\Psi|^2');
% <z(t)/z0>
% Zt=<z(t)/z0> = Sum_mu [z(t)*|Psit|^2*dz]
Zt=z*Mp2*dz/zo;
%plot of Zt
figure(4)
plot(tar,Zt)
set(gca,'FontSize',15);
xlabel('t');
ylabel('\langlez(t)\rangle/z_0');
title('expectation value of Z(t)');










