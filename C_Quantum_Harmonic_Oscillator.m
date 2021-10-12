clear all; 
close all;
hbar=1; 
m=1; 
%size of matrix
N=500; 
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
zo = sqrt(hbar/(m*omega))

%harmonic oscillator potential
V=(m/2*omega^2)*z.^2; 
figure(1);
plot(z*1e6,V)
xlabel('Distance z [micrometer]'); ylabel('Potential V [eV]')
title('Harmonic Oscillator Potential')

%Harmonic oscillator harmiltonian
H=cz*(diag(2*ones(N,1))+diag(-1*ones(N-1,1),1)+diag(-1*ones(N-1,1),-1))+diag(V);

%H(1,   1) = 0; H(  1, 2) = 0; H(2, 1) = 0; % So that f(0) = 0.

%H(N, N-1) = 0; H(N-1, N) = 0; H(N, N) = 0; % So that f(L) = 0.
%eigen value and eigenvector
[A,B]= eigs(H,N,'sm');
A=A./sqrt(dz);
%probability amplitude plot
figure(2)
plot(z*1e6, (A(:,1:4)))
xlabel('Distance [nmicrometer]'); ylabel('|\Psi|')
title('amplitude vs distance')

% Probability density plot
%figure(3)
%plot(z*1e6, conj(A(:,1:4)).*A(:,1:4))
%set(gca,'FontSize',15);
%axis tight
%xlabel('Distance [micrometer]'); ylabel('\Psi * \Psi')
%title('Probability density vs distance')

% Part (e)
j=1:N;
jmax=(zmax/(sqrt(2)*zo))^2;
%Energy eigen values
Ej= diag(B);
rj=((Ej/omega)-(j'-1/2))./(j'-1/2);
figure(4);
plot(j/jmax,rj);
set(gca,'FontSize',15);
axis tight
xlabel('j/jmax'); ylabel('rj')
