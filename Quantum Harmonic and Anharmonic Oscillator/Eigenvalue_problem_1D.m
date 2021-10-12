%EIGENVALUE PROBLEMS IN MATLAB
clc;
clear all;
close all;
%Define matrix
M=[1 -0.1; -0.1 1];
% [V,D]= eigs(M,N) returns a (N × N) diagonal matrix D of 
%the N eigenvalues and a (N × N) matrix V whose columns are 
%the corresponding eigenvectors.
[V,D]= eigs(M,2);
%orthogonality condition
V'*V;
dot(V(:,1),V(:,2));
% Eigen vector will change as
A=V(:,1);
B=V(:,2);
%conjugate transpose of eigen vectors
A1=A';
B1=B';
% expanision coeff. (d)
% let arbitary vecttor be [0 1]
V1=[0 1];
c=V1*V;
