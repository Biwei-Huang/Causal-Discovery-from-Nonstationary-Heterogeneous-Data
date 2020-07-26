% random fourier feature to appropriate the kernel

function [ Z ] = transformFeatures(X)
%TRANSFORMFEATURES Transforms data to the random Fourier feature space
%
%   Input: 
%   X - n x p data matrix (each row is a sample) p is the dimension of the variable
%   Omega - p x D matrix of random Fourier directions (one for each
%   dimension of a sample x)
%   beta - 1 x D vector of random angles
%
%   Output:
%   Z - D x n matrix of random Fourier features

% sample random Fourier directions and angles
[T, p] = size(X);
if(T>=1000)
    D = 1000; % RFF dimension
else
    D = 500; 
end
Omega = randn(p,D); % RVs defining RFF transform
beta = rand(1,D)*2*pi; 

Z = cos(bsxfun(@plus,X*Omega,beta))*sqrt(2/D);
Z = Z';
