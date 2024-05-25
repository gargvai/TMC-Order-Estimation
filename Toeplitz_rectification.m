function Rtoep = Toeplitz_rectification(Rhat);

%Enforcing toeplitz structure in the matrix

%Code used for the paper
% V. Garg, P. Giménez-Febrer, A. Pagès-Zamora, and I. Santamaria,
%“Source enumeration via Toeplitz matrix completion,”
%in IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Barcelona, Spain, May 2020.

%Input: The given covariance matrix
%Output: Covariance matrix with Toeplitz structure

M = size(Rhat,1); %Number of elements in array

rm = zeros(1,M);
        for im = 0:M-1
            for ii = 1:M-im
            rm(im+1) = rm(im+1) + Rhat(ii,ii+im);
            end
            rm(im+1) = rm(im+1)/(M-im);
        end
        Rtoep = toeplitz(rm);
