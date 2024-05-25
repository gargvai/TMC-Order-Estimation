%Code for the paper
% V. Garg, P. Giménez-Febrer, A. Pagès-Zamora, and I. Santamaria,
%“Source enumeration via Toeplitz matrix completion,”
%in IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Barcelona, Spain, May 2020.

%Algorithm: Order Estimation using TMC

%Inputes:
  %Rhat = The covariance Matrix (after applying Toeplitz Rectification)
  % pmax, mu and alpha = As explained in the paper
  % Unlike explained in the paper, using a while loop here to make the algorithm cost effective.
  % The loop might stop before reaching the value pmax is some cases to save time

%Output: Khat, the required rank/order

%% change p and check norm of difference present and previous generated matrices

clear d
M = size(Rhat,1); %Number of elements in array
Rmc = zeros(M,M); % initializing as p=0 gives zero matrix
dnorm(1) = 0;
count = 1; % initializing
Cond = 1; %condition for while loop
while Cond == 1
    Rold = Rmc;
    Rmc = Matrix_Completion_toepCons(Rhat,count,mu,alpha);
    d(count) = norm(Rold - Rmc,'fro');
    dnorm(count+1) = norm(Rmc,'fro');
    if d(count)/d(1)<1e-20 || count == pmax
        Cond = 0;
    end
    count = count+1;
end
[~ Khat] = max(d(1:end-1)./d(2:end));


