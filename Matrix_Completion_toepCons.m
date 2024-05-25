function Rmc = Matrix_Completion_toepCons(Rhat,p,mu,alpha)

%Code for the paper
% V. Garg, P. Giménez-Febrer, A. Pagès-Zamora, and I. Santamaria,
%“Source enumeration via Toeplitz matrix completion,”
%in IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Barcelona, Spain, May 2020.

%Matrix Completion algorithm with Toeplitz constraint

%Inputs:
%Rhat: Covariance Matrix
%p: the rank of the matrix (low rank)
%mu and alpha: regularization parametes described in the paper


M = size(Rhat,1); %Number of elements in array
R0 = Rhat-diag(diag(Rhat)); % Zeros in Diagonal Terms
Rm = R0./R0;
%% Initializing factors
W = rand(M,p);

  C= 1;
  Ratio = 1;
  itr = 1;
  while C == 1
  for ii = 1:M
        idx = find(real(Rm(ii,:)) == 1);
        sum1 = 0; sum2 = 0;
        for mm = 0:M-1-ii
            sum1 = sum1+ (W(ii+1,:)*W(ii+mm+1,:)'*W(ii+mm,:));
            sum2 = sum2 + (W(ii+mm,:)'*W(ii+mm,:));
        end
        W(ii,:) = (R0(ii,:)*W + R0(:,ii)'*W + alpha*sum1) * inv(2*W(idx,:)'*W(idx,:) + mu*eye(p) + alpha*sum2);
  end
  itr = itr+1;
  Ratio_old = Ratio;
  Rmc = W*W';

%converge
 Ratio = norm((Rmc-diag(diag(Rmc)))-(Rhat-diag(diag(Rhat))),'fro')/norm((Rhat-diag(diag(Rhat))),'fro');
  if abs(Ratio_old-Ratio)/Ratio_old < 1e-5 || itr > 100
      C = 0;
  end
  end
