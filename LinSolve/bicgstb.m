function [x, err,iter] = bicgstb(A, b, M1, M2, max_it, tol)
%  [x, err] = bicgstb(A, b, M1, M2, max_it, tol)
%
% bicgstab.m solves the linear system Ax=b using the 
% BiConjugate Gradient Stabilized Method with preconditioning.
%
% input   A        matrix
%         x        initial guess vector
%         b        right hand side vector
%         M 1, M2  pre and post preconditioner matrices
%         max_it   maximum number of iterations
%         tol      error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm

% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Adam Pidlisecky and Eldad Haber
% Last update, July 2006

  iter = 0;                                          

% initialization
 
  x = zeros(length(b),1); 
  bnrm2 = norm( b );
  alpha = 0;
  
  r = b;
  error = norm( r )/bnrm2; err(1)=error;

  omega  = 1.0;
  r_tld = r;
  ff = 0;

  for iter = 1:max_it,                              

     rho   = ( r_tld'*r );                          
     if ( rho == 0.0 ) break, end;

     if ( iter > 1 ),
        beta  = ( rho/rho_1 )*( alpha/omega );
        p = r + beta*( p - omega*v );
     else
        p = r;
     end;
 
     p_hat = M1\p; p_hat = M2\p_hat;
     v = A * p_hat;

     alpha = rho / ( r_tld'*v );
     s = r - alpha*v;
     if ( norm(s)/bnrm2 < tol ),                         
        x = x + alpha*p_hat;
        resid = norm( s ) / bnrm2;
        ff = 1;
        break;
     end

     s_hat = M1 \ s; s_hat = M2\s_hat;
     t = A*s_hat;

     omega = ( t'*s) / ( t'*t );
     x = x + alpha*p_hat + omega*s_hat;            
     r = s - omega*t;
     error = norm( r ) / bnrm2;
    % fprintf('it %d,  res = %e\n', iter, error);            
     err(iter+1) = error;
     if ( error <= tol ), ff = 1; break, end;
     if ( omega == 0.0 ), ff = 1; break, end;
     rho_1 = rho;

  end
 
if ff == 0, 
       fprintf(' bcg res = %e,  tol = %e\n',norm(A*x-b)/norm(b), tol);
end;
