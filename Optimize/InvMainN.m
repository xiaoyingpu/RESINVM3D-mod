function [xc,d,MTX,para] = InvMainN(MTX,para,x0)
% [xc,d,MTX] = InvMainN(MTX,para,x0)
% 
%% Approximate Unconstrained Gauss-Newton with Armijo rule line search 
%
% Input: MTX = Structure of global (see generateMTX)
%        para = parameter structure
%        x0 = initial iterate
%
% Output: xc = solution
%         d = predicted data  
%
% Copyright (c) 2007 by the Society of Exploration Geophysicists.
% For more information, go to http://software.seg.org/2007/0001 .
% You must read and accept usage terms at:
% http://software.seg.org/disclaimer.txt before use.
% 
% Revision history:
% Original SEG version by Eldad Haber and Adam Pidlisecky, 
% Last update August 2006

%%%Screen output during inversion contains the following information:
%%iteration number; current objective function value; data misfit;
%%datanorm; modelnorm; relative gradient norm(current gradient/original gradient). 

%------ Initialize -----------------------------------------

disp('Initializing');

% Extract active cells from the ref model
mref = shiftdim(MTX.mref(find(MTX.ACTIVE)));



%  Scale the data
scl = mean(mean(abs(MTX.dobs)));
fprintf('Data scaling constant = %e\n\n',scl);

MTX.dobs = MTX.dobs/scl;
MTX.RHS = MTX.RHS/scl;
para.e = para.e/scl;
% Calculate the data weighting matrix
WdW = data_weight(MTX,para);
MTX.WdW = WdW;

% Extract data 
dobs = mkvc(MTX.dobs);


if exist('x0'),
    xc = x0(:);   
else,
    disp('Starting model = Ref model');
    xc = mref(:);     
end;

% allocate some numbers
misfit = []; gc = 1; normg0 = 1;

%set up the display 
disp('==================== START OPTIMIZATION ========================');
disp(' ');



itc = 0; 
% Start Gauss-Newton loop
while(norm(gc)/normg0 > para.tol & itc < para.maxit & norm(gc)>1e-20)
   % iteration count
   itc = itc+1;
   
   %%% ---- Function and gradient evaluation -------
   % Update model
   MTX = updateMTX(MTX,xc,para.ilutol);
   % Solve forward problem for the potnetial field everywhere in the volume
   if ~exist('MTX.U'); MTX.U  = Asol(MTX,MTX.RHS,para.inintol); end;

   % Calculate the (predicted) data from u
   d = Qu(MTX.OBS,MTX.U,MTX.SRCNUM);
 
   % Calculate the value of the objective function
   %%The data component
   fd = 0.5*(d-dobs)'*(WdW)*(d-dobs);
   %%the model component
    %automatically determine a beta guess
      if itc == 1 &  exist('x0') & isempty(para.BETA);
          
          para.BETA = 0.05*(fd./( (xc-mref)'*MTX.WTW*(xc-mref))); 
      elseif itc == 1  & isempty(para.BETA);
          %%Temporarily assign a value, to be corrected later
          para.BETA = 0;
      end;
   fm = 0.5*para.BETA*( (xc-mref)'*MTX.WTW*(xc-mref));

   %%Add them for the total Objective function
   fc =fd+fm;

   % Evaluate the gradient
   % model objective function gradient
      
   grad_fm = para.BETA*MTX.WTW*(xc-mref);

   % data objective function gradient (a little more complicated)
   QTd = Qtu(MTX.OBS, WdW*(d - dobs), MTX.SRCNUM );
   lm = ATsol(MTX, -QTd, para.inintol);
   %% use lm to obtain the gradient of the data component
   grad_fd = calcGvT(xc,MTX.U,MTX,lm);

   % Combine the gradients
   gold = gc;
   gc = grad_fd + grad_fm;
   
   %%%% Store some quantities for later use
   misfitold = misfit;
   misfit = sqrt((d-dobs)'*WdW*(d-dobs))/sqrt(dobs'*WdW*dobs);
   if itc == 1, normg0 = norm(gc); f0 = fc; mis0 = misfit; end;

    
   
astr = sprintf('%12s    %8s    %8s    %10s  %10s   %8s',...
               'iter', 'OBJ', 'Misfit', 'DataNorm', 'ModelNorm','RelGrad');

fprintf('%s\n',astr);
   %display parameters
   ostr = sprintf('%12s %12s %12s %12s %12s %12s',...
          [int2str(itc),'/',int2str(para.maxit)],...
          sprintf('%6.2f%c',100*fc/f0,'%'),...
          sprintf('%10.3e',misfit),...
          sprintf('%10.3e',fd),...
          sprintf('%10.3e',fm),...
          sprintf('%10.3e',norm(gc)/normg0));
  
          fprintf('%s\n',ostr);

   
   % mu_LS is mu for the line search

   if misfit < 1e-5
      disp(' misfit below tolerance');   
      return;
   end;

   % Approximate the direction
   % A conjugate gradient solver for
   % (J'*J + para.BETA*WTW) s = -gc
   MTX.mc = xc;
   s = ipcg(MTX, para.BETA, -gc, para.intol, para.ininintol, para.init); 
      
   % Test for convergence
   if max(abs(s)) < 1e-3,
      fprintf('    max_s = %e,  norm(g) = %e\n', max(abs(s)), norm(gc));
      fprintf('STEP size too small CONVERGE  ');  return; 
   end;

   % Try the step 
   mu_LS = 1; 
   iarm = 0;     
   % Line search
   while 1,
      xt = xc + mu_LS*s;
      %%%% Evaluate the new objective function
      MTX = updateMTX(MTX,xt,para.ilutol);
      MTX.U = Asol(MTX,MTX.RHS,para.inintol);
      
      d = Qu(MTX.OBS,MTX.U,MTX.SRCNUM);
   
      fd = 0.5*(d-dobs)'*WdW*(d-dobs);
      
      %automatically determine a beta guess
      if itc == 1 & para.BETA ==0;
          para.BETA = 0.5*(fd./( (xt-mref)'*MTX.WTW*(xt-mref))) 
      end;
      fm = 0.5*para.BETA*( (xt-mref)'*MTX.WTW*(xt-mref));       
      ft = fd+fm;
  
      fgoal = fc - para.alp*mu_LS*(s'*gc);

      ostr = sprintf(' %12s %12s %12s %12s %12s %12s',...
      [int2str(itc),'.',int2str(iarm+1)],...
      sprintf('%6.2f%c',100*ft/f0,'%'),...
      sprintf('%10.3e',misfit),...
      sprintf('%10.3e',fd),...
      sprintf('%10.3e',fm),...
      sprintf('%10.3e',[]));
    
      fprintf('%s\n',ostr);
 
      if ft < fgoal, 
        break,
      else
   	    iarm = iarm+1;
        mu_LS = mu_LS/2;    
      end;  
      
      if(iarm > 5)
           disp(' Line search FAIL EXIT(0)');     
           return;             
		  end
      fgoal = fc - para.alp*mu_LS*(s'*gc);
   end  % end line search
   
   % Update model
   xc = xt; 

   ss = ['iter.',int2str(itc),'.mat'];
   save(ss,'xc');
   
   misfitnew = misfit;
   misfitold = misfitnew;
   xc = xt;
end
%%Create the final model (insert the active cells into the right place)
xtemp = MTX.mref;
xtemp(find(MTX.ACTIVE)) = xc;
xc = xtemp;


