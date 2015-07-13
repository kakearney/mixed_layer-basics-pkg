function [y, dy, Splits, Diag, failflag] = ode4splits(odefun,tspan,y0,varargin)
%ODE4SPLITSNONNEG Like ode4, but returns intermediate components, no negative
%
% [y, dy, Splits, Diag, failflag] = ode4splits(odefun,tspan,y0,p1, p2, ...)
%
% This function extends the ode4 function to calculate diagnostic variables
% and additive components as the ODE is solved.  For the additive
% components aspect, I assume that the ODE function is of the form dy/dt =
% dy1 + dy2 + dy3 + ....  This is particularly designed for biological
% modules in mixed_layer, where the total change is a sum of processes
% (production, grazing, predation loss, etc).
%
% Input variables:
%
%   odefun:     function handle to ode, of form [db,Splitdb,Diag] =
%               fun(t,b,P) 
%
%   tspan:      vector of time steps to integrate over
%
%   y0:         vector or 2D array of initial conditions
%
%   p#:         additional parameters to pass to the ODE function
%
% Output variables:
%
%   y:          new values at each time step
%
%   dy:         dy/dt over each time step
%
%   Splits:     structure of dy additive components at each time step
%
%   Diag:       structure of diagnostic variables at each time step
%
%   failflag:   true if any component becomes negative within the
%               Runge-Kutta calculations.



if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

% try
f0 = feval(odefun,tspan(1),y0,varargin{:});
% catch
%   msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
%   error(msg);  
% end  

% y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  

szy = size(y0);
nt = length(tspan);

y = zeros([szy nt]);
dy = zeros([szy nt-1]); % defined for every interval
y(:,:,1) = y0;


% 
% neq = length(y0);
% N = length(tspan);
% Y = zeros(neq,N);
% F = zeros(neq,4);

% Y(:,1) = y0;


for i = 2:nt
    
  ti = tspan(i-1);
  hi = h(i-1);
  yi = y(:,:,i-1);

  [f1, S1, D1] = feval(odefun,ti,yi,varargin{:});
  [f2, S2, D2] = feval(odefun,ti+0.5*hi,yi+0.5*hi*f1,varargin{:});
  [f3, S3, D3] = feval(odefun,ti+0.5*hi,yi+0.5*hi*f2,varargin{:}); 
  [f4, S4, D4] = feval(odefun,tspan(i),yi+hi*f3,varargin{:});
 
  dy(:,:,i-1) = (hi/6)*(f1 + 2*f2 + 2*f3 + f4);
  y(:,:,i) = yi + dy(:,:,i-1);
  
  if i == 2
      fld = fieldnames(S1);
      Tmp = structfun(@(x) zeros(size(x)), S1, 'uni', 0);
      Splits(1:nt-1) = Tmp; 
      
  end
  
  for is = 1:length(fld)
      Splits(i-1).(fld{is}) = hi/6.*(S1.(fld{is}) + 2.* S2.(fld{is}) + 2.*S3.(fld{is}) + S4.(fld{is}));
  end
  
  Diag(i-1) = D1;
  
  nonnegcheck = yi >= 0 & ...
                yi+0.5*hi*f1 >= 0 & ...
                yi+0.5*hi*f2 >= 0 & ...
                yi+hi*f3 >= 0 & ...
                y(:,:,i) >= 0;
  if ~all(nonnegcheck(:))
      failflag = true;
      return;
  end
  
      
end

failflag = false;



