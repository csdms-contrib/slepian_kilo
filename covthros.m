function varargout=covthros(th,params,k,scl)
% [covF,F]=covFros(th,params,k,scl)
%
% Calculates the entries of the theoretical unblurred covariance matrix of
% the estimate, indeed the inverse Fisher matrix, hopefully close to the
% expectation of the Hessian of the actual simulations. As seen in Olhede &
% Simons (2013) for the CORRELATED Forsyth loading model.
%
% INPUT:
%
% th       The six-parameter vector (true or estimated) [scaled]:
%          th(1)=D    Isotropic flexural rigidity 
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=r    The sub-surface to surface initial correlation coefficient
%          th(4)=s2   The first Matern parameter, aka sigma^2
%          th(5)=nu   The second Matern parameter
%          th(6)=rho  The third Matern parameter
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
% k        Wavenumber(s) at which the Fisher matrix is evaluated [1/m]
% scl      The vector with any scalings applied to the parameter vector
%
% OUTPUT:
%
% covF     The full-form covariance matrix between the parameters
% F        The SCALED full-form Fisher matrix
%
% Last modified by fjsimons-at-alum.mit.edu, 06/26/2026

% Default scaling is none
defval('scl',ones(size(th)))

% Scale up the parameter vector for the proper likelihood and score
th=th.*scl;

% First, the Fisher matrix at each wavenumber, unwrapped
mcF=Fisherkros(k,th,params);

% Take the expectation and put the elements in the right place
for ind=1:length(mcF)
  cF(ind)=nanmean(mcF{ind});
end

% The full-form matrix
F=trilosi(cF);

% We really should be writing FISHIROS
% F=fishiros(...

% Returns the unscaled covariance matrix and the scaled Fisher matrix
disp('I am assuming that your wavenumbers are the entire plane')
covF=fish2cov(F,scl,length(k(~~k))*2);

% Output
varns={covF,F};
varargout=varns(1:nargout);
