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
% covF     The covariance matrix between the parameters
% F        The SCALED Fisher matrix
%
% Last modified by fjsimons-at-alum.mit.edu, 06/22/2015

% Default scaling is none
defval('scl',ones(size(th)))

% Scale up the parameter vector for the proper likelihood and score
th=th.*scl;

% First, the Fisher matrix at each wavenumber, unwrapped
mcF=Fisherkros(k,th,params);

% Take the expectation and put the elements in the right place
for ind=1:length(mcF)
  mcF{ind}=nanmean(mcF{ind});
end

% The upper half of the full Fisher matrix
% These will become the variances
F(1,1)=mcF{1};
F(2,2)=mcF{2};
F(3,3)=mcF{3};
F(4,4)=mcF{4};
F(5,5)=mcF{5};
F(6,6)=mcF{6};

% These will be the covariances of D with others
F(1,2)=mcF{7};
F(1,3)=mcF{8};
F(1,4)=mcF{9};
F(1,5)=mcF{10};
F(1,6)=mcF{11};

% These will be the covariances of f2 with others
F(2,3)=mcF{12};
F(2,4)=mcF{13};
F(2,5)=mcF{14};
F(2,6)=mcF{15};

% These will be the covariances of r with others
F(3,4)=mcF{16};  
F(3,5)=mcF{17};
F(3,6)=mcF{18};  

% These will be the covariances of s2 with others
F(4,5)=mcF{19};
F(4,6)=mcF{20};

% These will be the covariances of nu with others
F(5,6)=mcF{21};

% Returns the unscaled covariance matrix and the scaled Fisher matrix
disp('I am assuming that your wavenumbers are the entire plane')
[covF,F]=fish2cov(F,scl,length(k(~~k))*2);

% Output
varns={covF,F};
varargout=varns(1:nargout);
