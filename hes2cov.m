function covH=hes2cov(H,df)
% covH=HES2COV(H,df)
%
% Turns a Hessian matrix into a covariance matrix according to asymptotic
% maximum-likelihood theory, see Simons & Olhede (2013).
%
% INPUT:
% 
% H       A full-form Hessian matrix, e.g. from FMINUNC/FMINCON, or HESSIOSL
% df      The number of degrees of freedom needed for final scaling, in
%         our case, this will be the number of independent wavenumbers
%         that went into the construction of the Whittle likelihood
%
% OUTPUT:
%
% covH    The covariance matrix based on this Hessian
%
% SEE ALSO:
%
% FISH2COV, TRILOSI, COVTHOSL
%
% Last modified by fjsimons-at-alum.mit.edu, 10/13/2016

% Calculate this version of the covariance matrix
covH=inv(-H)/df;
