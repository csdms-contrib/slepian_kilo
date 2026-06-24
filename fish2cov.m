function covF=fish2cov(F,scl,df)
% covF=FISH2COV(F,scl,df)
%
% Turns a Fisher matrix into a covariance matrix according to asymptotic
% maximum-likelihood theory, see Simons & Olhede (2013).
%
% INPUT:
% 
% F        A full-form Fisher matrix, e.g. from FISHIOSL
% df       The number of degrees of freedom needed for final scaling, in
%          our case, this will be the number of independent wavenumbers
%          that went into the construction of the Whittle likelihood
%
% OUTPUT:
%
% covF     The covariance matrix based on this Fisher matrix
%
% SEE ALSO:
%
% HES2COV, TRILOSI, COVTHOSL
%
% Last modified by fjsimons-at-alum.mit.edu, 10/31/2016

% Calculate this version of the covariance matrix
covF=inv(F)/df;
