function mcF=Fisherkros(k,th,params,xver)
% mcF=Fisherkros(k,th,params,xver)
%
% Calculates the entries in the Fisher matrix for Olhede & Simons (2013)
% for the CORRELATED initial-loading model prior to wavenumber averaging. 
%
% INPUT:
%
% k        Wavenumber(s) at which this is to be evaluated [1/m]
% th       The six-parameter vector with elements:
%          th(1)=D    Isotropic flexural rigidity 
%          th(2)=f2   The sub-surface to surface initial loading ratio 
%          th(3)=r    The sub-surface to surface initial correlation coefficient
%          th(4)=s2   The first Matern parameter, aka sigma^2 
%          th(5)=nu   The second Matern parameter 
%          th(6)=rho  The third Matern parameter 
% params   A structure with AT LEAST these constants that are known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
%
% OUTPUT:
%
% mcF      The 21-column Fisher matrix, in this order:
%          [1]  FDD    [2]  Ff2f2 [3] Frr
%          [4]  Fs2s2  [5]  Fnunu [6] Frhorho
%          [7]  FDf2   [8]  FDr   [9]  FDs2  [10] FDnu  [11] FDrho
%          [12] Ff2r   [13] Ff2s2 [14] Ff2nu [15] Ff2rho
%          [16] Frs2   [17] Frnu  [18] Frrho
%          [19] Fs2nu  [20] Fs2rho
%          [21] Fnurho
%
% EXAMPLE:
% 
% [~,~,th0,p,k]=simulros([],[],[],1);
% mcF=Fisherkros(k,th0,p,1);
%
% Last modified by fjsimons-at-alum.mit.edu, 01/08/2013

defval('xver',1)

% Extract the parameters from the input
D=th(1);
f2=th(2);
r=th(3);
DEL=params.DEL;
g=params.g;

% The number of parameters to solve for
np=length(th);
% The number of unique entries in an np*np symmetric matrix
npp=np*(np+1)/2;

% First the auxiliary quantities
phi=phios(k,D,DEL,g);
xi = xios(k,D,DEL,g);
% Note that this has a zero at zero wavenumber
pxm=(phi.*xi-1);

% First compute the "means" parameters
m=mAros(k,th,params,phi,xi,pxm,xver);

% Forcefully set f2 to a positive number even if it means a throw back
f2=abs(f2);
% Precompute some factors
f=sqrt(f2);

% Then compute the various entries in the order of the paper
lk=length(k(:));
% Some of them depend on the wave vectors, some don't
mcF=cellnan([npp 1],...
	    [lk 1 1 1 repmat(lk,1,7) 1 1 lk lk 1 repmat(lk,1,5)],...
	    repmat(1,1,npp));

% FDD, eq. (A34)
warning off MATLAB:divideByZero
fax=2*k(:).^8/g^2*dpos(DEL,-2,-2)./pxm.^2/f2./(1-r^2);
warning on MATLAB:divideByZero
mcF{1}=fax.*(2*f*dpos(DEL,1,1)*[f-3*r^2*f-r*f2-r]+...
	     f2*dpos(DEL,2,0)*[2+f2-r^2+2*r*f]+...
	     dpos(DEL,0,2)*[1+2*f2-r^2*f2+2*r*f]);
% Ff2f2, eq. (A35)
mcF{2}=(2-r^2)/2/f2^2/(1-r^2);

% Frr eq. (A36)
mcF{3}=2*(1+r^2)/(1-r^2)^2;

% Fthsths, eq. (131)
for j=4:np
  mcF{j}=2*m{j}.^2;
end

% FDf2, eq. (A37) might have gone via new ddTros but this is simpler
warning off MATLAB:divideByZero
fax=k(:).^4*dpos(DEL,-1,-1)./pxm/(1-r^2)/g/f^3;
warning on MATLAB:divideByZero
mcF{np+1}=fax*(2*f*dpos(DEL,0,1)-r^2*f*[dpos(DEL,1,0)+dpos(DEL,0,1)]...
	       -r*f2*dpos(DEL,1,0)+r*dpos(DEL,0,1));

% FDr, eq. (A38)
warning off MATLAB:divideByZero
fax=2*k(:).^4*dpos(DEL,-1,-1)./pxm/(1-r^2)/g/f;
warning on MATLAB:divideByZero
mcF{np+2}=fax*(f2*dpos(DEL,1,0)+dpos(DEL,0,1)...
	       -r*f*[dpos(DEL,1,0)+dpos(DEL,0,1)]);

% Fthlths, eq. (134)
for j=1:3
  mcF{j+8}=2*m{1}.*m{j+3};
end

% Ff2r, eq. (A39)
mcF{12}=-r/f2/(1-r^2);

% All further combinations Fthlths and Fthsthsp, eqs (134)-(135)
jcombo=nchoosek(1:np,2);
for j=7:length(jcombo)
  mcF{np+j}=2*m{jcombo(j,1)}.*m{jcombo(j,2)};
end

% Verification step
if xver==1
  % The "linear" tracecheck has already happed in MAROS
  % This should be twice the negative sum of the eigenvalues of
  % the product of this...
  [M,N,O]=dTros(k,th,params,phi,xi,pxm);
  % and this...
  [invTo,detTo,Lo]=Tros(k,th,params,phi,xi,pxm);
  % ... which we'll be again checking at a random wavenumber
  tracecheck(Lo,{M N O},{m{1:3}},9,1)
  % ... and then these also for ueber-completeness
  for ind=4:np
    tracecheck(Lo,{-repmat(m{ind},[size(invTo,1)/size(m{ind},1)],3).*invTo},...
	       {m{ind}},9,1)
  end
  % And then we check the sum of the squares of the eigenvalues which
  % hasn't happened so far
  tracecheck(Lo,{M N O},{mcF{1:3}},9,2)
end
