function varargout=dTros(k,th,params,phi,xi,pxm)
% [dToinvdD,dToinvdf2,dToinvdr]=dTros(k,th,params,phi,xi,pxm)
%
% Calculates first derivatives of isotropic Toinv under the CORRELATED
% initial-loading model as in Olhede & Simons (2013)

%
% k          Wavenumber(s) at which this is to be evaluated [1/m]
% th         The parameter vector with THREE elements (rest ignored)
%            th(1)=D    Isotropic flexural rigidity 
%            th(2)=f2   The sub-surface to surface initial loading ratio 
%            th(3)=r   The sub-surface to surface initial correlation coefficient
% params     A structure with AT LEAST these constants that are known:
%            DEL   surface and subsurface density contrast [kg/m^3]
%            g     gravitational acceleration [m/s^2]
% phi        Optionally precalculated phi, see PHIOS
% xi         Optionally precalculated xi, see XIOS
% pxm        Optionally precalculated (phi*xi-1)
%
% OUTPUT:
%
% dToinvdD    The first derivative of Toinv with respect to D
% dToinvdf2   The first derivative of Toinv with respect to f2
% dToinvdr    The first derivative of Toinv with respect to r
%
% Last modified by fjsimons-at-alum.mit.edu, 03/19/2014

% Extract the parameters from the input
D=th(1);
f2=th(2);
r=th(3);
DEL=params.DEL;
g=params.g;

if nargin<6
  % First the auxiliary quantities
  phi=phios(k,D,DEL,g);
  xi = xios(k,D,DEL,g);
  % Note that this has a zero at zero wavenumber
  pxm=(phi.*xi-1);
end

defval('pxp',pxm+2);

% Forcefully set f2 to a positive number even if it means a throw back
f2=abs(f2);
% Precompute some factors
f=sqrt(f2);
ddxi=[DEL(1)+DEL(2).*xi];

% First dToinvdD, which is a symmetric matrix

% The one without r
warning off MATLAB:divideByZero
% Eq. (A14)
fax=-2/f2/D./(xi-1).^2;
warning on MATLAB:divideByZero
dTinvdD(:,1)=fax.*...
    [1+f2*dpos(DEL,2,-2)+f2*dpos(DEL,1,-1)*(xi-1)];
dTinvdD(:,2)=fax.*...
    [dpos(DEL,-1,1)+phi/2-1/2+f2*dpos(DEL,1,-1)+1/2*f2*(xi-1)];
dTinvdD(:,3)=fax.*...
    [f2+dpos(DEL,-2,2)+dpos(DEL,-1,1)*(phi-1)];

% You still MIGHT have exactly zero... as in MLEROS('demo7')
if r~=0
  % The one with r
  warning off MATLAB:divideByZero
  % Eq. (A16)
  fax=-2/r/f/D./(xi-1).^2;
  warning on MATLAB:divideByZero
  dDTinvdD(:,1)=fax.*...
      [2*dpos(DEL,1,-1)+xi-1];
  dDTinvdD(:,2)=fax.*...
      [1+phi/2+xi/2];
  dDTinvdD(:,3)=fax.*...
      [2*dpos(DEL,-1,1)+phi-1];
  
  % The total one
  % Eq. (A20)
  dToinvdD=(1-r^2)^(-1)*(dTinvdD-r^2*dDTinvdD);
else
  % Don't ruin this by doing 0*Inf which turns it all to NaN
  dToinvdD=dTinvdD;
end

if nargout>=2
  % Then dToinvdf2, which is also a symmetric matrix

  % The one without r as in Tros
  warning off MATLAB:divideByZero
  % Eq. (A15)
  fax=-dpos(DEL,-2,0)*ddxi.^2/f2^2./pxm.^2;
  warning on MATLAB:divideByZero
  dTinvdf2(:,1)=fax;
  dTinvdf2(:,2)=fax*dpos(DEL,-1,1).*xi;
  dTinvdf2(:,3)=fax*dpos(DEL,-2,2).*xi.^2;

  % You still MIGHT have exactly zero... as in MLEROS('demo7')
  if r~=0
    % The one with r as in Tros
    warning off MATLAB:divideByZero
    % Eq. (A9)
    fax=dpos(DEL,-1,-1)*ddxi.^2/r/f./pxm.^2;
    warning on MATLAB:divideByZero
    invDT=[fax.*(               2*phi) ...
	   fax.*(  dpos(DEL,-1,1)*pxp) ...
	   fax.*(2*dpos(DEL,-2,2)*xi )];
    % Eq. (A9)
    dDTinvdf2=-invDT/2/f2;
    % The total one
    % Eq. (A21)
    dToinvdf2=(1-r^2)^(-1)*(dTinvdf2-r^2*dDTinvdf2);
  else
    dToinvdf2=dTinvdf2;
  end
else
  % Don't ruin this by doing 0*Inf which turns it all to NaN
  dToinvdf2=NaN;
end

if nargout>=3
  
  % You still MIGHT have exactly zero... as in MLEROS('demo7')
  if r~=0
    % Then dToinvdr, which is also a symmetric matrix
    warning off MATLAB:divideByZero
    % Eq. (A7)
    fax=DEL(1)^(-2)*ddxi.^2/f2./pxm.^2;
    warning on MATLAB:divideByZero
    % First the part without the correlation
    invT=[fax.*(                   1+f2*dpos(DEL,2,-2)*phi.^2) ...
	  fax.*(dpos(DEL,-1,1)*xi   +f2*dpos(DEL,1,-1)*phi   ) ...
	  fax.*(dpos(DEL,-2,2)*xi.^2+f2)];
    % The total one
    % Eq. (A22)
    dToinvdr=2*r/(1-r^2)^2*[invT-(1+r^2)/2*invDT];
  else
    % Don't ruin this by doing 0*Inf which turns it all to NaN
    dToinvdr=0;
  end
else
  dToinvdr=NaN;
end

% Output
varns={dToinvdD,dToinvdf2,dToinvdr};
varargout=varns(1:nargout);

