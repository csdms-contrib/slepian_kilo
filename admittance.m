function varargout=admittance(Te,f,zm)
% [Q,k]=admittance(Te,f,zm)
% admittance(Te,f,zm)
%
% INPUT:
%
% Te    Elastic thickness      >0
% f     Loading fraction
% zm    Depth to the interface >0
%
% A serious function to calculate the admittance of two stochastic
% loading processes at the surface and depth zm, in a lithosphere with
% elastic thickness Te, and a loading fraction that doesn't depend on the
% orientation.
%
% Last modified by fjsimons-at-alum.mit.edu, 12/06/2006

% Define default values first, all in SI units
defval('Te',40*1e3);
defval('f',[0 0.5 1 2 5]);
defval('zm',35*1e3);
defval('k',logspace(-3,-0.8,200)/1000)
defval('D1',2670);
defval('D2',670);
defval('E',1e11);
defval('v',0.25);
defval('xver',1)

% Gravity 
g=fralmanac('GravAcc');
G=fralmanac('GravCst');

% Turcotte and Schubert (3-72)
D=(E*Te.^3)/(12*(1-v^2)); % Flexural Rigidity [Pa*m^3 or N*m]
disp(sprintf('D= %8.3e',D))

% Create grid on which to calculate admittance
[K,F]=meshgrid(k,f);

% Forsyth Eqs. (3) and (6)
xai=1+D.*K.^4/D2/g;
phi=1+D.*K.^4/D1/g;

% Bouguer admittance in accord with Forsyth's Eqs (11)-(12)
% If the sign of the exp is wrong, big mistake...
Q=-2*pi*G*D1*exp(-K*abs(zm)).*...
  (xai.^-1+phi.*F.^2*D1^2*D2^-2.*xai.^-2)./...
  (1+F.^2*D1^2*D2^-2.*xai.^-2);

% Convert to mgal/m
Q=Q/1e-5;

% OUTPUT or PLOT?
if nargout
  varnames={'Q','k'}
  for index=1:nargout
    varargout{index}=eval(varnames{index});
  end
else
  clf
  ah=gca;
  % Plot this up
  semilogx(k*1000,Q(:,:))
  set(ah,'ydir','rev')
  xl(1)=xlabel('wavenumber (rad/km)');
  yl(1)=ylabel('Bouguer-topography admittance (mgal/m)');
  axis tight
  yli=[0.02 -0.2];
  ylim(sort(yli))
  set(ah,'ytick',sort([yli(1) 0 -2*pi*G*D1/1e-5 yli(2)]))
  set(ah,'ygrid','on')
  set(ah,'xgrid','on')
  pos=[0.0102   -0.0402
       0.0121   -0.0527
       0.0141   -0.0652
       0.0172   -0.0857
       0.0252   -0.1087];
  for index=1:length(f)
    [bh(index),th(index)]=boxtex...
	([pos(index,1) pos(index,2)],ah,num2str(f(index)),10,[],0.8);
  end
  [bh(index+1),th(index+1)]=boxtex([0.003 0.0098],ah,'f',10,[],0.8);
  set(th,'horizontala','center','FontS',12)
  figdisp([],1)
  longticks(ah)
  ax=xtraxis1d(ah);
  xl(2)=xlabel('wavelength (km)');
  longticks(ax,2)
  set([ah gca xl yl],'FontS',12)
  set([xl yl],'FontS',15)
end



