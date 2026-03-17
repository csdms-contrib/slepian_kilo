function olhedesimons1r(th,params)
% OLHEDESIMONS1R(th,params)
%
% Accurate cartoon of gravity/topography/flexure relations under the
% published scaling and with initial-load correlation
%
% th       The true parameter vector with elements:
%          D     Isotropic flexural rigidity [Nm]
%          f2    The sub-surface to surface initial loading ratio
%                the cases where f2=0 have been properly handled
%          r     The subsurface-to-surface initial loading correlation
%          s2    The first Matern parameter, aka sigma^2
%          nu    The second Matern parameter
%          rho   The third Matern parameter
% params   A structure with constants that are (assumed to be) known:
%          DEL   Surface and subsurface density contrast [kg/m^3]
%          g     Gravitational acceleration [m/s^2]
%          z2    The positive depth to the second interface [m]
%          dydx  Sampling interval in the y and x directions [m m]
%          NyNx  Number of samples in the y and x directions
%          blurs 0 Don't blur likelihood (N=1 has the same effect)
%                N Blur likelihood using the [default: N=2] resampled Fejer window
% 
% Last modified by fjsimons-at-alum.mit.edu, 04/03/2014

% Propagate the BLUROS change but need refresher on the DIFER statements.

% Generate random topography and gravity etc

% Make it work for f2=Inf
% Check out the mean values also - it drifts see pathological cases of Matern!
% Put in the free-air anomaly also
% Check isostasy re free-air anomaly

% Flexural rigidity (in Nm)
defval('D',7e22);

% Loading ratio (dimensionless)
defval('f2',0.4); 

% Initial-loading correlation (dimensionless)
defval('r',-0.75);

% The three Matern parameters
defval('s2',0.0125);
defval('nu',2);
defval('rho',2e4);

% The interface depth
defval('z2',8000);

% Put them into a vector together if you don't have input
defval('th',[D f2 r s2 nu rho]);

% Here is the extra verification parameter
defval('xver',0)

% Young's modulus 
defval('E',1.4e11);
% Poisson's ratio
defval('v',0.25);
% Conversion to Te [m]
Te=(th(1)*12*(1-v^2)/E)^(1/3);

disp(sprintf('T_e = %i km',round(Te/1000)))

% Extract the relevant lithospheric parameters again
D=th(1);

% Density contrast, gravitational acceleration, interface depth
fields={'DEL','g','z2','dydx','NyNx','blurs'};
% Note if you override one you have to override them all
defstruct('params',fields,...
	  {[2670 630],9.81,z2,[20 20]*1e3,[128 128],2});

% Extract the relevant loading parameters again
DEL=params.DEL;
g  =params.g;
NyNx=params.NyNx;
% Index location of the profile
hNy=params.NyNx(1)/2;

% Should be able to check this for Airy icebergs with D=0 and f2=0 
airyratio=params.DEL(2)/params.DEL(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the gravity and topography - optional blurring
[Hx,Gx,th0,params,k,Hk,Gk]=simulros(th,params);

blurs=params.blurs;
% Do the blurring here also
if blurs~=0
  % Blurs IS the refinement parameter; make new wavenumber grid
  dydx=params.dydx;
  k2=knum2(blurs*NyNx,[(blurs*NyNx(1)-1)*dydx(1) (blurs*NyNx(2)-1)*dydx(2)]);
  phi=phios(k2,D,DEL,g);
  xi=xios(k2,D,DEL,g);
  % Blur the product (not the product of blurs!)
  pxm=bluros(phi.*xi-1,params);
  % Blur the terms
  phi=bluros(phi,params);
  xi =bluros(xi ,params);
else
  phi=phios(k,D,DEL,g);
  xi = xios(k,D,DEL,g);
  pxm=(phi.*xi-1);
end

d1pd2xi=DEL(1)    +DEL(2)*xi;
d1phpd2=DEL(1)*phi+DEL(2);

% Extract the observed FINAL topography and gravity
H01k=Hk(:,1);
H02k=Hk(:,2);

% Our eqs (30)
dcon1=[phi.*d1pd2xi/DEL(2)          d1pd2xi/DEL(1)];
dcon2=[     d1phpd2/DEL(2)      xi.*d1phpd2/DEL(1)];

% Don't attempt the matrix product here, rearranging is very slow
% Find the INITIAL applied loads in the wavenumber domain
H1k=[dcon1(:,1).*H01k+dcon1(:,2).*H02k]./pxm;
H2k=[dcon2(:,1).*H01k+dcon2(:,2).*H02k]./pxm;

% At the zero wavenumber we're always in trouble - make it zero
% Maybe instead I should restore the mean?
% Discussion with Sofia Jan 31 2011 on the Lindeberg balance in the
% variance terms which implies that omitting one wavenumber shouldn't
% have an effect on the resulting estimate. Order statistics.
kzero=k==0;
H1k(kzero)=0;
H2k(kzero)=0;

% Check with APLILOD retroactively (which is never blurred)
if xver==1 && blurs~=1
  [a,b]=aplilod(H01k,H02k,k(:),D,DEL);
  difer(a-H1k,7,[],NaN)
  difer(b-H2k,7,[],NaN)
end

% Now extract the EQUILIBRIUM loads that produced the final ones
% subscripts are first loading then response
H11k= H1k*DEL(2).*xi ./d1pd2xi;
H12k=-H1k*DEL(1)     ./d1pd2xi;
H22k= H2k*DEL(1).*phi./d1phpd2;
H21k=-H2k*DEL(2)     ./d1phpd2;

% Not sure about this yet
H11k(kzero)=0;
H12k(kzero)=0;
H22k(kzero)=0;
H21k(kzero)=0;
H01k(kzero)=0;
H02k(kzero)=0;

% Check that you're adding up - after blurring central wavenumber will be
% messed up so this will need to be a special case also
difer(H1k- [H11k-H12k],9,[],NaN)
difer(H2k- [H22k-H21k],9,[],NaN)
difer(H01k-[H11k+H21k],9,[],NaN)
difer(H02k-[H12k+H22k],9,[],NaN)

% Coherence and admittance - note position of averaging!
[~,Ctops]=radavg(v2s([Hk(:,1).*conj(Gk)]));
[~,Cbot1]=radavg(v2s(Hk(:,1).*conj(Hk(:,1))));
[~,Cbot2]=radavg(v2s(Gk.*conj(Gk)));
coh=abs(Ctops).^2./Cbot1./Cbot2;
adm=Ctops./Cbot1;

% FINAL loads in the space domain
H01x=Hx(:,1);
H02x=Hx(:,2);

% INITIAL loads in the space domain
H1x=tospace(H1k,NyNx);
H2x=tospace(H2k,NyNx);

% Report on the correlation that we observe in this way
[R,P,Rm,Rp]=corrcoef(H1x(:),H2x(:));
disp(sprintf(...
    '\nObserved correlation coefficient %.3f < %.3f < %.3f with p-value %i%s\n',...
	     Rm(2),R(2),Rp(2),round(P(2)*100),'%'))

% EQUILIBRIUM loads in the space domain
H11x=tospace(H11k,NyNx);
H22x=tospace(H22k,NyNx);
H12x=tospace(H12k,NyNx);
H21x=tospace(H21k,NyNx);

% How does it stack up in the space domain? Nearly constant offsets?
difer(H1x- [H11x-H12x],9,[],NaN)
difer(H2x- [H22x-H21x],9,[],NaN)
difer(H01x-[H11x+H21x],9,[],NaN)
difer(H02x-[H12x+H22x],9,[],NaN)

% Start the plot
clf
[ah,ha,H]=krijetem(subnum(3,3));

% Go outside of the box
x=linspace(0,params.dydx(2)*NyNx(2)/1000,NyNx(2));
x(1)=x(1)-[x(2)-x(1)];
x(end)=x(end)+[x(end)-x(end-1)];
xl='distance (km)';
yl='elevation (km)';

% Initial topography at the surface
[p,pl]=plotit(x,H1x,0,hNy,ah(1),z2,xl,yl);
% Initial topography at the subsurface
[q,ql]=plotit(x,0,H2x,hNy,ah(7),z2,xl,yl);

% Equilibrium topography from the surface loading
[rt,rl]=plotit(x,H11x,H12x,hNy,ah(2),z2,xl,yl);
% Equilibrium topography from the subsurface loading
[s,sl]=plotit(x,H21x,H22x,hNy,ah(8),z2,xl,yl);

% Final topography at the surface and the subsurface
[t,tl,ylow]=plotit(x,H01x,H02x,hNy,ah(6),z2,xl,yl);

% Perform a last check of consistency
difer(H1x-[H11x-H12x],8,[],NaN)
difer(H2x-[H22x-H21x],8,[],NaN)

% Cosmetics, note the axis limits conversion 
%set(ah,'xlim',minmax(x),'ylim',[ylow 1.25*max(max([H1x H2x]))/1000])
set(ah,'xlim',minmax(x),'ylim',[ylow -ylow/2])
set(ah,'xtick',0:800:max(x),'ytick',[ylow ylow/2 0 -ylow/2])
longticks(ah)
nolabels(ah([2 8]),2)
delete([rl(2) sl(2)]); rl=rl(1); sl=sl(1);

movev(ah(1:2),-.1)
movev(ah(7:8),+.1)
moveh(ah(2),-.025)
moveh(ah(8),-.025)

% On the last plot, also the Bouguer gravity anomaly [mgal]
bgrav=rindeks(v2s(Gx),hNy)*1e5;
axi=ah(6);
axes(axi)
hold on
% Scale the anomaly to the layer, somewhat thicker though
info=2;
infl=info*max(abs(rindeks(v2s(H02x),hNy)))/1000;
yhi=max(get(axi,'ylim'));
gravpo=scale(bgrav,[-1 1]*infl);
% Maximum range and interval that you wish to consider 
gravtix=[-50:10:50];
gravtix=[gravtix(indeks(find(min(bgrav)<gravtix),1)) ...
	 gravtix(indeks(find(gravtix<max(bgrav)),'end'))];
% Symmetrize on smallest
gravtix=sort([[-1 1]*min(abs(gravtix)) 0]);
g=plot(x,yhi+infl+gravpo,'k');
gravti=  yhi+infl+interp1(bgrav,gravpo,gravtix);
% Add ticks for the labels
gravti=[-z2/1000 0 gravti];
gravtix=[0 0 gravtix];

longticks(axi)
set(axi,'ylim',[ylow [yhi+2*infl]*1.1])
%gl=ylabel();
axih=getpos(axi,4);
set(axi,'xtick',0:800:max(x))
shrink(axi,1,[yhi-ylow]/range(ylim))
[xax,xlx,ylx]=xtraxis(axi,[],[],[],gravti,gravtix,'mgal');
hh=get(xax,'ytickl');
hh(1,:)=32;
hh(2,:)=32;
set(xax,'ytickl',hh);

% Center y label
movev(ylx,gravti(end-1)-getpos(ylx,2))
%set(ah,'xgrid','on')
longticks(xax)
delete(ah([3 4 5 9]))

tx(1)=text(x(1),2*ylow+4,sprintf('%s = %g ; %s = %g ; %s = %g',...
				 'D',D,'f^2',f2,'r',r));
tx(2)=text(x(1),2*ylow-0,sprintf('%s = %g','\sigma^2',s2));
tx(3)=text(x(1),2*ylow-3,sprintf('%s = %g','\nu',nu));
tx(4)=text(x(1),2*ylow-6,sprintf('%s = %g','\rho',rho));

set([ah([1 2 7 8 6]) xax tx ylx pl ql rl sl tl xlx ylx],'FontS',8)

axes(ah(1))
h(1)=text('string','$\mathcal{H}_{1}$',...
	  'interpreter','latex','pos',[x(end)+range(x)/20 0]);
axes(ah(2))
h(2)=text('string','$\mathcal{H}_{11}$',...
	  'interpreter','latex','pos',[x(end)+range(x)/20 0]);
axes(ah(2))
h(3)=text('string','$\mathcal{H}_{12}$',...
	  'interpreter','latex','pos',[x(end)+range(x)/20 -z2/1000]);
axes(ah(7))
h(4)=text('string','$\mathcal{H}_{2}$',...
	  'interpreter','latex','pos',[x(end)+range(x)/20 -z2/1000]);
axes(ah(8))
h(5)=text('string','$\mathcal{H}_{21}$',...
	  'interpreter','latex','pos',[x(end)+range(x)/20 0]);
axes(ah(8))
h(6)=text('string','$\mathcal{H}_{22}$',...
	  'interpreter','latex','pos',[x(end)+range(x)/20 -z2/1000]);
axes(ah(6))
h(7)=text('string','$\mathcal{H}_{\circ 1}$',...
	  'interpreter','latex','pos',[x(end)+range(x)/20 0]);
axes(ah(6))
h(8)=text('string','$\mathcal{H}_{\circ 2}$',...
	  'interpreter','latex','pos',[x(end)+range(x)/20 -z2/1000]);
% Activate the last axis again!
axes(xax)

fig2print(gcf,'portrait')
figdisp
%figdisp([],[],[],1,'pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,l,ylow]=plotit(x,tops,bots,hin,ax,offs,xl,yl)
col1=grey(7);
col2=grey(5);

% Plots a cross section at index
axes(ax)
if length(tops)==1
  tops=zeros(1,size(v2s(bots),2));
else 
  tops=rindeks(v2s(tops),hin);
end

if length(bots)==1
  bots=zeros(1,length(tops));
else 
  bots=rindeks(v2s(bots),hin);
end

% Convert everything to km
tops=tops/1000;
bots=bots/1000;
% Fake the offset to a fraction of the true thickness
offs=offs/1000;

% Lower axis limit
ylow=-2*offs;

% Lowermost layer
verybots=zeros(1,length(bots))+ylow;
p(1)=fill([x fliplr(x)],[tops fliplr(bots)-offs],col1);
hold on
p(2)=fill([x fliplr(x)],[bots-offs verybots],col2);
hold off

% Labels 
l(1)=xlabel(xl);
l(2)=ylabel(yl);
