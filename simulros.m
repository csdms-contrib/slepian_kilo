function varargout=simulros(th0,params,xver)
% [Hx,Gx,th0,params,k,Hk,Gk,Sb,Lb]=SIMULROS(th0,params,xver)
%
% Simulates data under the CORRELATED two-layer Forsyth model as used
% by Olhede & Simons with the primary spectra S11 etc being the
% initial-loading ones. 
%
% S11 is the power spectral density of the INITIAL load at interface 1.
%
% INPUT:
%
% th0      The true parameter vector with elements:
%          th0(1)=D    Isotropic flexural rigidity 
%          th0(2)=f2   The sub-surface to surface initial loading ratio 
%          th0(3)=r    The sub-surface to surface initial correlation coefficient
%          th0(4)=s2   The first Matern parameter, aka sigma^2, in field units^2 
%          th0(5)=nu   The second Matern parameter 
%          th0(6)=rho  The third Matern parameter, in the units of the grid 
% params   A structure with constants that are (assumed to be) known:
%          DEL   surface and subsurface density contrast [kg/m^3]
%          g     gravitational acceleration [m/s^2]
%          z2    the positive depth to the second interface [m]
%          dydx  sampling interval in the y and x directions [m m]
%          NyNx  number of samples in the y and x directions
%          blurs  0 No wavenumber blurring
%                 1 No wavenumber blurring, effectively
%                 N Fejer convolutional  BLUROS  on an N-times refined grid
%                -1 Fejer multiplicative BLUROSY using exact procedure
%               Inf Simulate using SGP circulant embedding
%          kiso  [immaterial whether this is here or not]
%          quart 1 quadruple, then QUARTER [default]
%                0 size as is, watch for periodic correlation behavior
% xver     1 for extra verification, 0 if not needed
%
% OUTPUT:
%
% Hx       Real matrix of spatial-domain observations [m m], see Hk
% Gx       Real vector with Bouguer gravity anomaly data [m/s^2]
% th0      The true parameter vector pertaining to this data
% params   The structure with the known knowns, see above
% k        Wavenumber(s) suitable for the data sets returned [rad/m]
% Hk       A complex matrix of Fourier-domain observations, namely
%          final surface and subsurface topography [m]
% Gk       Complex vector with Fourier-domain Bouguer anomaly [m/s^2]
% Sb       The spectral matrix that you've used in this process
% Lb       The Cholesky decomposition of the spectral matrix which you
%          might use to evaluate the fit later on, in which case you
%          don't need any of the previous output
%
% EXAMPLE:
%
% simulros('demo1')
% simulros('demo2')
% simulros('demo3')
%
% SEE ALSO:
%
% MLEROS, LOADING, SIMULOS
%
% Last run on MATLAB Version: 9.7.0.1190202 (R2019b)
%
% Last modified by fjsimons-at-alum.mit.edu, 08/25/2026

% Check how it behaves when NOT a power of two! FFT should still be exact
% so wouldn't matter. Implement the windowing!

% Here is the true parameter vector and the only variable that gets used 
defval('th0',[7e22 0.4 -0.75 0.0025 2 2e4]);

% If not a demo...
if ~isstr(th0)
  % Supply the needed parameters, keep the givens, extract to variables
  fields={'DEL','g','z2','dydx','NyNx','blurs','kiso','quart'};
  defstruct('params',fields,...
	    {[2670 630],9.81,35000,[20 20]*1e3,[128 128],3,NaN,0});
  struct2var(params)

  % Here is the extra verification parameter
  defval('xver',1)

  % To create the gravity observation also
  G=fralmanac('GravCst');
  
  if xver==1
    % Dump to screen
    osdisp(th0,params)
  end

  % First make the wavenumbers, given the data size and the data length
  [k,dci,dcn,kx,ky]=knums(params);
  
  if xver
    % this should make sense as the spacing in wavenumber domain
    dkydkx=2*pi./NyNx./dydx;
    diferm(unique(diff(ky))-dkydkx(1))
    diferm(unique(diff(kx))-dkydkx(2))
  end

  % Bypass this procedure altogether to go for circulant embedding
  if ~isinf(params.blurs)
      % Now construct the whole-spectral matrix
      Z1=randgpn(k,dci,dcn);
      Z2=randgpn(k,dci,dcn);
      % This cramps the style But still. Should I put in a zero at dci? 
      % disp(sprintf('Z1: mean %+6.3f ; stdev %6.3f',...
      % 	       mean(Z1(:)),std(Z1(:))))
      % disp(sprintf('Z2: mean %+6.3f ; stdev %6.3f',...
      % 	       mean(Z2(:)),std(Z2(:))))

      % This handles all the cases 0, 1, >1 except Inf and -1
      if blurs>=0
          % FJS should make the gravity BEFORE blurring says SCO
          % S=S11.*T;
          % SG2=[2*pi*G*DEL(2)*exp(-k(:).*z2)]   .*S(:,2);
          % SG3=[2*pi*G*DEL(2)*exp(-k(:).*z2)].^2.*S(:,3);
          % Sb=bluros([S(:,1) SG2 SG3],params,xver);
          
          % Make the bivariate lithospheric coupling matrix
          [~,~,L,T]=Tros(knums(params,1),th0,params);
          % We need the (blurred) power spectrum - the theoretical quantity
          Sb=maternosp(th0,params,xver,[],T);
      elseif blurs<0
          % What if it is negative then deal with that
          % Need a way to convert initial-topography space based covariance
          % to final-interface space-based covariance i.e. through solving the
          % space-domain equation
          error('EXACT BLURRING Not ready yet')
      end

      % And then we do the Cholesky decomposition of the blurred lithospheric-spectral matrix, explicitly
      Lb=[sqrt(Sb(:,1)) Sb(:,2)./sqrt(Sb(:,1)) ...
	  sqrt((Sb(:,1).*Sb(:,3)-Sb(:,2).^2)./Sb(:,1))];
      
      % Should make sure that this is real! Why wouldn't it be?
      Lb=realize(Lb);

      if xver==1
          cholcheck(Lb,Sb,6,1)
          cholcheck(Lb,Sb,6,2)
          % If it turns out to not have been blurred we did have an analytical way
          if ismember(blurs,[0 1])
              diferm(Lb,sqrt(maternos(k,th0)).*L);
          end
      end

      % Blurred or unblurred, go on

      % And put it all together, unwrapped over k and over x
      Hk=[Lb(:,1).*Z1(:) [Lb(:,2).*Z1(:)+Lb(:,3).*Z2(:)]];
      
      % Without the L we should probably get the equilibrium topographies
      % So for Airy icebergs we should just have a simple scaling?
      if th0(1)==0 && th0(2)==0
          % This should be 1
          difer(Lb(:,1)-1,[],[],NaN)
          % This should have some relation to airyratio
          airyratio=DEL(1)/DEL(2);
          difer(Lb(:,2)+airyratio,[],[],NaN)
          % This should be 0
          difer(Lb(:,3),[],[],NaN)
      end

      % FJS need to change this to blur the chi rather than chi the blur 03/19/2014
      % With this sign convention the depth is positive
      Gk=2*pi*G*DEL(2)*exp(-k(:).*z2).*Hk(:,2);

      % And go to the space domain - unitary transform
      Hx(:,1)=tospace(Hk(:,1),params);
      Hx(:,2)=tospace(Hk(:,2),params);
      Gx     =tospace(Gk     ,params);

      if xver==1
          % Check Hermiticity before transformation, absolute tolerance
          hermcheck(reshape(Hk(:,1),NyNx))
          hermcheck(reshape(Hk(:,2),NyNx))
          hermcheck(reshape(k      ,NyNx))
          % Check unitarity of the transform; relative tolerance
          diferm(Hk(:,1)-tospec(Hx(:,1),params),[],9-round(log10(mean(abs(Hk(:,1))))));
          diferm(Hk(:,2)-tospec(Hx(:,2),params),[],9-round(log10(mean(abs(Hk(:,2))))));
          diferm(Gk     -tospec(Gx     ,params),[],9-round(log10(mean(abs(Gk     )))));
      end
  else
      % Make the Matern covariance OBJECT as required - vectorized
      % Calling it with Inf
      % Calling SGP then correlate then transform then times M then back to x
      keyboard
  end

  % Return the output if requested
  defval('Hk',[])
  defval('Sb',[])
  defval('Lb',[])
  varns={Hx,Gx,th0,params,k,Hk,Gk,Sb,Lb};
  varargout=varns(1:nargout);
elseif strcmp(th0,'demo1')
  svnm=th0;
  [Hx,Gx,th0,p,k,Hk,Gk]=feval(mfilename);
  struct2var(p)

  clf
  kelicol
  [ah,ha]=krijetem(subnum(2,2));
  [tl(1),cb(1),xc(1),xa(1)]=plotit(ah(1),Hx(:,1),size(k),...
		       'final surface','topography (%s)','km');
  [tl(2),cb(2),xc(2),xa(2)]=plotit(ah(2),Hx(:,2),size(k),sprintf(...
      'final subsurface z_2 = %i km',round(z2/1000)),'topography (%s)','km');
  [tl(3),cb(3),xc(3),xa(3)]=plotit(ah(3),Gx*1e5,size(k),sprintf(...
      'gravity anomaly with %s%s = %i kg m^{-3}',...
      '\Delta','\rho',DEL(2)),'Bouguer anomaly (%s)','mgal');

  % Maybe we can calculate the coherence here also?
  %   NW=3;
  %   [FX,FY,SX,SY,SXY,COH2]=...
  %     mtm3(reshape(Hx(:,1),NyNx),reshape(Gx,NyNx),NW,max(round(2*NW)-2,1));

  % Cosmetics
  they=linspace(1,NyNx(1),5);
  thex=linspace(1,NyNx(2),5);
  spunkm=(NyNx-1).*dydx/1000;
  set(ah,'YLim',they([1 end])+[-1 1]/2,...
	 'XLim',thex([1 end])+[-1 1]/2,...
	 'YTick',they,...
	 'XTick',thex,...
	 'YTickLabel',-spunkm(1)/2+(they-1)*dydx(1)/1000,...
	 'XTickLabel',-spunkm(2)/2+(thex-1)*dydx(2)/1000)
  longticks([ah cb])
  nolabels(ah(1:2),1)
  nolabels(ha(3:4),2)

  % Plot the parameters here
  axes(ah(4))
  axis([-1 1 -1 1])
  nolabels(ah(4)); noticks(ah(4)); box on; axis image
  
  xof=-0.75;
  tx(1)=text(xof, 0.75,sprintf('%s = %12.3g','D',       th0(1)));
  tx(2)=text(xof, 0.50,sprintf('%s = %12.3g','f^2',     th0(2)));
  tx(5)=text(xof, 0.25,sprintf('%s = %12.3g','r',       th0(3)));

  tx(3)=text(xof, 0.00,sprintf('%s = %12.3g','\sigma^2',th0(end-2)));
  tx(4)=text(xof,-0.25,sprintf('%s = %12.3g','\nu',     th0(end-1)));
  tx(5)=text(xof,-0.50,sprintf('%s = %12.3g','\rho',    th0(end)));
  moveh(ah(4),-getpos(cb(2),3))

  fig2print(gcf,'portrait')
  figdisp([],svnm,[],1,'pdf')
elseif strcmp(th0,'demo2')
  % Try to get close to the example we had in 2000
  NyNx=[100 100];
  D=2.4722e+23;
  f2=0.5;
  s2=0.001;
  nu=0.5;
  rho=30000;
  z2=15000;
  % Perform the simulation... not complete
  [Hx,Gx,th0,k,Hk,Gk,params]=simulros;
elseif strcmp(th0,'demo3')
  params.blurs=2;
  params.NyNx=[128 128];
  [Hx,Gx,th0,params,k,Hk,Gk,Sb,Lb]=feval(mfilename,[],params); 
  Lk=Lkros(k,th0,params,Hk); 
  imagesc(decibel(v2s(Lk)))
  axis image
end

% Plotting routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tls,cbs,xcb,xa]=plotit(aha,dats,nm,stronk,strink,unid)
axes(aha)
imagesc(reshape(dats,nm)); 
axis image
limc=halverange(dats,95,NaN);
set(aha,'clim',limc)
tls=title(stronk);
xa=xlabel(sprintf('mean %+6.3f ; stdev %6.3f %s',...
		  mean(dats(:)),std(dats(:)),unid));
cbs=colorbar('ver');
xcb=ylabel(sprintf(strink,unid));
set(cbs,'ylim',limc)
