%md# Spring-mass system example
%md
%md[![Octave script](https://img.shields.io/badge/script-url-blue)](https://github.com/ONSAS/ONSAS.m/blob/master/examples/springMass/springMass.m)
%md
%mdIn this example a simple spring-mass system is considered. The notation and the analytical solution are based on chapter 3 from
%mdthe book Dynamics of Structures by Ray W. Clough and Joseph Penzien, Third Edition, 2003.
%md
%md First, the path to the ONSAS folder is added and scalar parameters are set.
% add path
close all, clear all; addpath( genpath( [ pwd '/../../src'] ) );
% scalar parameters for spring-mass system
k        = 39.47 ;
c        = 0.1   ;
m        = 1     ;
p0       = 40    ; % amplitude of applied load
u0       = 0.0   ; % initial displacement
%mdparameters for the equivalent truss model
l   = 1   ;
A   = 0.1 ;
rho = m * 2 / ( A * l ) ;
E   = k * l /   A       ;
%md
omegaN       = sqrt( k / m );
omegaBar     = 4*omegaN ;
xi           = c / m  / ( 2 * omegaN ) ;
nodalDamping = c ;
freq         = omegaN / (2*pi)      ;
TN           = 2*pi / omegaN        ;
dtCrit       = TN / pi              ;
%md
%md## Analytic solution
%md$$
%md  u(t) =
%md     ( A_c \cos( \omega_D  t ) + B \sin( \omega_D t ) ) e^{ -\xi \omega_N t } +
%md    G_1  \cos( \bar{\omega} t ) + G_2 \sin( \bar{\omega} t )
%md$$
%md
if (c == 0) && (p0 == 0) % free undamped
  myAnalyticFunc = @(t)   (   u0 * cos( omegaN * t )  ) ;
  analyticCheckTolerance = 2e-1 ;
else
  beta   = omegaBar / omegaN ;  omegaD = omegaN * sqrt( 1-xi^2 ) ;
  G1 = (p0/k) * ( -2 * xi * beta   ) / ( ( 1 - beta^2 )^2 + ( 2 * xi * beta )^2 ) ;
  G2 = (p0/k) * (  1      - beta^2 ) / ( ( 1 - beta^2 )^2 + ( 2 * xi * beta )^2 ) ;
  if u0 < l
    Ac = u0 - G1 ;
    B  =  (xi*omegaN*Ac - omegaBar*G2 ) / (omegaD);
  else
    error('this analytical solution is not valid for this u0 and l0');
  end
  myAnalyticFunc = @(t) ...
     ( Ac * cos( omegaD   * t ) + B  * sin( omegaD   * t ) ) .* exp( -xi * omegaN * t ) ...
    + G1  * cos( omegaBar * t ) + G2 * sin( omegaBar * t ) ;
    analyticCheckTolerance = 5e-2 ;
end

%md## Numerical solution
%md
%md### Numerical case 1: truss element model with newmark method
%md
%md#### Materials
%md
materials(1).hyperElasModel  = '1DrotEngStrain' ;
materials(1).hyperElasParams = [ E 0 ] ;
materials(1).density         = rho ;
%md
%md### Elements
%md
elements(1).elemType = 'node' ;
elements(2).elemType = 'truss';
elements(2).elemCrossSecParams = {'circle', [sqrt(4*A/pi) ] } ;
elements(2).massMatType = 'lumped' ;
%md
%md### Boundary conditions
%md
boundaryConds(1).imposDispDofs =  [ 1 3 5 ] ;
boundaryConds(1).imposDispVals =  [ 0 0 0 ] ;
%
boundaryConds(2).imposDispDofs =  [ 3 5 ] ;
boundaryConds(2).imposDispVals =  [ 0 0 ] ;
boundaryConds(2).loadsCoordSys = 'global'                  ;
boundaryConds(2).loadsTimeFact = @(t) p0*sin( omegaBar*t )                    ;
boundaryConds(2).loadsBaseVals = [ 1 0 0 0 0 0 ] ;
%md
%md### Initial conditions
initialConds.nonHomogeneousInitialCondU0    = [ 2 1 u0    ] ;
%md
%md### mesh
%md
mesh.nodesCoords = [  0  0  0 ; ...
                      l  0  0 ] ;
mesh.conecCell = { } ;
mesh.conecCell{ 1, 1 } = [ 0 1 1 0   1   ] ;
mesh.conecCell{ 2, 1 } = [ 0 1 2 0   2   ] ;
mesh.conecCell{ 3, 1 } = [ 1 2 0 0   1 2   ] ;
%md
%md and the following parameters correspond to the iterative numerical analysis settings
analysisSettings.methodName    = 'newmark' ;
analysisSettings.deltaT        =   0.005  ;
analysisSettings.finalTime      =   1.2*TN   ;
analysisSettings.stopTolDeltau =   1e-8 ;
analysisSettings.stopTolForces =   1e-8 ;
analysisSettings.stopTolIts    =   10   ;
%md
otherParams.problemName = 'springMass' ;
global exportFirstMatrices
exportFirstMatrices = true      ;
%md
[matUsNewmark, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
times = 0:analysisSettings.deltaT:(analysisSettings.finalTime+analysisSettings.deltaT) ;
valsAnaly = myAnalyticFunc(times) ;
valsNewmark = matUsNewmark(6+1,:) ;
%md
%md
%md### Numerical case 2: nodal mass model with $\alpha$-HHT method
%md
exportFirstMatrices = true      ;
materials(1).density  = 0 ;
materials(2).nodalMass = [m m m] ;
mesh.conecCell{ 2, 1 } = [ 2 1 2 0   2   ] ;
%md
analysisSettings.methodName    = 'alphaHHT' ;
analysisSettings.alphaHHT      =   0   ;
%md
otherParams.plotsFormat = 'vtk' ;
%md
[matUsHHT, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;
valsHHT  = matUsHHT(6+1,:) ;
%md
verifBooleanNewmark =  ( ( norm( valsAnaly - valsNewmark ) / norm( valsAnaly ) ) <  analyticCheckTolerance ) ;
verifBooleanHHT     =  ( ( norm( valsAnaly - valsHHT     ) / norm( valsAnaly ) ) <  analyticCheckTolerance ) ;
verifBoolean = verifBooleanHHT && verifBooleanNewmark ;
%md
%md## Plots
figure
hold on, grid on
plot(times, valsAnaly,'b-x')
plot(times, valsNewmark,'r-o')
plot(times, valsHHT,'g-s')
labx = xlabel('time');   laby = ylabel('\lambda(t)') ;
legend( 'analytic', 'truss-Newmark','nodalMass-HHT', 'location','northoutside')
print('output/springMassCheck.png','-dpng')
if exist('../../docs/src/assets/')==7
  % printing plot also to docs directory
  disp('printing plot also to docs directory')
  print('../../docs/src/assets/springMassCheck.png','-dpng')
end
%md
%md```@raw html
%md<img src="../../assets/springMassCheck.png" alt="plot check" width="500"/>
%md```
%md
