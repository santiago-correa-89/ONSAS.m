%md# A Linear Beam Element Analysis with Uniform Load Example
%md
%md## Previous definitions
close all, clear all;
addpath( genpath( [ pwd '/../../src' ] ) ) ; % add ONSAS directory to path
%md
%md scalar auxiliar parameters
E = 25e6 ; nu = 0.17 ; Lx = 1 ; Ly = 1 ; Lz = 1 ;
b = 0.45 ; % cross-section svgwidth
qx = 0  ; % applied uniform load in global X axis
qy = 1  ; % applied uniform load in global Y axis
qz = 0  ; % applied uniform load in global Z axis
%md
%md## MEBI parameters: Material-Element-BoundaryConditions-InitialConditions
%md
%md### Materials
materials.hyperElasModel = 'linearElastic' ;
materials.hyperElasParams = [ E, nu ] ;
%md### Elements
elements(1).elemType  = 'node'  ;
elements(2).elemType  = 'frame' ;
elements(2).elemTypeGeometry = [2, b, b] ;
elements(2).elemTypeParams   =  1 ;
%md### BoundaryConditions
% Supports
boundaryConds(1).imposDispDofs = [ 1 2 3 4 5 6 ] ;
boundaryConds(1).imposDispVals = [ 0 0 0 0 0 0 ] ;
% Loads
boundaryConds(2).loadsCoordSys = 'global' ;
boundaryConds(2).loadsBaseVals = [ qx 0 qy 0 qz 0] ;
%md### InitialConditions
%md empty struct
initialConds = struct() ;
%md
%md## Mesh
%md Mesh nodes
mesh.nodesCoords = ...
					[ 0 		0	 		0			; ...
						Lx/2  Ly/2  Lz/2 	; ...
						Lx   	Ly		Lz 		] ;
%md Conec Cell
mesh.conecCell = { } ;
%md Node elements. It is important to observe that only nodes with a boundary condition assigned are required to be included in the connectivity cell
mesh.conecCell{1,1} = [ 0 1 1 0 	1 ] ;
mesh.conecCell{2,1} = [ 0 1 1 0 	3 ] ;
%md Frame elements
mesh.conecCell{3,1} = [ 1 2 2 0   1 2 ] ;
mesh.conecCell{4,1} = [ 1 2 2 0   2 3 ] ;

% Analysis settings
analysisSettings = struct() ;

otherParams.problemName = 'linearBeamElement_uniformLoad' ;

[matUs, loadFactorsMat] = ONSAS( materials, elements, boundaryConds, initialConds, mesh, analysisSettings, otherParams ) ;

