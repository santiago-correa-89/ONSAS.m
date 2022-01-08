% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

%md function that constructs the assembled Fext vector for one given BC

function fext = elem2NodalLoads ( Conec, indBC, elements, boundaryCond, Nodes )

  % declare output fext
  nnodes      = size( Nodes, 1)      ;
  fext        = zeros( 6*nnodes, 1 ) ;

  % get element indexes with current BC
  elemsWithBC = find( Conec(:,3) == indBC ) ;

  % extract BC load base vals and coord system
  loadCoordSys = boundaryCond.loadsCoordSys ;
  loadvals     = boundaryCond.loadsBaseVals ;

  loadedNodes = []                   ;

  % loop in elements with current BC
  for elemInd = 1:length( elemsWithBC )

    elem      = elemsWithBC( elemInd )       ;
    elemInd   = Conec( elem, 2 )             ;
    elemType  = elements( elemInd ).elemType ;

    %md nodal loads
    if strcmp( elemType, 'node') % node

      if strcmp( loadCoordSys, 'global' )
        nodes     = Conec( elem, 4+1 ) ;
      else
        error(' only global flag in load by now.');
      end
      elemNodeLoadsMatrix = loadvals ;

    %md truss
    elseif strcmp( elemType , 'truss') ; %
      error(' not yet.');

    %md frame
    elseif strcmp( elemType , 'frame') ; %
    
			if strcmp( loadCoordSys, 'global' )
        nodes     = Conec( elem, 4+(1:2) ) ;
      else
        error(' only global flag in load by now.');
      end
      
      % Local directions 
      % -----------------------------------------------------------------------------
      elemLength = sqrt ( sum( ( Nodes(nodes(2),:) - Nodes(nodes(1),:) ).^2, 2 ) ) ;
      
      dx = Nodes(nodes(2),1) - Nodes(nodes(1),1) ;
      dy = Nodes(nodes(2),2) - Nodes(nodes(1),2) ;
      dz = Nodes(nodes(2),3) - Nodes(nodes(1),3) ;
      
      dxdl = ( dx ) ./ elemLength ; 
			dydl = ( dy ) ./ elemLength ; 
			dzdl = ( dz ) ./ elemLength ;
      
      lxy = sqrt( dx^2 + dy^2 ) ;
      
      exL = [ dxdl dydl dzdl ]' ;
      
      if dy > 1e-5*elemLength 
				eyL = [ dx -dy 0 ]' / lxy ;
			else	 
				eyL = [ 0 1 0 ]' ; % Convention adopted
			end	
			
			ezL = cross(exL, eyL) ;	
			
			xG = [1 0 0]' ;
			yG = [0 1 0]' ;
			zG = [0 0 1]' ;
			 
			% Load vectors
			% -----------------------------------------------------------------------------	
      % qx global
      cosXX = exL.*xG ; cosYX = eyL.*xG ; cosZX = ezL.*xG ;
      qxG_L = qx*[ cosXX cosYX cosZX ] ;
      
      % qy global
      cosXY = exL.*yG ; cosYY = eyL.*yG ; cosZY = ezL.*yG ;
      qyG_L = qy*[ cosXY cosYY cosZY ] ;
    
			% qz global
      cosXZ = exL.*zG ; cosYZ = eyL.*zG ; cosZZ = ezL.*zG ;
      qzG_L = qz*[ cosXZ cosYZ cosZZ ] ;
    
			qxL = qxG_L(1)+qyG_L(1)+qzG_L(1) ; 
			qyL = qxG_L(2)+qyG_L(2)+qzG_L(2) ; 
			qzL = qxG_L(3)+qyG_L(3)+qzG_L(3) ; 
    
			nodal_qxL = qxL * elemLength * [ 1/2 0 0 0 0 0	;  
																			 1/2 0 0 0 0 0 ] ;
			nodal_qyL = qyL * elemLength * [ 0 0 1/2 0 0  elemLength/12 ; 
																			 0 0 1/2 0 0 -elemLength/12 ] ;
			nodal_qzL = qzL * elemLength * [ 0 0 0 -elemLength/12 1/2 0 ;
																			 0 0 0  elemLength/12 1/2 0 ] ;
			
    
			elemNodeLoadsMatrix = nodal_qxL + nodal_qyL + nodal_qzL ; 
			
    %md edge
    elseif strcmp( elemType , 'edge') ; %
      nodes          = Conec( elem, 4+(1:2) ) ;
      % vector from node 1 to node 2
      orientedVector = Nodes( nodes(2),:) - Nodes( nodes(1),:) ;

      lengthElem = norm( orientedVector ) ;
      thickness  = elements( elemInd ).elemTypeGeometry ;

      factor = lengthElem * thickness * 0.5 ;

      % check for plane state loads
      assert( sum( loadvals( [ 2 4 5 6 ] )==0 )==4,'error in loads of edge' )

      if strcmp( loadCoordSys, 'global' )
        Fx =   loadvals( 1 ) * factor ;
        Fy =   loadvals( 3 ) * factor ;
        Fz = 0 ;
      elseif strcmp( loadCoordSys, 'local' )
        % consider a 90 degrees rotation of the oriented vector of the line element
        Fx = - orientedVector( 2 ) / lengthElem * factor ;
        Fy =   orientedVector( 1 ) / lengthElem * factor ;
        Fz = 0 ;
      end % if global/local system

      elemNodeLoadsMatrix = ones( length(nodes), 1 )*[Fx 0 Fy 0 Fz 0] ;

      assert( size( elemNodeLoadsMatrix, 2)==6,"error, maybe missing thickness")



    %md triangle tension
    elseif strcmp( elemType , 'triangle') ; %

      nodes = Conec( elem, 4+(1:3) ) ;

      areaElem = 0.5 * norm( cross( ...
        Nodes( nodes(2),:) - Nodes( nodes(1),:) , ...
        Nodes( nodes(3),:) - Nodes( nodes(1),:) ...
        ) ) ;

      if strcmp( loadCoordSys, 'global' )

        Fx = loadvals( 1 ) * areaElem / 3 ;
        Fy = loadvals( 3 ) * areaElem / 3 ;
        Fz = loadvals( 5 ) * areaElem / 3 ;

        assert( sum( loadvals( [ 2 4 6 ] ) == 0 ) == 3, ...
          'error only pressure loads, not moments, create an issue!' );

      elseif strcmp( loadCoordSys, 'local' ) % local coordinates load

        dofsaux = nodes2dofs( nodes , 6 ) ;
        dofs    = dofsaux(1:2:length(dofsaux)) ;
        % compute the normal vector of the element
        normal  = cross( ...
          Nodes( nodes(2),:) - Nodes( nodes(1),:) , ...
          Nodes( nodes(3),:) - Nodes( nodes(1),: ) ) ;
        % and normalize it
        n = normal / norm( normal ) ;

        Fx = n(1) * loadvals( 5 ) * areaElem / 3 ;
        Fy = n(2) * loadvals( 5 ) * areaElem / 3 ;
        Fz = n(3) * loadvals( 5 ) * areaElem / 3 ;

        assert( sum( loadvals( [ 1 2 3 4 6 ] ) == 0 ) == 5, ...
          'error only normal pressure loads in local coords, create an issue!' );

      else
        loadCoordSys
        error(' loadsCoordSys field must be local or global.');
      end % if global/local system

      elemNodeLoadsMatrix = ones( length(nodes), 1 ) * [ Fx 0 Fy 0 Fz 0 ] ;

    end %if elemTypes
    %mdadd loads to matrix of loaded nodes
    loadedNodes = [ loadedNodes ; ...
                    nodes'  elemNodeLoadsMatrix ] ;
  end % for elements

  %md convert to assembled fext vetor
  if exist( 'loadedNodes' ) ~= 0
    for i=1:size( loadedNodes ,1)
      aux = nodes2dofs ( loadedNodes(i,1), 6 ) ;
      fext( aux ) = fext( aux ) + loadedNodes(i,2:7)' ;
    end
  end
