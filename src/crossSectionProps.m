
function [ Area, J, Iyy, Izz, Jrho ] = crossSectionProps ( elemCrossSecParams, rho )

if strcmp( elemCrossSecParams{1}, 'generic' ) %general section

    elemCrossSecParamsVec = elemCrossSecParams{2} ;
    Area = elemCrossSecParamsVec( 1 ) ;
    J    = elemCrossSecParamsVec( 2 ) ;
    Iyy  = elemCrossSecParamsVec( 3 ) ;
    Izz  = elemCrossSecParamsVec( 4 ) ;

    if length( elemCrossSecParamsVec ) > 5
        Jrho =  diag( elemCrossSecParamsVec( 5:7 ) ) ;
    else
        Jrho = rho * diag( [ J Iyy Izz ] ) ;
    end

elseif strcmp( elemCrossSecParams{1}, 'rectangle' )
    elemCrossSecParamsVec = elemCrossSecParams{2} ;
    Area = elemCrossSecParamsVec( 1 ) * elemCrossSecParamsVec( 2 )          ;
    Iyy  = elemCrossSecParamsVec( 1 ) * elemCrossSecParamsVec( 2 )^3 / 12.0 ;
    Izz  = elemCrossSecParamsVec( 2 ) * elemCrossSecParamsVec( 1 )^3 / 12.0 ;

    % torsional constant from table 10.1 from Roark's Formulas for Stress and Strain 7th ed.
    a = .5 * max( elemCrossSecParamsVec(1:2) ) ;
    b = .5 * min( elemCrossSecParamsVec(1:2) ) ;

    J = a * b^3 * ( 16/3 - 3.36 * b/a * ( 1 - b^4 / ( 12*a^4 ) ) ) ;

    Jrho = rho * diag( [ J Iyy Izz ] ) ;

elseif strcmp( elemCrossSecParams{1}, 'circle' )
    diameter = elemCrossSecParams{2} ;
    Area = pi*diameter^2/4             ;
    Iyy  = pi*diameter^4/64            ;
    Izz  = Iyy                         ;
    J    = Iyy + Izz                   ;
    Jrho = rho * diag( [ J Iyy Izz ] ) ;
else
  error(' section type not implemented yet, please create an issue')
end
