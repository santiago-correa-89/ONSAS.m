% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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

% script for storing solver's input/output data structures.

modelCurrSol = struct ( 'timeIndex'          , timeIndex,...
                        'currTime'           , currTime,...
                        'U'                  , U , ...
                        'Udot'               , Udot , ...
                        'Udotdot'            , Udotdot , ...
                        'Stress'             , Stress , ...
                        'convDeltau'         , convDeltau, ...
                        'systemDeltauMatrix' , systemDeltauMatrix, ...
                        'nKeigneg'           , nKeigneg, ...
                        'nKeigpos'           , nKeigpos, ...
                        'factorCrit'         , factorCrit , ...
                        'timeStepStopCrit'   , timeStepStopCrit , ...
                        'timeStepIters'      , timeStepIters     ... 
                      ) ;

BCsData  = struct( 'constantFext'      , constantFext       , ...
                   'loadFactorsFunc'   , loadFactorsFunc    , ...
                   'variableFext'      , variableFext       , ...
                   'currLoadFactor'    , currLoadFactor     , ...
                   'nextLoadFactor'    , nextLoadFactor     , ...
                   'neumdofs'          , neumdofs           , ...
                   'KS'                , KS                 , ...
                   'userLoadsFilename' , userLoadsFilename    ...
                 ) ;

modelProperties = struct( 'coordsElemsMat'           , coordsElemsMat           ...
                        , 'Conec'                    , Conec                    ...
                        , 'materialsParamsMat'       , materialsParamsMat       ...
                        , 'elementsParamsMat'        , elementsParamsMat        ...
                        , 'crossSecsParamsMat'       , crossSecsParamsMat       ...
                        , 'nNodes'                   , nNodes                   ...
                        , 'nElems'                   , nElems                   ...
                        , 'booleanScreenOutput'      , booleanScreenOutput      ...
                        , 'booleanCSTangs'           , booleanCSTangs           ...
                        , 'booleanConsistentMassMat' , booleanConsistentMassMat ...
                        , 'nodalDispDamping'         , nodalDispDamping         ...
                        , 'numericalMethodParams'    , numericalMethodParams    ...
                        , 'stabilityAnalysisBoolean' , stabilityAnalysisBoolean ...
                        , 'problemName'              , problemName              ...
                        , 'outputDir'                , outputDir                ...
                        ) ;
