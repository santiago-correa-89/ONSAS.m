
close all, clear all %#ok

if isunix, dirSep = '/'; else dirSep = '\'; end
addpath( [ pwd  dirSep '..' dirSep  'src' dirSep ] ); octaveBoolean = isThisOctave ;

keyfiles = { 'staticVonMisesTruss/onsasExample_staticVonMisesTruss.m';
            'uniformCurvatureCantilever/onsasExample_uniformCurvatureCantilever.m' ;
            'linearPlaneStrainExample/onsasExample_linearPlaneStrain.m' ;
            'uniaxialExtension/onsasExample_uniaxialExtension.m'  } ;

current  = 1 ;   verifBoolean = 1 ;  testDir = pwd ;

while current <= length(keyfiles) && verifBoolean == 1

  % run current example
  fprintf([' === running script: ' keyfiles{current} '\n' ]);

  % save key files data to avoid clear all commands
  save( '-mat', 'exData.mat', 'current', 'keyfiles', 'dirSep', 'testDir' );

  run( [ pwd dirSep '..' dirSep 'examples' dirSep keyfiles{current} ] ) ;

  if verifBoolean
    status = 'PASSED';
  else
    status = 'FAILED';
  end

  % reload key files data and increment current
  load('exData.mat') ;

  fprintf([' === test problem %2i:  %s === \n\n'], current, status );

  current = current + 1 ;
  delete('exData.mat');
  cd ( testDir )
end

if verifBoolean ==1
  fprintf('test PASSED!\n')
else
  error('test examples not passed.')

end
