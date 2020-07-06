addpath('C:\fretAnalysisObject');

%%

myfret = fretAnalysisObject('C:\Users\Jozsef\Dropbox\Data for Jozsef\For FRET Histograms\170816_SNAPmGlu2-CyAC_100 uM L-Glu_LT20\#02Ch2\');

myfret.getInterferences(); 
myfret.parse_segResFinal();
myfret.calculateFret();

save( 'fretAnalysis#06ch2.mat', '-v7.3' )