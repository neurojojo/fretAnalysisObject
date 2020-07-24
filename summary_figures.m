for N = 1:6; myfret1{N} = fretAnalysisObject(sprintf('C:\\#wes_smFRET\\#0%iCh2\\',N),'Apo'); end;
for N = 1:6; myfret2{N} = fretAnalysisObject(sprintf('C:\\Users\\Jozsef\\Dropbox\\Data for Jozsef\\For FRET Histograms\\170816_SNAPmGlu2-CyAC_100 uM L-Glu_LT20\\#0%iCh2\\',N),'Glu_100'); end; 
for N = 1:6; myfret3{N} = fretAnalysisObject(sprintf('C:\\Users\\Jozsef\\Dropbox\\Data for Jozsef\\For FRET Histograms\\170901_SNAPmGlu2-CyAC_15 uM L-Glu\\#0%iCh2\\',N),'Glu_15'); end; 
for N = 1:6; myfret4{N} = fretAnalysisObject(sprintf('C:\\Users\\Jozsef\\Dropbox\\Data for Jozsef\\For FRET Histograms\\170816_WA_SNAPf-mGlu2_CyAC_10 uM LY341495-antagonist\\#0%iCh2\\',N),'Antagonist'); end; 

%%

for N = 1:6;
    myfret1{N}.getInterferences(); 
    myfret1{N}.removeInterference();
    myfret1{N}.calculateFret('free');
    myfret1{N}.calculateFret('immobile');
end
close all

for N = 2:6
    myfret2{N}.getInterferences(); 
    myfret2{N}.removeInterference();
    myfret2{N}.calculateFret('free');
    myfret2{N}.calculateFret('immobile');
end


close all
for N = 1:6;
    myfret3{N}.getInterferences(); 
    myfret3{N}.removeInterference();
    myfret3{N}.calculateFret('free');
    myfret3{N}.calculateFret('immobile');
end



close all
for N = 1:6;
    myfret4{N}.getInterferences(); 
    myfret4{N}.removeInterference();
    myfret4{N}.calculateFret('free',1);
    myfret4{N}.calculateFret('immobile',1);
end
%%
for N = 1:6;
    myfret{N}.getInterferences(); 
    myfret{N}.removeInterference();
    %myfret{N}.calculateFret('free');
end
%%

for N = 1:6
    myfret{N}.calculateFret('free',1);
    myfret{N}.calculateFret('immobile',1);
end

%% Figure

alldata = [cellfun( @(x) mean( x.fretTraces.Fret.Calculated_free ), myfret1, 'ErrorHandler', @(x,y) nan )';...
cellfun( @(x) mean( x.fretTraces.Fret.Calculated_free ), myfret2, 'ErrorHandler', @(x,y) nan )';...
cellfun( @(x) mean( x.fretTraces.Fret.Calculated_free ), myfret3, 'ErrorHandler', @(x,y) nan )';...
cellfun( @(x) mean( x.fretTraces.Fret.Calculated_free ), myfret4, 'ErrorHandler', @(x,y) nan )'];
lbls = [ 1*ones(6,1); 3*ones(6,1); 2*ones(6,1); 4*ones(6,1) ];

figure('color','w');
h=boxplot( alldata, lbls ); hold on;

scatter( lbls, alldata , 64, 'Jitter', 'on', 'JitterAmount', 0.1, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerFaceAlpha', 0.5 )


all_lines = findobj(ax,'Type','Line');arrayfun( @(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines )
all_lines = findobj(gca,'Type','Line');arrayfun( @(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines )
myboxes = findobj(gca,'Tag','Box')arrayfun( @(box) patch( box.XData, box.YData, 'm', 'FaceAlpha', 0.5), myboxes(1:5) )
myboxes = findobj(gca,'Tag','Box'); arrayfun( @(box) patch( box.XData, box.YData, 'm', 'FaceAlpha', 0.5), myboxes(1:5) )
%% FREE

thisfretset = myfret4;

collect_Fret_values = arrayfun( @(N) thisfretset{N}.fretTraces.Fret.Calculated_free', [1:6] , 'UniformOutput', false, 'ErrorHandler', @(x,y) nan);

bins = [0:0.025:1];

results = cell2mat( cellfun( @(cell_) (histc( cell_, bins )./numel(cell_))', collect_Fret_values, 'UniformOutput', false ) )
figure('color','w'); 


counts = histc(cell2mat(collect_Fret_values),bins);
b=bar( bins, mean(results,2) ); hold on;

title( sprintf('%s - All cells (free)',regexprep(thisfretset{N}.experimentName,'_',' ')), 'color', 'k', 'fontweight', 'normal', 'fontsize', 24); box off;

hold on;
mygm = fitgmdist( cell2mat(collect_Fret_values)',1 );
counts = pdf( mygm, bins' ) * mean(diff(bins));
line( bins', counts, 'color', 'k', 'linewidth', 3 )
errorbar( bins, mean(results,2), std(results,[],2)/sqrt(6), 'linewidth', 2 )
set(gca,'color',[1,1,1],'xcolor',[0,0,0],'ycolor',[0,0,0],'linewidth',4,'fontsize',18,'TickDir','out');
saveas(gcf, sprintf('%s_summary_free.png',thisfretset{N}.experimentName) )

%IMMOBILE

thisfretset = myfret4;

collect_Fret_values = arrayfun( @(N) thisfretset{N}.fretTraces.Fret.Calculated_immobile', [1:6] , 'UniformOutput', false, 'ErrorHandler', @(x,y) nan);

bins = [0:0.025:1];

results = cell2mat( cellfun( @(cell_) (histc( cell_, bins )./numel(cell_))', collect_Fret_values, 'UniformOutput', false ) )
figure('color','w'); 


counts = histc(cell2mat(collect_Fret_values),bins);
b=bar( bins, mean(results,2) ); hold on;

title( sprintf('%s - All cells (immobile)',regexprep(thisfretset{N}.experimentName,'_',' ')), 'color', 'k', 'fontweight', 'normal', 'fontsize', 24); box off;

hold on;
mygm = fitgmdist( cell2mat(collect_Fret_values)',1 );
counts = pdf( mygm, bins' ) * mean(diff(bins));
line( bins', counts, 'color', 'k', 'linewidth', 3 )
errorbar( bins, mean(results,2), std(results,[],2)/sqrt(6), 'linewidth', 2 )
set(gca,'color',[1,1,1],'xcolor',[0,0,0],'ycolor',[0,0,0],'linewidth',4,'fontsize',18,'TickDir','out');
saveas(gcf, sprintf('%s_summary_immobile.png',thisfretset{N}.experimentName) )



%%
experiment=1

ch1_bleach_first = arrayfun( @(x) (myfret{experiment}.fretTraces.Ch2.traceMetadata(x).endOfTrace-myfret{experiment}.fretTraces.Ch1.traceMetadata(x).endOfTrace)>10, [1:myfret{experiment}.Ntracks] );
ch1_post_bleach = arrayfun( @(N) nanmedian( myfret{experiment}.fretTraces.Ch1.int_clean(N,[myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace:myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace+100]) ), [1:myfret{experiment}.Ntracks], 'ErrorHandler', @(x,y) nan ) 
ch2_post_bleach = arrayfun( @(N) nanmedian( myfret{experiment}.fretTraces.Ch2.int_clean(N,[myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace:myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace+100]) ), [1:myfret{experiment}.Ntracks], 'ErrorHandler', @(x,y) nan ) 

mytbl = table( [1:numel(ch1_post_bleach)]', ch1_bleach_first', ch1_post_bleach', ch2_post_bleach', sign((ch2_post_bleach'-ch1_post_bleach')), ((ch2_post_bleach'-ch1_post_bleach')), ...
    'VariableNames', {'Idx','Ch1_first','Ch1_pb','Ch2_pb','Ch2_over_Ch1','Ch2_sub_Ch1'})
mytbl = sortrows(mytbl,'Ch2_over_Ch1'); 
mytbl_sub = mytbl( mytbl.Ch1_first==true, : );
mytbl_sub = mytbl_sub( and( mytbl_sub.Ch2_pb > 0, mytbl_sub.Ch2_over_Ch1 > 0) , : )
mytbl_sub = sortrows(mytbl_sub,'Ch2_sub_Ch1')


%%
N=233
figure; plot( myfret{experiment}.fretTraces.Ch1.int_clean(N,:) ); hold on; plot( myfret{experiment}.fretTraces.Ch2.int_clean(N,:) );

xlim([0,myfret{experiment}.fretTraces.Ch2.traceMetadata(N).endOfTrace+200])
title( sprintf('%i %i', myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace, myfret{experiment}.fretTraces.Ch2.traceMetadata(N).endOfTrace) )
