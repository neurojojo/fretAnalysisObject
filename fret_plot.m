x = table2array(myfret.traceInterference{2});




%%
n=1
figure;
subplot(2,1,1);
plot( myfret.fretTraces.Ch1.int(n,:) )
ylim = get(gca,'YLim');
line( [myfret.fretTraces.Ch1.traceMetadata( n ).endOfTrace, myfret.fretTraces.Ch1.traceMetadata( n ).endOfTrace], ...
    ylim, 'color', 'k' )

subplot(2,1,2);
plot( myfret.fretTraces.Ch2.int(n,:) )
ylim = get(gca,'YLim');
line( [myfret.fretTraces.Ch2.traceMetadata( n ).endOfTrace, myfret.fretTraces.Ch2.traceMetadata( n ).endOfTrace], ...
    ylim, 'color', 'k' )

%%
N=1;
Y0 = myfret.fretTraces.Ch2.int(N,:);
Y = myfret.fretTraces.Ch2.int(N,:);
figure;
cutoff=0.5;
hold on;
crit = find( myfret.traceInterference{N}(:,1) < cutoff);
Y(crit) = nan;
%plot( Y0 ); 
hold on; plot(Y)