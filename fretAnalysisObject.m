classdef fretAnalysisObject < handle
    
    
    % USAGE %
    %
    % Load in data:
    % myfret = fretAnalysisObject('C:\#wes_smFRET\#01Ch2\','Experiment name'); 
    % myfret.getInterferences(); 
    % myfret.removeInterference();
    % %In the next command, omit the 1 in order to avoid re-producing figures%
    % myfret.calculateFret('free',1); 
    % myfret.calculateFret('immobile',1);
    % myfret.make_cellViewTraces();
    
    properties
        images
        fretTraces
        ROImatrix
        dcmss
        experimentName
        
        filename
        cellnum
        
        Ntracks
        Ntimes
        
        DCMSSmatrix
        
        % Details about traces %
        Ch1
        Ch2
        traceInterference
        numPixels
        
        % Interference detection parameters %
        dw = 3;
        
        % Baseline fitting parameters %
        windowsize = 100; % VERY FLEXIBLE
        auROC_min = 0.3; % Needed to describe where there is a "stable" period
        NbaselinePts = 100; % These are the number of points needed to sample the baseline

    end
    
    methods
        
        function obj = fretAnalysisObject( experimentDir, experimentName )
            
            if exist( fullfile( experimentDir, 'fretTracesDif.mat' ) )
                    fprintf('\n----------------------------------------------\n');
                    fprintf('Loading traces from fretTracesDif.mat\n');
                    fprintf('----------------------------------------------\n');
                    fretTracesFile = fullfile( experimentDir, 'fretTracesDif.mat' );
            else
                    fprintf('\n----------------------------------------------\n');
                    fprintf('No appropriate fret Traces files were located\n');
                    fprintf('----------------------------------------------\n');
                   return 
            end
            
            obj.experimentName = experimentName;
            
            tmp = load( fretTracesFile, 'fretTraces' );
            obj.filename = fretTracesFile;
            obj.cellnum = regexp( obj.filename, '(?<=#)[0-9]{2}', 'match' );
            obj.fretTraces = tmp.fretTraces;
            
            img_files = dir(sprintf('%s\\ImageData\\*.tif',experimentDir));
            myimgs = arrayfun( @(i) imread( sprintf('%s\\ImageData\\%s',experimentDir,img_files(i).name) ), [1:size(img_files,1)], 'UniformOutput', false );
            obj.images = myimgs;
            
            obj.Ntracks = numel(obj.fretTraces.Ch2.traceMetadata);
            obj.Ntimes = numel(obj.images);
            
        end
        
        
        % A function used to get the bleaching times for the acceptor
        % channel %
        %
        % It will use dynamic referencing of the structure fields to
        % populate information for both ch1 and ch2 %
        function getBleach( obj, ch )
            
            for i = 1:size( obj.fretTraces.(ch).traceMetadata, 2 )
                
                obj.(ch).bleachInfo(i) = obj.getBleach_single( ch, i );
                
            end
            
        end
        
        function output = getBleach_single( obj, traceid, varargin )
            
            % CONSTANTS %
            windowsize = obj.windowsize; % VERY FLEXIBLE
            auROC_min = obj.auROC_min; % Needed to describe where there is a "stable" period
            NbaselinePts = obj.NbaselinePts; % These are the number of points needed to sample the baseline

            % BASES %
            % (1) obj.fretTraces.(ch).int(traceid,:) (raw data) is in full
            % frame basis [1:4000] and same for every trace
            % (2) z is in idxes basis [startOfTrace to endOfTrace+baseline]
            % (3) auroc is in idxes basis
            % (4) bleachpoint_frame is in frame basis [1:4000]
            % (5) bleachpoint_idxes is in idxes basis 
            
            [Ch1.start_, Ch1.end_, Ch1.baseline_] = deal( obj.fretTraces.Ch1.traceMetadata(traceid).startOfTrace, obj.fretTraces.Ch1.traceMetadata(traceid).endOfTrace, obj.fretTraces.Ch1.traceMetadata(traceid).lenBaseline );
            [Ch2.start_, Ch2.end_, Ch2.baseline_] = deal( obj.fretTraces.Ch2.traceMetadata(traceid).startOfTrace, obj.fretTraces.Ch2.traceMetadata(traceid).endOfTrace, obj.fretTraces.Ch2.traceMetadata(traceid).lenBaseline );
            
            % min(50,windowsize) is intended to catch situations where the 
            % DO NOT USE AUROC
            % aurocs = arrayfun( @(x) obj.auroc( z, x, min(50,windowsize) ), [1:numel(z)] );
            % [~,b] = min(aurocs);
            
            % Baseline calculation %
            % First, calculate auROCS for pre- and post-bleaching
            % [auroc_fret,auroc_baseline] = deal( aurocs(1:b), aurocs(b+1:end ) );
            
            % Second, split F data into pre- and post-bleaching
            z2 = obj.fretTraces.Ch2.int_clean(traceid,:); z2(z2==0)=nan; % Removing interferences creates 0
            z3 = obj.fretTraces.Ch2.int(traceid,:); % No clean intensity for post-bleaching is calculated
            Ch2.F_postbleach = nanmedian( z3(Ch2.end_: min( numel(z3), Ch2.end_+100 ) ) );
            
            z1 = obj.fretTraces.Ch1.int_clean(traceid,:); % Acceptor intensity
            Ch1.F_postbleach = nanmedian( z1(Ch2.end_: min( numel(z1), Ch2.end_+100 ) ) );
           
            % Next, stepwise adjustment of traces
            % Remove the offset after bleaching occurs
            z_1 = obj.fretTraces.Ch1.int_clean(traceid,:) - Ch1.F_postbleach;
            z_2 = obj.fretTraces.Ch2.int_clean(traceid,:) - Ch2.F_postbleach;
            
            Ch2.F_postAcceptor = nanmedian( z_2(Ch1.end_ : Ch2.end_) ); % Goes until donor (Ch2) bleaches
            Ch1.F_postAcceptor = nanmedian( z_1(Ch1.end_ : Ch2.end_) ); % Goes until donor (Ch2) bleaches
            
            Ch1.F_prebleach = nanmean( z1(Ch1.start_ : Ch1.end_) );
            Ch2.F_prebleach = nanmean( z2(Ch2.start_:Ch2.end_) );
            
            
            
            if ~isempty( varargin )
                figure('color','w'); 
                plot( z1 ); 
                line( [Ch1.start_, Ch1.end_], [Ch1.F_prebleach, Ch1.F_prebleach], 'color', 'r' )
                line( [Ch1.end_, Ch2.end_], [Ch1.F_postAcceptor, Ch1.F_postAcceptor], 'color','g','linewidth',2 );
                line( [Ch2.end_, numel(z1)], [Ch1.F_postbleach,Ch1.F_postbleach], 'color','k','linewidth',2 );
                set(gca,'TickDir','out'); box off;
                set(gcf,'Position',[35,558,1400,240]);
                xlim([Ch1.start_,numel(z1)]);
                
                figure('color','w'); 
                plot( z2 ); hold on; plot( Ch2.end_:numel(z3), z3(Ch2.end_:numel(z3)) );
                line( [Ch2.start_, Ch1.end_], [Ch2.F_prebleach, Ch2.F_prebleach], 'color', 'r' )
                line( [Ch1.end_, Ch2.end_], [Ch2.F_postAcceptor, Ch2.F_postAcceptor], 'color','g','linewidth',2 );
                line( [Ch2.end_, numel(z2)], [Ch2.F_postbleach,Ch2.F_postbleach], 'color','k','linewidth',2 );
                set(gca,'TickDir','out'); box off;
                set(gcf,'Position',[35,558,1400,240]);
                xlim([Ch2.start_,numel(z2)]);
                
            end
            
            output.Ch1 = Ch1;
            output.Ch2 = Ch2;
            
        end
        
        % This function is run only once both channels have had bleach
        % correction applied to them
        %
        % It will save outputs into obj.Ch1.int, obj.Ch2.int, and
        % obj.da_ratio
        function correctBleach( obj )
            
            bleachpoint_idxes_1   = arrayfun( @(x) x.bleachpoint_idxes , obj.Ch1.bleachInfo )';
            bleachpoint_idxes_2   = arrayfun( @(x) x.bleachpoint_idxes , obj.Ch2.bleachInfo )';
            
            prebleach_baseline_1   = arrayfun( @(x) x.prebleach_baseline , obj.Ch1.bleachInfo )';
            prebleach_baseline_2   = arrayfun( @(x) x.prebleach_baseline , obj.Ch2.bleachInfo )';
            
            postbleach_baseline_1   = arrayfun( @(x) x.postbleach_baseline , obj.Ch1.bleachInfo )';
            postbleach_baseline_2   = arrayfun( @(x) x.postbleach_baseline , obj.Ch2.bleachInfo )';
            
            Ncorrected = size( obj.Ch1.bleachInfo, 2 );
            
            da_ratio = arrayfun( @(traceId) (obj.Ch2.bleachInfo(traceId).prebleach_baseline-obj.Ch2.bleachInfo(traceId).postbleach_baseline)./(obj.Ch1.bleachInfo(traceId).prebleach_baseline-obj.Ch1.bleachInfo(traceId).postbleach_baseline), ...
                [1:Ncorrected] );
            corrected_Ch1 = arrayfun( @(x) da_ratio(x)*(obj.fretTraces.Ch1.int(x,:) - obj.Ch1.bleachInfo(x).postbleach_baseline), [1:Ncorrected], 'UniformOutput', false )';
            corrected_Ch2 = arrayfun( @(x) (obj.fretTraces.Ch2.int(x,:) - obj.Ch2.bleachInfo(x).postbleach_baseline), [1:Ncorrected], 'UniformOutput', false )';
            
            obj.Ch1.int = cell2mat( corrected_Ch1 );
            obj.Ch2.int = cell2mat( corrected_Ch2 );
            obj.da_ratio = da_ratio;
            
        end
        
        function printSave( obj )
            
            for i = 1:obj.Ntracks; obj.printSave_single(i,1); pause(0.05); close all; end
            
        end
        
        function printSave_single( obj, traceId, varargin )
           
            figure('color','w'); 
            ax_ = arrayfun(@(x) subplot(1,2,x,'NextPlot','add','TickDir','out'), [1:2]);
            plot( ax_(1), obj.fretTraces.Ch1.int(traceId,:) ); hold on; plot( ax_(1), obj.fretTraces.Ch2.int(traceId,:))
            plot( ax_(2), obj.Ch1.int(traceId,:) ); hold on; plot( ax_(2), obj.Ch2.int(traceId,:))
            set( gcf, 'position', [160,540,1400,250] );
            xlims = cell2mat(arrayfun(@(this_ax) this_ax.XLim, ax_, 'UniformOutput', false ));
            xlims = [ min(xlims), max(xlims) ]; 
            ylims = cell2mat(arrayfun(@(this_ax) this_ax.YLim, ax_, 'UniformOutput', false ));
            ylims = [ min(ylims), max(ylims) ]; 
            arrayfun( @(this_ax) set(this_ax,'XLim',xlims,'YLim',ylims), ax_ );
            
            mycluster = split( obj.filename, '\' ); experiments = regexprep( mycluster( numel(mycluster)-2 ), '_', ' ' );
            suptitle( sprintf('%s Trace %03i of %i',experiments{1},traceId,obj.Ntracks) );
            set(ax_(1),'Title',text(0,0,'Uncorrected','FontName','Arial','FontWeight','normal'))
            ax_(1); legend({'Ch1','Ch2'})
            set(ax_(2),'Title',text(0,0,'Corrected','FontName','Arial','FontWeight','normal'))
            ax_(2); legend({sprintf('Ch1 Prebleach=%1.2f Postbleach=%1.2f',obj.Ch1.bleachInfo(traceId).prebleach_quality, obj.Ch1.bleachInfo(traceId).postbleach_quality),...
                sprintf('Ch2 Prebleach=%1.2f Postbleach=%1.2f',obj.Ch2.bleachInfo(traceId).prebleach_quality, obj.Ch2.bleachInfo(traceId).postbleach_quality)});
            
            arrayfun( @(this_ax) set(this_ax,'XLabel',text(0,0,'Frame'),'YLabel',text(0,0,'F')), ax_ );
            
            savedir = regexp( obj.filename, '.*(?<=\\$)', 'match' );
            savedir = sprintf( '%s\\Savefiles\\', savedir{1} );
            
            if exist(savedir)==0; mkdir(savedir); end;
            if nargin>2
                if varargin{1}
                    saveas( gcf, sprintf( '%s%s_Trace_%i_of_%i.png', savedir, regexprep(experiments{1},' ','_'), traceId, obj.Ntracks ) );
                end
            end
            
        end
        
        function checkCorrection( obj, traceId )
            
            [ch1_baseline, ch2_baseline] = deal( obj.Ch1.baseline(traceId), obj.Ch2.baseline(traceId) );
            da_ratio = obj.Ch2.F0(traceId)./obj.Ch1.F0(traceId);
            
            myfig = figure('color','w');
            ax_ = arrayfun( @(x) subplot( 3,2,x,'NextPlot','add','Tag',sprintf('Ax_%i',x),'Parent',myfig ), [1:6] );
            
            
            plot(ax_(1), ch1_int); title('Ch1_intensity'); set(ax_(1),'Title',text(0,0,'Uncorrected Ch1'));
            plot(ax_(2), ch2_int); title('Ch2_intensity'); set(ax_(2),'Title',text(0,0,'Uncorrected Ch2'));
            
            plot(ax_(3), da_ratio * ( ch1_int - ch1_baseline ) );  set(ax_(3),'Title',text(0,0,'DA & BL corrected Ch1'));
            plot(ax_(4), ch2_int - ch2_baseline );  set(ax_(4),'Title',text(0,0,'DA & BL corrected Ch2'));
            
            plot(ax_(5), da_ratio * ( ch1_int - ch1_baseline ) );
            plot(ax_(5), ch2_int - ch2_baseline ); set(ax_(5),'Title',text(0,0,'Corrected Ch1 & Ch2'));
            
            set(ax_(6),'Title',text(0,0,'Empty'));
        end
        
        
        % A function used to calculate correlations  %
        function getInterferences( obj )
           
            [imrows,imcols] = deal(size(obj.images{1},1), size(obj.images{1},2) );
            template = padarray( padarray( fspecial('gaussian',3,1), 1, -5 )', 1, -5 );
            periphery = padarray( padarray( zeros(3,3), 1, 1 )', 1, 1 );
            
            obj.traceInterference = nan( obj.Ntracks, obj.Ntimes );
            obj.numPixels = nan( obj.Ntracks, obj.Ntimes );
                
            for N = 1:obj.Ntracks
                    
                dw=1;
                
                idx = [ obj.fretTraces.Ch2.traceMetadata(N).startOfTrace:obj.fretTraces.Ch2.traceMetadata(N).endOfTrace ]; 
                score_ = [];
                background = [];
                numpixels = [];
                nan_removed_x = fix( obj.removeNaNs( obj.fretTraces.Ch2.x(N,:) ) );
                nan_removed_y = fix( obj.removeNaNs( obj.fretTraces.Ch2.y(N,:) ) );
                
                for frame_ = idx
                    
                    x = nan_removed_x(frame_);
                    y = nan_removed_y(frame_);
                    
                    background(:,:,frame_) = obj.images{frame_}( [max(y-2*dw,1):min(y+2*dw,imcols)], [max(x-2*dw,1):min(x+2*dw,imrows)] );
                    
                end
                
                % Convolve the data using a 3d gaussian
                background2 = imgaussfilt3(background, [0.1,0.1,10] )+imgaussfilt3(background, [0.1,0.1,3] );
                % Compare each frame in the convolved movie to a template
                for frame_ = idx
                    score_(frame_) = corr2( background2(:,:,frame_), template );
                end
                
                periphery_frames = find( score_ > obj.min_corr );
                periphery_v_t = fliplr(arrayfun( @(x) mean2(background(:,:,x).*periphery), periphery_frames ));
                obj.fretTraces.Ch2.background(N) = nanmean( periphery_v_t( find( periphery_v_t>0, 100 ) ) ); % Mean of the first 100 points
                obj.fretTraces.Ch1.background(N) = nanmedian( obj.fretTraces.Ch1.int(N,min(4000,obj.fretTraces.Ch2.traceMetadata(N).endOfTrace+[0:100])) ); % Mean of the first 100 points
                
                %figure; plot( periphery_v_t, 'o' ); title(periphery_mean)
                obj.traceInterference(N,idx) = score_(idx);
               
            end
        end
        
        function removeInterference( obj )
            
            obj.fretTraces.Ch1.int_clean = obj.fretTraces.Ch1.int - repmat( obj.fretTraces.Ch1.background, obj.Ntimes, 1 )';
            info = fitgmdist( obj.traceInterference( obj.traceInterference>0.5 ), 1 );
            
            threshold = info.mu - 3*info.Sigma;
            fprintf('Using threshold of %1.2f\n',threshold)
            obj.fretTraces.Ch2.int_clean = obj.fretTraces.Ch2.int .* (obj.traceInterference>threshold)  - repmat( obj.fretTraces.Ch2.background, obj.Ntimes, 1 )';
            
            % figure; hist( log( myfret.flatten( myfret.fretTraces.Ch2.int .* (isnan(myfret.traceInterference)==0)
        end
        
        function calculateFret( obj, movType, varargin )
            
            
            if isempty(varargin); vis = 0; else vis=1; end;
            
            % Calculate all fret values %
            if ~or( strcmp(movType,'free'), strcmp(movType,'immobile') ); fprintf("\nInput error! Input either \'free\' or \'immobile\'\n\n"); return; end
            
            fret = (obj.fretTraces.Ch1.int_clean)./(obj.fretTraces.Ch1.int_clean+obj.fretTraces.Ch2.int_clean);
            if strcmp(movType,'immobile'); fret = fret .* (obj.fretTraces.Diff.Ch1.idlCnf) + fret .* (obj.fretTraces.Diff.Ch1.idlImb); end
            if strcmp(movType,'free'); fret = fret .* (obj.fretTraces.Diff.Ch1.idlFre); end
            fret(find(fret>=1 | fret<=0 ))=nan;
            fret_truncated = arrayfun( @(x) fret( x, find(  fret(x,:)>0 & isnan( fret(x,:) ) == 0 , 500 ) ), [1:size(fret,1)], 'UniformOutput', false );
            
            % Calculate prebleach F values and see which may be unlikely %
            fret_by_particle = cellfun( @(x) nanmean( x(x>0) ), fret_truncated );
            outputs = arrayfun( @(x) obj.getBleach_single( x ), [1 : obj.Ntracks] );
            
            % Ch1_ch2 prebleach F is a two column array with rows equal to
            % the number of tracks under analysis
            
            ch1_ch2_prebleach_F = [arrayfun( @(x) x.Ch1.F_prebleach, outputs )', arrayfun( @(x) x.Ch2.F_prebleach, outputs )'];
            
            % 'VariableNames', {'Min_Ch1','Max_Ch1','Min_Ch2','Max_Ch2'} )
            % Rows are each of the groups produced by dbscan
            %
            % The r variable is a first-pass filter based on the range of
            % values we believe the Ch1 and Ch2 could take
            %
            % Note that this is only used to create a sample for the Gaussian
            % fit in the next step (it is not an actual threshold on the data)
            %
            % These parameters should be inferred but for now they're hard
            % coded
            [r,~] = find( sum( or( ch1_ch2_prebleach_F<10, ch1_ch2_prebleach_F>600), 2 ) == 0 );
            
            if vis; figure; scatter( ch1_ch2_prebleach_F(r,1), ch1_ch2_prebleach_F(r,2), 50 ); end
            
            % Calculate a Gaussian fit for the thresholded data
            
            fret_by_cell_gm = fitgmdist( [ch1_ch2_prebleach_F(r,1), ch1_ch2_prebleach_F(r,2)], 1 );
            
            ll_data = pdf( fret_by_cell_gm, [ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2)] );
            result_ = ll_data > .5*10^-5;
            starttime = arrayfun( @(x) x.endOfTrace, obj.fretTraces.Ch1.traceMetadata )';
            npoints = cellfun( @(x) numel(x), fret_truncated );
            
            figure; imagesc( fret( result_, : ) )
                        
            % A mesh grid to span the values of the times in the movie data
            % and the extracted tracks under analysis
            %[fret_time,fret_trace] = meshgrid( [1:obj.Ntimes],[1:obj.Ntracks] );
            %x = fret( result_, : ).*fret_time( result_, : ); 
            %y = fret( result_, : ).*fret_trace( result_, : );
            
            if vis; figure('color','k'); 
            %scatter( ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2), 50, result_, 'filled', 'MarkerFaceAlpha', 0.5 );
        
            ax_ = arrayfun(@(x) subplot(2,2,x), [1:4] );
            subplot(2,2,1);scatter( ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2), 50, result_, 'filled', 'MarkerFaceAlpha', 0.5 );
            subplot(2,2,2);scatter( ch1_ch2_prebleach_F(find(result_==1),1), ch1_ch2_prebleach_F(find(result_==1),2), 50, starttime(find(result_==1)), 'filled', 'MarkerFaceAlpha', 0.5 );
            subplot(2,2,3);scatter( ch1_ch2_prebleach_F(find(result_==1),1), ch1_ch2_prebleach_F(find(result_==1),2), 50, npoints(find(result_==1)), 'filled', 'MarkerFaceAlpha', 0.25 );
            
            myst = suptitle( sprintf('Cell %s (%s)', obj.cellnum{1}, movType ) );
            myst.Color=[1,1,1]; myst.FontWeight='normal'; myst.FontSize= 24;
            set(gcf,'Colormap',jet(40)); 
            arrayfun(@(x) set(x,'color',[0,0,0],'XLabel',text(0,0,'Acceptor F'),'YLabel',text(0,0,'Donor F'),'XLim',[0,600],'YLim',[0,2000],'xcolor',[1,1,1],'ycolor',[1,1,1],'linewidth',4,'fontsize',18,'TickDir','out'), ax_); end
            
            % Filter by start time %
            fret = cell2mat( arrayfun( @(x) cell2mat(fret_truncated(x))', intersect( find(starttime<4000),find( result_ == 1) ), 'UniformOutput', false ) );
            
            % Create the Gaussian fit for the fret data
            mygm = fitgmdist( fret(:),1 );
            
            if vis; figure('color','w');
            bins = [0:0.025:1];
            counts = histc(fret(:),bins);
            bar( bins, counts ); hold on;
            scatter( bins, counts, 60, 'k', 'filled' );
            
            title( sprintf('Cell %s (%s)', obj.cellnum{1}, movType ), 'color', 'k', 'fontweight', 'normal', 'fontsize', 24 ); box off;
            set(gca,'color',[1,1,1],'xcolor',[0,0,0],'ycolor',[0,0,0],'linewidth',4,'fontsize',18,'TickDir','out');
            
            counts = pdf( mygm, bins' ) * mean(diff(bins)) * numel(fret(:));
            line( bins', counts, 'color', 'k', 'linewidth', 3 );
            filename = sprintf('c:\\#wes_smfret\\figures\\%s_Histogram_Cell_%s_%s.png', obj.experimentName, obj.cellnum{1}, movType );
            saveas(gcf, filename);
            end
            
            obj.fretTraces.Fret.(sprintf('Calculated_%s',movType)) = fret;
            obj.fretTraces.Fret.gmfit = mygm;
            obj.fretTraces.Fret.kept_particles_logical = result_;
            
        end
        
        function displayInfo( obj, trace_id )
            
            % FIGURE Distance between donor and acceptor %
                figure('color','w', 'position', [50         460        1400         330]); 
                subplot(1,3,1); plot( diag( pdist2( [obj.fretTraces.Ch2.x(trace_id,:)', obj.fretTraces.Ch2.y(trace_id,:)'],...
                    [obj.fretTraces.Ch1.x(trace_id,:)', obj.fretTraces.Ch1.y(trace_id,:)'] ) ), 'color', 'k' );
                title('Distance between Ch1 and Ch2')
                xlabel('Frame');ylabel('Distance (pixels)');
                xlim([-500,4000]);
                ylim_ = get(gca,'YLim');

                % Start of trace
                text( obj.fretTraces.Ch1.traceMetadata(trace_id).startOfTrace, 0.01, 'traceMetadata(n).startOfTrace','Rotation', 90, 'Color', 'k');
                line( [obj.fretTraces.Ch1.traceMetadata(trace_id).startOfTrace,obj.fretTraces.Ch1.traceMetadata(trace_id).startOfTrace], ...
                    ylim_, 'color', 'r', 'linestyle', '--' )

                text( obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace*1.2, 0.01, 'traceMetadata(n).endOfTrace', 'Rotation', 90, 'Color', 'k');
                line( [obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace,obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace], ...
                    ylim_, 'color', 'r', 'linestyle', '--' )

                text( obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace+obj.fretTraces.Ch1.traceMetadata(trace_id).lenBaseline, 0.01, ...
                    sprintf('traceMetadata(n).endOfTrace\n+traceMetadata(n).lenBaseline'), 'Rotation', 90, 'Color', 'k');
                line( [obj.fretTraces.Ch1.traceMetadata(trace_id).lenBaseline+obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace, ... 
                    obj.fretTraces.Ch1.traceMetadata(trace_id).lenBaseline+obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace], ...
                    ylim_, 'color', 'r', 'linestyle', '--' )

                subplot(1,3,2); plot( obj.fretTraces.Ch2.x(trace_id,:)', obj.fretTraces.Ch2.y(trace_id,:)');
                    hold on; plot( obj.fretTraces.Ch1.x(trace_id,:)', obj.fretTraces.Ch1.y(trace_id,:)' );
                xlabel('x position');ylabel('y position'); legend({'Ch1','Ch2'})
                title('Position');

                subplot(1,3,3); 
                color_ = zeros( 1 , size(obj.fretTraces.Fret.idlTotal,2) );
                % What is idealized
                idealized_x = find(obj.fretTraces.Fret.idlTotal(trace_id,:)==1);
                trace_x = [obj.fretTraces.Ch2.traceMetadata(trace_id).startOfTrace:obj.fretTraces.Ch2.traceMetadata(trace_id).endOfTrace];
                pre_baseline_x = setxor( [ max(1,obj.fretTraces.Ch2.traceMetadata(trace_id).startOfTrace-obj.fretTraces.Ch2.traceMetadata(trace_id).lenBackset):obj.fretTraces.Ch2.traceMetadata(trace_id).endOfTrace ], trace_x );
                post_baseline_x = [ obj.fretTraces.Ch2.traceMetadata(trace_id).endOfTrace:obj.fretTraces.Ch2.traceMetadata(trace_id).endOfTrace+obj.fretTraces.Ch2.traceMetadata(trace_id).lenBaseline ];
                
                %color_( idealized_x ) = 1;
                color_( trace_x ) = 2;
                color_( pre_baseline_x ) = 3;
                color_( post_baseline_x ) = 4;
                
                x = [1:4000];
                y = obj.fretTraces.Ch2.int(trace_id,x);
                
                x(end)=nan; y(end)=nan; color_(end) = nan;
                
                patch( x, y, color_, 'edgecolor', 'flat' );
                hold on;
                plot( idealized_x, obj.fretTraces.Ch2.int(trace_id,idealized_x), 'k--' )
                
                title('Intensity');
                ylim_ = get(gca,'YLim');
                xlim([-500,4000]);
                
                % Start of trace
                line( [obj.fretTraces.Ch1.traceMetadata(trace_id).startOfTrace,obj.fretTraces.Ch1.traceMetadata(trace_id).startOfTrace], ...
                    ylim_, 'color', 'r', 'linestyle', '--' )
                line( [obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace,obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace], ...
                    ylim_, 'color', 'r', 'linestyle', '--' )
                line( [obj.fretTraces.Ch1.traceMetadata(trace_id).lenBaseline+obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace, ... 
                    obj.fretTraces.Ch1.traceMetadata(trace_id).lenBaseline+obj.fretTraces.Ch1.traceMetadata(trace_id).endOfTrace], ...
                    ylim_, 'color', 'r', 'linestyle', '--' )
                xlim([-500,4000]);
                rng(5);
                myparula = parula(100); myparula = myparula(randperm(100),:);
                set(gcf,'colormap',myparula);
                suptitle(sprintf('Trace %i',trace_id))
                
                
            % END FIGURE Distance between donor and acceptor %
            
        end
        
        
        function output = parse_segResFinal( obj )
            
            % Parse DC-MSS output 
            
            % Each of these are segments that are at least 20 frames in
            % length, they are specified by their ids and state
            ids =       obj.fretTraces.Diff.Ch1.segResFinal(:,23); % ids column
            seg    =    obj.fretTraces.Diff.Ch1.segResFinal(:,[1:2]); % These are the indices of the segment in the frame basis
            state  =    obj.fretTraces.Diff.Ch1.segResFinal(:,3);
            Nframes =   seg(:,2)-seg(:,1);
            Ntraces =   numel(unique(ids));
            
            % Collect seg and state to create a vector for all the visited
            % diffusion states %
            tbl1 = table( ids,seg,state,Nframes, 'VariableNames', {'Ids','Seg','State','Nframes'} );
            tbl1.segs_extended = rowfun( @(id,seg,state,Nframes) repmat(state,1,seg(2)-seg(1)+1), tbl1, 'OutputFormat', 'cell' );
            tbl1.segs_xvalues = rowfun( @(id,seg,state,Nframes,x) [seg(1):seg(2)], tbl1, 'OutputFormat', 'cell' );
            
            % Create a table using the static method
            % fretTracksTable.flatten( )
            table_of_segs = varfun( @(x) { obj.flatten(x) }, tbl1, 'GroupingVariables', {'Ids'}, 'InputVariables', {'segs_extended'} );
            
            tmp = arrayfun(@(x) obj.flatten( tbl1(tbl1.Ids==x,:).segs_extended ), unique(ids), 'UniformOutput', false);
            output = tmp;
            
            % Use tmp above to calculate a matrix for DC MSS %
            DCMSSmatrix = nan( obj.Ntracks, obj.Ntimes );
            for i = 1:Ntraces
                DCMSSmatrix(i, obj.fretTraces.Ch2.traceMetadata(i).startOfTrace+[1:numel(output{i})] ) = output{i};
            end
            
            obj.DCMSSmatrix = DCMSSmatrix([1:obj.Ntracks],[1:obj.Ntimes]);
            
        end
       
        function getFretValues( obj )
            
            Ntracks = size( obj.fretTraces.Ch2.x, 1 );
            donorRange_array = nan( Ntracks, 2 );
            
            % First calculate the donor range for every track %
            for mol = 1:Ntracks

                lenBackset   = obj.fretTraces.Ch2.traceMetadata(mol).lenBackset;
                startOfTrace = obj.fretTraces.Ch2.traceMetadata(mol).startOfTrace;
                endOfTrace   = obj.fretTraces.Ch2.traceMetadata(mol).endOfTrace;
                lenBaseline  = obj.fretTraces.Ch2.traceMetadata(mol).lenBaseline;

                traceStart   = startOfTrace - lenBackset;
                traceEnd     = endOfTrace   + lenBaseline;

                donorRange_array(mol,:) = [traceStart, traceEnd];
            end
            
            % Next modify donorRange for each cell to exclude overlapping %
            inxn = {}; % Indices of interactions (ranges from 1 to numel in the session)
            diffusing = {}; % Indices of interactions (ranges from 1 to numel in the session)
            fraction_inside = [];

            for i = 1:Ntracks

                inxn{i} = find(obj.ROImatrix(i,:)<2);
                diffusing{i} = find(fretTracksObject.DCMSSmatrix(i,:)==1);
                donorrange = donorRange_array(i,1):donorRange_array(i,2);

                fraction_inside(i,:) = [ numel(inxn{i}), donorRange_array(i,2)-donorRange_array(i,1) ] ;

            end

            Ntimes = 120;

            % Calculate FRET %
            acceptor   = arrayfun( @(x) {fretTracksObject.Ch1.intensity(x,                 setxor( diffusing{i}, [donorRange_array(x,1):donorRange_array(x,2)] )      )}, [1:fretTracksObject.Ntracks], 'UniformOutput', true );
            donor      = arrayfun( @(x) {fretTracksObject.Ch2.intensity(x,                 setxor( diffusing{i}, [donorRange_array(x,1):donorRange_array(x,2)] )      )}, [1:fretTracksObject.Ntracks], 'UniformOutput', true );
            idl        = arrayfun( @(x) {fretTracksObject.fret.idlTotal(x,           setxor( diffusing{i}, [donorRange_array(x,1):donorRange_array(x,2)] )      )}, [1:fretTracksObject.Ntracks], 'UniformOutput', true );

            acceptor_filtered   = arrayfun( @(x) {fretTracksObject.Ch1.intensity(x,        setxor( [diffusing{i}, inxn{x}], [donorRange_array(x,1):donorRange_array(x,2)] )      )}, [1:fretTracksObject.Ntracks], 'UniformOutput', true );
            donor_filtered      = arrayfun( @(x) {fretTracksObject.Ch2.intensity(x,        setxor( [diffusing{i}, inxn{x}], [donorRange_array(x,1):donorRange_array(x,2)] )      )}, [1:fretTracksObject.Ntracks], 'UniformOutput', true );
            idl_filtered        = arrayfun( @(x) {fretTracksObject.fret.idlTotal(x,  setxor( [diffusing{i}, inxn{x}], [donorRange_array(x,1):donorRange_array(x,2)] )      )}, [1:fretTracksObject.Ntracks], 'UniformOutput', true );

            total = {};
            total_filtered = {};

            for i = 1:fretTracksObject.Ntracks
               total{i} = idl{i} .* (acceptor{i}./(acceptor{i}+donor{i})); 
               total{i} = total{i}( 1:min(numel(total{i}),Ntimes) );
               total{i} = total{i}(total{i}>-1);

               total_filtered{i} = idl_filtered{i} .* (acceptor_filtered{i}./(acceptor_filtered{i}+donor_filtered{i})); 
               total_filtered{i} = total_filtered{i}( 1:min(numel(total_filtered{i}),Ntimes) );
               total_filtered{i} = total_filtered{i}(total_filtered{i}>=0);

            end


            tohist = cell2mat( total ); tohist = tohist( tohist<=1 );
            tohist_filtered = cell2mat( total_filtered ); tohist_filtered = tohist_filtered( tohist_filtered<=1 );
            tohist_excluded = setxor( tohist, tohist_filtered ); tohist_excluded = tohist_excluded( tohist_excluded<=1 );

            bins = [0:0.05:1];

            tohist_scatter = histc( tohist, bins ) / sum( tohist );
            tohist_filtered_scatter = histc( tohist_filtered, bins ) / sum( tohist_filtered ) ;
            tohist_excluded_scatter = histc( tohist_excluded, bins ) / sum( tohist_excluded ) ;

            %figure('color','w','position', [450,500,1400,470] ); 
            %subplot(1,3,1); plot( bins, tohist_scatter ); xlim([0,1]); ylim([0,.25]); title('Diffusing + Idealized')
            %subplot(1,3,2); plot( bins, tohist_filtered_scatter ); xlim([0,1]); ylim([0,.25]); title('Diffusing + Idealized + Filtered')
            %subplot(1,3,3); plot( bins, tohist_excluded_scatter ); xlim([0,1]); ylim([0,.25]); 
            %title(sprintf('Excluded by filtering\n(N=%i FRET values over %i traces)',numel(tohist_excluded),fretTracksObject.Ntracks))
            %suptitle(fretTracksObject.fretfile)

            figure('color','w','position', [450,500,1400,470] ); 
            subplot(1,3,1); bar( bins, tohist_scatter ); xlim([0,1]); ylim([0,.5]); 
            % PLOTTING FIT %
            [mu,sigma] = normfit( tohist );
            Y = normpdf( bins, mu, sigma );
            Y = rescale(Y,.9*max(tohist_scatter));
            hold on; plot( bins, Y, 'k', 'linewidth', 3 );
            % END PLOTTING %
            title(sprintf('Diffusing + Idealized (mu=%1.4f)',mu))

            subplot(1,3,2); bar( bins, tohist_filtered_scatter ); xlim([0,1]); ylim([0,.5]); 
            % PLOTTING FIT %
            [mu,sigma] = normfit( tohist_filtered );
            Y = normpdf( bins, mu, sigma );
            Y = rescale(Y,.9*max(tohist_scatter));
            hold on; plot( bins, Y, 'k', 'linewidth', 3 );
            % END PLOTTING %
            title(sprintf('Diffusing + Idealized + Filtered (mu=%1.4f)',mu))

            subplot(1,3,3); bar( bins, tohist_excluded_scatter ); xlim([0,1]); ylim([0,.5]); 
            % PLOTTING FIT %
            [mu,sigma] = normfit( tohist_excluded );
            Y = normpdf( bins, mu, sigma );
            Y = rescale(Y,.9*max(tohist_scatter));
            hold on; plot( bins, Y, 'k', 'linewidth', 3 );
            % END PLOTTING %

            title(sprintf('Excluded by filtering\n(mu=%1.3f, N=%i FRET values over %i traces)',mu,numel(tohist_excluded),fretTracksObject.Ntracks))
            suptitle(fretTracksObject.fretfile)
            
            
            
        end
        
        function watch( obj, N, videoFlag, frames )
            
                mydata = obj.fretTraces.Ch2.x(N,:);
                
                figure;
                idx = [ obj.fretTraces.Ch2.traceMetadata(N).startOfTrace:obj.fretTraces.Ch2.traceMetadata(N).endOfTrace]; 
                
                if videoFlag;
                videoObj = VideoWriter( sprintf('c:\\temp\\particle_%i',N), 'MPEG-4');
                open(videoObj);
                end
                
                if isempty(frames); frames = idx; end;
                
                figure('color','w');
                for frame_ = frames;% idx

                    x=fix( obj.fretTraces.Ch2.x(N,frame_) );
                    y=fix( obj.fretTraces.Ch2.y(N,frame_) );
                    dw = 10;
                    
                    if or( isnan(x), isnan(y) )
                        continue
                    end
                    
                    imagesc( obj.images{frame_}( y+[-dw:dw], x+[-dw:dw] ), [0, 100] );
                    title( sprintf('Frame: %i',frame_) )
                    
                    if videoFlag
                    frame = getframe(gcf);
                    writeVideo(videoObj,frame);
                    end
                    
                    pause(0.1);

                end
                
                if videoFlag; close(videoObj); end;
                
        end
        
        
        function fretTraces = make_cellViewTraces( obj )
            
            fretTraces = obj.fretTraces;
            for ch = {'Ch1','Ch2'}
                fretTraces.(ch{1}).int_original = fretTraces.(ch{1}).int;
                fretTraces.(ch{1}).int = fretTraces.(ch{1}).int_clean; 
            
                fretTraces.(ch{1}) = rmfield(fretTraces.(ch{1}),'int_clean');
                fretTraces.(ch{1}).int( ~obj.fretTraces.Fret.kept_particles_logical, : ) = nan;
            end
            
        end
        
    end
    
    
    methods(Static)
        
        
        function output = findNormDistPoints( mydata, mylabels )
            
            svm = fitcsvm( mydata, mylabels );
            sv = svm.SupportVectors;
            %figure; plot( sv(:,1), sv(:,2), 'o' )

            output = polyfit( sv(:,1), sv(:,2) , 1 );
            [ xmin,xmax ] = deal( -output(2)/output(1)*100, max(mydata(:,1)) )
            vertices = [ xmin, polyval(output,xmin) ; xmax, polyval(output,xmin); xmax, polyval(output,xmax) ]

            in_poly = inpolygon( mydata(:,1), mydata(:,2), vertices(:,1), vertices(:,2) );

            figure; plot( sv(:,1), sv(:,2), 'o' ); lsline
            hold on;  plot( mydata( in_poly, 1 ), mydata( in_poly, 2 ), 'o' )
            hold on;  plot( mydata( ~in_poly, 1 ), mydata( ~in_poly, 2 ), 'ko' )

            [x1,x2] = deal( mydata( ~in_poly, 1 ), mydata( ~in_poly, 2 ) );
            GMModel = fitgmdist ( [x1,x2], 1 );

            gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
            g = gca;
            fcontour(gmPDF,[g.XLim g.YLim])
            
            y = gmPDF( mydata(:,1), mydata(:,2) );
            [ gm_x, gm_y ] = deal( linspace( min(y), max(y), 30 ), histc( y, linspace( min(y), max(y), 30 ) ) );
            %figure; bar( linspace( min(y), max(y), 20 ), histc( y, linspace( min(y), max(y), 20 ) ) );
            [~,b] = min( gm_y(1:end-1) );
            threshold_ = 1.5*10^-6;
            
            %plot( vertices(:,1), vertices(:,2), 'r--' )
            figure; scatter( mydata(:,1), mydata(:,2), 20 );
            hold on; scatter( mydata(find(y>threshold_),1), mydata(find(y>threshold_),2), 20 );

            output = find(y>threshold_);
            
        end
        
        function output = nanzscore( matrix )
           
            output = ( matrix - nanmean( matrix, 1 ) )./ nanstd( matrix, [], 1 );
            
        end
        
        
        function output = flatten( array )
            
            if isa(array,'double')
                output = array(:);
                output = output( and( isnan(output)==0, output~=0 ) );
            end
            
            if isa(array,'cell')
                cell_flat = [];
                
                for i = 1:numel(array)
                   if numel(array)>1
                       cell_flat = [cell_flat, array{i}];  
                   else
                       cell_flat = array{1};
                   end
                end
                
                output = cell_flat;
            end
            
        end
        
        function txt = myupdatefcn3(~,event_obj,obj,inputdata)
            
            % Get the data that generated the point you clicked on %
            dataIdx = get(event_obj,'DataIndex');
            theData = inputdata(dataIdx,:);% Plot it! %
            figure('color','w'); plot( inputdata(dataIdx,:) );% Return a text label for your original figure %
            title( sprintf('Idx = %i, Min = %1.2f,Max = %1.2f', dataIdx, min(theData), max(theData) ) );
            txt = sprintf('Idx = %i, Min = %1.2f,Max = %1.2f', dataIdx, min(theData), max(theData) );

        end
        
        function output = numrois( myim, mean_, std_ ) % Zthresh is calculated on a per frame basis
            
            output.NumObjects = 0;
            output.NumPixels = 0;
            output.MeanIntensity = 0;
            
            if nargin==0
                return;
            end
            
            zscore2 = @(x) (x - mean_) ./ std_;
            zscoreparam = 2;
            bwareaopen_param=2;
            
            myim_z = zscore2(myim) > zscoreparam;
            im_features = bwconncomp( myim_z );
            
            %im_features = bwconncomp( bwareaopen( zscore2( myim ) > zscoreparam, bwareaopen_param) );
            %output = numel(im_features.PixelIdxList{1});
            
            %if im_features.NumObjects==0; output.NumObjects = 0; output.NumPixels = 0; output.MeanIntensity = 0; return; end
            
            output.NumObjects = im_features.NumObjects;
            output.NumPixels = nansum( cellfun( @(x) numel(x), im_features.PixelIdxList ) );
            output.MeanIntensity = nanmean(nanmean(myim));
                
        end
        
        function output = removeNaNs( input )
           
            nanvals = find( isnan(input)==1 );
            
            x = setdiff( [1:numel(input)], nanvals );
            v = input(x);
            xq = nanvals;
            
            y = interp1( x,v,xq );
            input(xq) = y;
            
            output = input;
            
        end
        
        function auroc_value = auroc( data, split, window )
            
            if or( split-window<1, split+window>numel(data) );
                auroc_value = nan;
                return;
            end
            
            left_window = split+[-window:0]';
            right_window = split+[1:window]';
            
            left = data(left_window); Nleft = numel(left);
            right = data(right_window); Nright = numel(right);

            [min_,max_] = deal( min(data), max(data) );
            increment = (max_ - min_)/100; 
            auroc_value = [];

            coded = [zeros(Nleft,1);ones(Nright,1)]; 

                %loop through thresholds
                for i = 1:100
                    threshold = i*increment + min_;
                    result = gt([left;right],threshold);
                    %calculate true positives and false positives
                    tp(i)=sum((result==1)&(coded==1))/Nleft;
                    fp(i)=sum((result==1)&(coded==0))/Nright;
                end

                % Calculate the Riemann sum underneath the curve (AUC)
                [a,b] = sort(fp); 
                auroc_value = trapz(sort(fp),tp(b));

        end
        
        
    end
    
end

