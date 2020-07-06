classdef fretAnalysisObject < handle
    
    % (1) Clear bleaching step (10%)
    % (2) Steady baseline
    %
    %
    
    % How to use %
    %
    % Load in data
    % myfret = fretAnalysisObject('C:\#wes_smFRET\#01Ch2\'); 
    
    % Look for multiple particles in a square
    % area of side length (obj.dw) centered on Ch2 x,y coordinates
    % 
    % myfret.getInterferences(); 
    % myfret.getBleach('Ch1'); 
    % myfret.getBleach('Ch2');
    % 
    % Combine segments from obj.fretTraces.Diff.Ch1.segResFinal
    % Column labels:
    % ids = 23; segstart = 1; segend = 2; state = 3
    %
    % myfret.parse_segResFinal();
    %
    % Calculates fret by combining  
    %
    % myfret.calculateFret();
    
    properties
        images
        fretTraces
        ROImatrix
        dcmss
        
        filename
        
        Ntracks
        Ntimes
        
        DCMSSmatrix
        
        % Details about traces %
        Ch1
        Ch2
        traceInterference = {};
        
        % Interference detection parameters %
        dw = 5;
        
        % Baseline fitting parameters %
        windowsize = 100; % VERY FLEXIBLE
        auROC_min = 0.3; % Needed to describe where there is a "stable" period
        NbaselinePts = 100; % These are the number of points needed to sample the baseline

    end
    
    methods
        
        function obj = fretAnalysisObject( experimentDir )
            
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
            
            tmp = load( fretTracesFile, 'fretTraces' );
            obj.filename = fretTracesFile;
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
        
        function output = getBleach_single( obj, ch, traceid, varargin )
            
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
            
            [start_, end_, baseline_] = deal( obj.fretTraces.(ch).traceMetadata(traceid).startOfTrace, obj.fretTraces.(ch).traceMetadata(traceid).endOfTrace, obj.fretTraces.(ch).traceMetadata(traceid).lenBaseline );
            idxes = [start_ : min(start_+end_+baseline_,numel(obj.fretTraces.(ch).int(traceid,:))) ];

            if and( isfield( obj.fretTraces.Ch1,'int_clean' ),isfield( obj.fretTraces.Ch2,'int_clean' ) );
                z = obj.fretTraces.(ch).int_clean(traceid,idxes)';
                % Interpolate nans
                z(find(isnan(z)==1)) = interp1( find(isnan(z)==0), z(find(isnan(z)==0)), find(isnan(z)==1) );
            else
                z = obj.fretTraces.(ch).int(traceid,idxes)';
            end

            % min(50,windowsize) is intended to catch situations where the 
            aurocs = arrayfun( @(x) obj.auroc( z, x, min(50,windowsize) ), [1:numel(z)] );
            [~,b] = min(aurocs);
            bleachpoint = b;

            % Baseline calculation %
            % First, calculate auROCS for pre- and post-bleaching
            [auroc_fret,auroc_baseline] = deal( aurocs(1:b), aurocs(b+1:end ) );
            
            % Second, split F data into pre- and post-bleaching
            F_prebleach = z(b-50:b);
            F_postbleach = z(b+1:end);
            
            % Third, accumulate postbleach points
            F_postbleach_steady = F_postbleach(find(abs(auroc_baseline-0.5)<(0.5-auROC_min), NbaselinePts ));

            baseline =                      nanmedian( F_postbleach_steady );
            F0 =                            nanmedian( F_prebleach );
            F0_quality =                    1 - nanstd( F_prebleach )./F0;
            baseline_quality =              1 - nanstd( F_postbleach_steady )./F0;
            
            output =                        struct();
            
            output.aurocs =                 aurocs;
            output.channel =                ch;
            
            output.parameter_windowsize =   windowsize;
            output.parameter_auROC_min =    auROC_min;
            
            output.bleachpoint_idxes =      bleachpoint;
            output.bleachpoint_frame =      idxes(bleachpoint);
            
            output.prebleach_baseline =     F0;
            output.postbleach_baseline =    baseline;
            output.prebleach_quality =      F0_quality;
            output.postbleach_quality =     baseline_quality;
            
            if ~isempty( varargin )
                figure('color','w'); 
                plot( z ); line( [bleachpoint,bleachpoint], [min(z),max(z)], 'color', 'r' )
                line( [1,bleachpoint], [F0+baseline,F0+baseline], 'color','k','linewidth',5 );
                line( [bleachpoint,numel(z)], [baseline,baseline], 'color','k','linewidth',5 );
                set(gca,'TickDir','out'); box off;
                set(gcf,'Position',[35,558,1400,240]);
            end
            
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
        
        
        % A function used to getInterferences amomg pixels %
        function getInterferences( obj )
           
            %ROImatrix = [];
            [imrows,imcols] = deal(size(obj.images{1},1), size(obj.images{1},2) );
            
            for N = 1:2%obj.Ntracks
                
                %tic;
                mydata = obj.fretTraces.Ch2.x(N,:);
                
                idx = [ obj.fretTraces.Ch2.traceMetadata(N).startOfTrace:obj.fretTraces.Ch2.traceMetadata(N).endOfTrace+obj.fretTraces.Ch2.traceMetadata(N).lenBaseline]; 
                
                for frame_ = idx
                    x=fix( obj.fretTraces.Ch2.x(N,frame_) );
                    y=fix( obj.fretTraces.Ch2.y(N,frame_) );
                    dw = obj.dw;
                    
                    try
                        background = obj.images{frame_}( [max(y-3*dw,1):min(y+3*dw,imcols)], [max(x-3*dw,1):min(x+3*dw,imrows)] );
                        [mean_,sd_] = deal( mean2(background), std2(background) );
                        if isempty(intersect( [y+[-dw:dw], x+[-dw:dw]], [0,256] ))
                            tmp_traceInterference(N).information(frame_) = obj.numrois( obj.images{frame_}( y+[-dw:dw], x+[-dw:dw] ), mean_, sd_ );
                            %ROImatrix(N,frame_) = 1;
                        else
                            tmp_traceInterference(N).information(frame_) = obj.numrois(); % Frames cannot be assessed; too close to edge of image
                            %ROImatrix(N,frame_) = nan; 
                        end
                    catch
                       %1
                    end
                    %2
                end
                
                % Once a trace is complete, store all the structured info
                % as a table for easier analysis
                obj.traceInterference{N} = struct2table( tmp_traceInterference(N).information );
                
                %t = toc;
                %fprintf('Time remaining: %1.2f\n', t/(N/obj.Ntracks))
            end
            
            %obj.ROImatrix = ROImatrix;
            
        end
        
        function removeInterference( obj )
            
            % Kmeans method %
            
            tmp_ = table2array( obj.traceInterference{1} );
            tmp = [ smooth(tmp_(:,1),5),...
                    smooth(hampel(tmp_(:,2),100,1),5),...
                    smooth(hampel(tmp_(:,3),100,1),5) ];
            tmp = obj.nanzscore(tmp);
            
            % Display using patch %
            
                
            obj.fretTraces.Ch1.int_clean = obj.fretTraces.Ch1.int;
            obj.fretTraces.Ch2.int_clean = obj.fretTraces.Ch2.int;
            
            for i = 1:size(obj.ROImatrix,1)
                ROImatrix_(i,:) = smooth( obj.ROImatrix(i,:), 5 );
            end
            obj.fretTraces.Ch1.int_clean( ROImatrix_>1 ) = nan;
            obj.fretTraces.Ch2.int_clean( ROImatrix_>1 ) = nan;
            
        end
        
        function calculateFret( obj )
            
            ROImatrix_ = obj.ROImatrix < 2;
            ROImatrix_ = ROImatrix_.*circshift( ROImatrix_, [0, 1] );%.*circshift( ROImatrix_, [0, 3] ).*circshift( ROImatrix_, [0, 4] );
            
            DCMSSmatrix_ = (obj.DCMSSmatrix==2);
            
            goodValues = ROImatrix_.*DCMSSmatrix_;
            
            Nvalues = 125;
            goodValues_array = arrayfun( @(x) find(goodValues(x,:)~=0,Nvalues), [1:size(goodValues,1)], 'UniformOutput', false );
            
            tmp = zeros( size( ROImatrix_ ) );
            for i = 1:size( ROImatrix_,1 );
                tmp = tmp + accumarray( [i*ones(numel(goodValues_array{i}),1), goodValues_array{i}'], ones(1,numel(goodValues_array{i})), size(ROImatrix_) );
            end
           
            fret = (obj.fretTraces.Ch1.int.*tmp)./(obj.fretTraces.Ch1.int.*tmp+obj.fretTraces.Ch2.int.*tmp);
            fret( or(fret<0, fret>1) ) = nan;
            
            lengths_by_particle = cellfun( @(x) numel(x), goodValues_array );
            
            fret_by_particle = nanmean(fret, 2);
            
            figure('color','w'); subplot(1,2,1);
            bins = [0:0.05:1];
            counts = histc(fret(:),bins);
            bar( bins, counts )
            title( nanmean(fret(:)) )
            subplot(1,2,2);
            plot( sort(fret_by_particle) );
            suptitle( obj.filename );
            
            datacursormode on
            fig = figure('color','w');
            dcm_obj = datacursormode(fig);
            scatter( [1:numel(lengths_by_particle)], fret_by_particle, 1+lengths_by_particle, 'filled' );
            set(dcm_obj,'UpdateFcn',{@fretAnalysisObject.myupdatefcn3, obj, fret});
            xlabel('Frame'); ylabel('FRET efficiency');
            
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
        
        function watch( obj, N, videoFlag )
            
                mydata = obj.fretTraces.Ch2.x(N,:);
                
                figure;
                idx = [ obj.fretTraces.Ch2.traceMetadata(N).startOfTrace:obj.fretTraces.Ch2.traceMetadata(N).endOfTrace+obj.fretTraces.Ch2.traceMetadata(N).lenBaseline]; 
                
                if videoFlag;
                videoObj = VideoWriter( sprintf('c:\\temp\\particle_%i',N), 'MPEG-4');
                open(videoObj);
                end
                
                figure('color','w');
                for frame_ = idx

                    x=fix( obj.fretTraces.Ch2.x(N,frame_) );
                    y=fix( obj.fretTraces.Ch2.y(N,frame_) );
                    dw = 10;
                    
                    if or( isnan(x), isnan(y) )
                        continue
                    end
                    
                    imagesc( obj.images{frame_}( y+[-dw:dw], x+[-dw:dw] ), [0, 100] );
                    title( sprintf('Particles: %i Frame: %i',obj.ROImatrix(N,frame_),frame_) )
                    
                    if videoFlag
                    frame = getframe(gcf);
                    writeVideo(videoObj,frame);
                    end
                    
                    pause(0.1);

                end
                
                if videoFlag; close(videoObj); end;
                
        end
        
    end
    
    
    methods(Static)
        
        
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
            
            output.NumObjects = nan;
            output.NumPixels = nan;
            output.MeanIntensity = nan;
            
            if nargin==0
                return;
            end
            
            zscore2 = @(x) (x - mean_) ./ std_;
            zscoreparam = 3;
            bwareaopen_param=2;
            myim_z = zscore2(myim) > zscoreparam;
            im_features = bwconncomp( myim_z );
            
            %im_features = bwconncomp( bwareaopen( zscore2( myim ) > zscoreparam, bwareaopen_param) );
            %output = numel(im_features.PixelIdxList{1});
            if im_features.NumObjects==0; return; end
            
            output.NumObjects = im_features.NumObjects;
            output.NumPixels = numel( im_features.PixelIdxList{1} );
            output.MeanIntensity = mean(myim(myim_z));
                
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

