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
        
        % Interference detection parameters %
        dw = 2;
        min_corr = 0.7;
        Intensity_Distribution_Parameters = struct('mu',[200 200],'Sigma',[20 0;0 20],'CovarianceType','full'); 
        Lowpass_Parameters = [0.1,0.1,10]; % x blur, y blur, and t blur (10 frames are used for more distant interference) 
        Highpass_Parameters = [0.1,0.1,3]; % (here 3 frames are used for more recent interferences)
                
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
        
        % A function used to calculate correlations  %
        function getInterferences( obj )
            
            % Create a template matrix which is basically this
            %
            %        -5 -5   -5  -5 -5
            %        -5 0.1 0.1 0.1 -5
            %        -5 0.1 0.2 0.1 -5
            %        -5 0.1 0.1 0.1 -5
            %        -5 -5   -5  -5 -5
            %
            [imrows,imcols] = deal(size(obj.images{1},1), size(obj.images{1},2) );
            template = padarray( padarray( fspecial('gaussian',3,1), 1, -5 )', 1, -5 );
            periphery = padarray( padarray( zeros(3,3), 1, 1 )', 1, 1 );
            
            % Placeholder matrix to be populated by correlation values %
            % --the correlations measure each tracks correspondence to the
            % Gaussian template                                        
            %
            obj.traceInterference = nan( obj.Ntracks, obj.Ntimes );
                
            for N = 1:obj.Ntracks % Run through the loop for each track %
                
                % Idx refers to all of the frames in which Ch2 (donor) is tracked %
                % This is specific to every track as they start and bleach
                % at different frames
                idx = [ obj.fretTraces.Ch2.traceMetadata(N).startOfTrace:obj.fretTraces.Ch2.traceMetadata(N).endOfTrace ]; 
                [score_,background] = deal([],[]);
                
                % The tracking data x and y contains NaNs which will cause
                % the correlations to fail, so they are removed using the
                % static methods shown below 
                [nan_removed_x,nan_removed_y] = deal( fix( obj.removeNaNs( obj.fretTraces.Ch2.x(N,:) ) ), fix( obj.removeNaNs( obj.fretTraces.Ch2.y(N,:) ) ) );
                
                for frame_ = idx    
                    [x,y] = deal( nan_removed_x(frame_), nan_removed_y(frame_) );               
                    background(:,:,frame_) = obj.images{frame_}( [max(y-obj.dw,1):min(y+obj.dw,imcols)], [max(x-obj.dw,1):min(x+obj.dw,imrows)] ); 
                end
                
                % Convolve the data using a 3d gaussian to capture interference as both a temporal
                % and spatial component contaminating nearby frames and pixels
                
                background2 = imgaussfilt3(background, obj.Lowpass_Parameters )+imgaussfilt3(background, obj.Highpass_Parameters );
                
                % Compare each frame in the convolved movie to a template
                for frame_ = idx; score_(frame_) = corr2( background2(:,:,frame_), template ); end
                
                % Fill the 'traceInterference' attribute of the object
                obj.traceInterference(N,idx) = score_(idx);
                
                % To find the background values for Ch1 and Ch2, we look at
                % 100 frames immediately following the disappearance of the
                % track (these background values are used to produce a baseline near zero after
                % bleaching has occurred)
                obj.fretTraces.Ch2.background(N) = nanmedian( obj.fretTraces.Ch2.int(N,min(4000,obj.fretTraces.Ch2.traceMetadata(N).endOfTrace+[0:500])) ); 
                obj.fretTraces.Ch1.background(N) = nanmedian( obj.fretTraces.Ch1.int(N,min(4000,obj.fretTraces.Ch1.traceMetadata(N).endOfTrace+[0:500])) ); 

            end
            
            % Once the interferences have been calculated, create a field
            % inside the Ch1 and Ch2 structures, which is called
            % "int_clean"
            %
            % This field is a matrix of size obj.Ntracks x obj.Ntimes
            % so it can be compared to the U-track output stored in the
            % structure's "int" field
            
            % Only include frames for which traceInterference exceeds threshold (meaning a frame well-matched to the Gaussian
            % template) and then subtract off the background
            
            obj.fretTraces.Ch1.int_clean = obj.fretTraces.Ch1.int - repmat( obj.fretTraces.Ch1.background, obj.Ntimes, 1 )';
            obj.fretTraces.Ch2.int_clean = obj.fretTraces.Ch2.int - repmat( obj.fretTraces.Ch2.background, obj.Ntimes, 1 )'; 
            obj.fretTraces.Ch2.int_clean = obj.fretTraces.Ch2.int_clean .* (obj.traceInterference > obj.min_corr );
            
        end
                
        function output = getBleach_single( obj, traceid, varargin )
            
            % The output of this function is only used in calculateFret to
            % assess the mean fluorescence intensity prior to bleaching
            %
            % These intensities values are calculated for each track
            % in both Ch1 and Ch2, and can then be used to construct
            % a 2D Gaussian for data quality inspection and filtering
            [bleach_Ch1, bleach_Ch2] = deal( struct(), struct() );
            
            [bleach_Ch1.start_, bleach_Ch1.end_, bleach_Ch1.baseline_] = deal( obj.fretTraces.Ch1.traceMetadata(traceid).startOfTrace, obj.fretTraces.Ch1.traceMetadata(traceid).endOfTrace, obj.fretTraces.Ch1.traceMetadata(traceid).lenBaseline );
            [bleach_Ch2.start_, bleach_Ch2.end_, bleach_Ch2.baseline_] = deal( obj.fretTraces.Ch2.traceMetadata(traceid).startOfTrace, obj.fretTraces.Ch2.traceMetadata(traceid).endOfTrace, obj.fretTraces.Ch2.traceMetadata(traceid).lenBaseline );
            
            [z1,z2] = deal( obj.fretTraces.Ch1.int_clean(traceid,:), obj.fretTraces.Ch2.int_clean(traceid,:) ); 
            z2(z2==0)=nan; % Removing interferences creates 0
            
            [bleach_Ch1.F_prebleach,bleach_Ch2.F_prebleach] = deal( nanmean( z1(bleach_Ch1.start_ : bleach_Ch1.end_) ), nanmean( z2(bleach_Ch2.start_ : bleach_Ch2.end_) ) );
            [output.Ch1,output.Ch2] = deal( bleach_Ch1, bleach_Ch2 );
            
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
            
            % Calculate a Gaussian fit for the thresholded data
            
            S = obj.Intensity_Distribution_Parameters;
            ch1_ch2_prebleach_F(ch1_ch2_prebleach_F<=0)=nan; % Eliminate rows with illogical values
            fret_by_cell_gm = fitgmdist( [ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2)], 1, 'Start', S );
            
            ll_data = pdf( fret_by_cell_gm, [ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2)] );
            result_ = ll_data > .5*10^-5;
            starttime = arrayfun( @(x) x.endOfTrace, obj.fretTraces.Ch1.traceMetadata )';
            npoints = cellfun( @(x) numel(x), fret_truncated );
            
            if vis; 
            figure('color','k'); 
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
        
        function output = nanzscore( matrix )
           
            output = ( matrix - nanmean( matrix, 1 ) )./ nanstd( matrix, [], 1 );
            
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
        
    end
    
end

