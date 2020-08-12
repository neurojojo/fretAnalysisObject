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
        donor_images
        acceptor_images
            
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
        donorInterference
        acceptorCorrelation
        
        % Background fluorescence %
        background_by_cell
        background_grid
        background_blocksize = 25;
        
        % Interference detection parameters %
        dw = 2;
        min_corr = 0.85;
        Intensity_Distribution_Parameters = struct('mu',[20 20],'Sigma',[10 0;0 50],'CovarianceType','full'); 
        Lowpass_Parameters = [0.1,0.1,10]; % x blur, y blur, and t blur (10 frames are used for more distant interference) 
        Highpass_Parameters = [0.1,0.1,3]; % (here 3 frames are used for more recent interferences)
        template = padarray( padarray( fspecial('gaussian',3,1), 1, -5 )', 1, -5 );
        
        % Save directory %
        savedir = '';
        
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
            
            fprintf('Loading donor images...\n');
            myimgs = arrayfun( @(i) imread( sprintf('%s\\ImageData\\%s',experimentDir,img_files(i).name) ), [1:size(img_files,1)], 'UniformOutput', false );
            obj.donor_images = myimgs;
            
            fprintf('Loading acceptor images...\n');
            myimgs = arrayfun( @(i) imread( sprintf('%s\\ImageData\\%s',regexprep(experimentDir,'Ch2','Ch1'),img_files(i).name) ), [1:size(img_files,1)], 'UniformOutput', false );
            obj.acceptor_images = myimgs;
            
            obj.Ntracks = numel(obj.fretTraces.Ch2.traceMetadata);
            obj.Ntimes = numel(obj.donor_images);
            
        end
        
        function calculateBackground( obj )
            
            blocksize = obj.background_blocksize;
            Nblocks = ceil( size( obj.donor_images{1}, 1 ) / blocksize(1) );
            Nframes = obj.Ntimes;
            
            myoutput = []; 
            fprintf('Processing background for Ch1...\n');
            for i = 1:Nframes; myoutput.Ch1(:,:,i) = blockproc( double( obj.acceptor_images{i} ), [blocksize blocksize], @(x) nanmean(x.data(:)),...
                    'PadPartialBlocks', true, 'PadMethod', 'symmetric' ); end;
            
            fprintf('Processing background for Ch2...\n');
            for i = 1:Nframes; myoutput.Ch2(:,:,i) = blockproc( double( obj.donor_images{i} ), [blocksize blocksize], @(x) nanmean(x.data(:)),...
                    'PadPartialBlocks', true, 'PadMethod', 'symmetric' ); end;
            
            grid_ = reshape( [1:(Nblocks*Nblocks)], Nblocks, Nblocks );
            grid_ = blockproc( grid_, [1,1], @(x) repmat( x.data, blocksize, blocksize ) );

            myoutput.Ch1 = reshape( myoutput.Ch1, Nblocks*Nblocks, Nframes );
            myoutput.Ch2 = reshape( myoutput.Ch2, Nblocks*Nblocks, Nframes );
            
            obj.background_by_cell = myoutput;
            obj.background_grid = grid_;
            
        end
        
        % A function used to calculate correlations  %
        function getInterferences( obj, varargin )
            
            % Create a template matrix which is basically this
            %
            %        0   0   0   0   0
            %        0  0.1 0.1 0.1  0
            %        0  0.1 0.2 0.1  0
            %        0  0.1 0.1 0.1  0
            %        0   0   0   0   0
            %
            [imrows,imcols] = deal(size(obj.donor_images{1},1), size(obj.donor_images{1},2) );
            template = obj.template; 
            periphery = padarray( padarray( zeros(3,3), 1, 1 )', 1, 1 );
            
            % Placeholder matrix to be populated by correlation values %
            % --the correlations measure each tracks correspondence to the
            % Gaussian template                                        
            %
            obj.donorInterference = nan( obj.Ntracks, obj.Ntimes );
                
            for N = 1:obj.Ntracks % Run through the loop for each track %
                
                % Idx refers to all of the frames in which Ch2 (donor) is tracked %
                % This is specific to every track as they start and bleach
                % at different frames
                idx = [ obj.fretTraces.Ch2.traceMetadata(N).startOfTrace:obj.fretTraces.Ch2.traceMetadata(N).endOfTrace ]; 
                [score_,donor_image] = deal([],[]);
                
                % The tracking data x and y contains NaNs which will cause
                % the correlations to fail, so they are removed using the
                % static methods shown below 
                [nan_removed_x,nan_removed_y] = deal( fix( obj.removeNaNs( obj.fretTraces.Ch2.x(N,:) ) ), fix( obj.removeNaNs( obj.fretTraces.Ch2.y(N,:) ) ) );
                [nan_removed_x_ch1,nan_removed_y_ch1] = deal( fix( obj.removeNaNs( obj.fretTraces.Ch1.x(N,:) ) ), fix( obj.removeNaNs( obj.fretTraces.Ch1.y(N,:) ) ) );
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % This section will register each frame of a particle to the template %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dw=obj.dw;
                
                register_template = template; register_template(register_template==min(register_template)) = 0;
                for frame_ = idx
                    
                    [x,y] = deal( nan_removed_x(frame_), nan_removed_y(frame_) );   
                    obj.donor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ) = obj.register( register_template , obj.donor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ) );            
                    donor_image(:,:,frame_) = obj.donor_images{frame_}( [max(y-dw,1):min(y+dw,imcols)], [max(x-dw,1):min(x+dw,imrows)] ) - ...
                        nanmean( obj.background_by_cell.Ch2( obj.background_grid(y,x), frame_ ) );
                    
                    [x,y] = deal( nan_removed_x_ch1(frame_), nan_removed_y_ch1(frame_) );  
                    if and(x>dw,y>dw) % x,y may be zero because this channel will bleach first
                        obj.acceptor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ) = obj.register( register_template ,...
                            obj.acceptor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ) );
                        obj.acceptorCorrelation(N,frame_) = corr2( obj.acceptor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ), template );
                        acceptor_image(:,:,frame_) = obj.acceptor_images{frame_}( [max(y-dw,1):min(y+dw,imcols)], [max(x-dw,1):min(x+dw,imrows)] ) - ...
                            nanmean( obj.background_by_cell.Ch1( obj.background_grid(y,x), frame_ ) );
                    else % If the channel is indeed blank 
                        acceptor_image(:,:,frame_) = nan(5,5);
                    end
                end
                
                % Convolve the data using a 3d gaussian to capture interference as both a temporal
                % and spatial component contaminating nearby frames and pixels
                
                background_donor = imgaussfilt3(donor_image, obj.Lowpass_Parameters )+imgaussfilt3(donor_image, obj.Highpass_Parameters );
                
                % Compare each frame in the convolved movie to a template
                for frame_ = idx; 
                    score_(frame_) = corr2( background_donor(:,:,frame_), template ); 
                end
                % Fill the 'donorInterference' attribute of the object
                obj.donorInterference(N,idx) = score_(idx);
                
                sum2 = @(x) sum(sum(x));
                obj.Ch1.adjusted_intensity(N,:) = nan(1,obj.Ntimes);
                obj.Ch1.adjusted_intensity(N,idx) = arrayfun( @(x) sum2(acceptor_image([2:4],[2:4],x)), idx ); %squeeze(donor_image(3,3,idx));
                
                obj.Ch2.adjusted_intensity(N,:) = nan(1,obj.Ntimes);
                obj.Ch2.adjusted_intensity(N,idx) = arrayfun( @(x) sum2(donor_image([2:4],[2:4],x)), idx ); %squeeze(donor_image(3,3,idx));
                
            end
            
            1
            
            % Once the interferences have been calculated, create a field
            % inside the Ch1 and Ch2 structures, which is called
            % "int_clean"
            %
            % This field is a matrix of size obj.Ntracks x obj.Ntimes
            % so it can be compared to the U-track output stored in the
            % structure's "int" field
            
            % Only include frames for which donorInterference exceeds threshold (meaning a frame well-matched to the Gaussian
            % template) and then subtract off the background
            
            obj.fretTraces.Ch1.int_clean = obj.Ch1.adjusted_intensity .* (obj.donorInterference > obj.min_corr );
            obj.fretTraces.Ch2.int_clean = obj.Ch2.adjusted_intensity .* (obj.donorInterference > obj.min_corr );
            
        end
                
        function output = getBleach_single( obj, traceid, varargin )
            
            % The output of this function is only used in calculateFret to
            % assess the mean fluorescence intensity prior to bleaching
            %
            % These intensities values are calculated for each track
            % in both Ch1 and Ch2, and can then be used to construct
            % a 2D Gaussian for data quality inspection and filtering
            [bleach_Ch1, bleach_Ch2] = deal( struct(), struct() );
            
            % The intensities which are used to calibrate the Gaussian
            % Mixture Model are taken from the interference-free imaging
            % data
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
            
            % A matrix of size Ntracks x Ntimes %
            fret = (obj.fretTraces.Ch1.int_clean)./(obj.fretTraces.Ch1.int_clean+obj.fretTraces.Ch2.int_clean);
            
            % Same matrix split into free and immobile segments %
            if strcmp(movType,'immobile'); fret = fret .* (obj.fretTraces.Diff.Ch1.idlCnf) + fret .* (obj.fretTraces.Diff.Ch1.idlImb); end
            if strcmp(movType,'free'); fret = fret .* (obj.fretTraces.Diff.Ch1.idlFre); end
            
            % A matrix of size Ntracks x Ntimes %
            fret(find(fret>=1 | fret<=0 ))=nan;
            
            fret_per_particle = arrayfun( @(x) fret( x, find(  fret(x,:)>0 & isnan( fret(x,:) ) == 0 , 150 ) ), [1:size(fret,1)], 'UniformOutput', false );
            acceptorCorr_per_particle = arrayfun( @(x) obj.acceptorCorrelation( x, find(  fret(x,:)>0 & isnan( fret(x,:) ) == 0 , 150 ) ), [1:size(fret,1)], 'UniformOutput', false );
                       
            fret_full = cell2mat( cellfun( @(x) accumarray( [1:numel(x)]', x',[150,1],'',nan ), fret_per_particle , 'UniformOutput', false) )';
            acceptorCorr_full = cell2mat( cellfun( @(x) accumarray( [1:numel(x)]', x',[150,1],'',nan ), acceptorCorr_per_particle , 'UniformOutput', false) )';
            
            % Calculate prebleach F values and see which may be unlikely %
            % Ch1_ch2 prebleach F is a two column array with rows equal to
            % the number of tracks under analysis
            prebleach_values_Ch1_Ch2 = arrayfun( @(x) obj.getBleach_single( x ), [1 : obj.Ntracks] );
            ch1_ch2_prebleach_F = [arrayfun( @(x) x.Ch1.F_prebleach, prebleach_values_Ch1_Ch2 )', arrayfun( @(x) x.Ch2.F_prebleach, prebleach_values_Ch1_Ch2 )'];
            
            ch1_prebleach_F_full = repmat( ch1_ch2_prebleach_F(:,1), 1, obj.Ntimes );
            ch2_prebleach_F_full = repmat( ch1_ch2_prebleach_F(:,2), 1, obj.Ntimes );
            % fret_idxes will index the timepoint at which a fret value was
            % observed
            fret_idxes = arrayfun( @(x) find(  fret(x,:)>0 & isnan( fret(x,:) ) == 0 , 150 ), [1:size(fret,1)], 'UniformOutput', false );
            ch1_intensity =  arrayfun( @(x) obj.fretTraces.Ch1.int_clean( x, find(  fret(x,:)>0 & isnan( fret(x,:) ) == 0 , 150 ) ), [1:obj.Ntracks], 'UniformOutput', false );
            ch2_intensity =  arrayfun( @(x) obj.fretTraces.Ch2.int_clean( x, find(  fret(x,:)>0 & isnan( fret(x,:) ) == 0 , 150 ) ), [1:obj.Ntracks], 'UniformOutput', false );
            
            % Analysis on a per track basis %
            [fret_idxes, ch1_intensity, ch2_intensity] = deal( cell2mat( fret_idxes )', cell2mat( ch1_intensity )', cell2mat( ch2_intensity )' );
            
            % Ch1 and Ch2 intensity as a function of frame %
            % To reveal time dependency of Ch1/Ch2 particle intensity %
            
            ch1_intensity_time = accumarray( fret_idxes, ch1_intensity, [4000,1], @nanmean );
            ch2_intensity_time = accumarray( fret_idxes, ch2_intensity, [4000,1], @nanmean );
            figure; subplot(1,2,1); scatter( find(ch1_intensity_time~=0), ch1_intensity_time(ch1_intensity_time~=0) ); ylim([0,1000])
            subplot(1,2,2); scatter( find(ch2_intensity_time~=0), ch2_intensity_time(ch2_intensity_time~=0) ); ylim([0,1000])
            
            
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
            
            % Adjust based on the data %
            % Precondition based on very loose assumptions %
            boolean_pass = and(and(ch1_intensity>0,ch2_intensity>0),and(ch1_intensity<600,ch2_intensity<600));
            ch1_intensity( ~boolean_pass ) = nan;
            ch2_intensity( ~boolean_pass ) = nan;
            
            [ch1_mean,ch2_mean] = deal( nanmean(nanmean( obj.fretTraces.Ch1.int_clean ) ),...
                nanmean(nanmean( obj.fretTraces.Ch2.int_clean ) ) );
            
            condition = and( and(obj.fretTraces.Ch1.int_clean~=0,isnan(obj.fretTraces.Ch1.int_clean)==0),...
                 and(obj.fretTraces.Ch2.int_clean~=0,isnan(obj.fretTraces.Ch2.int_clean)==0) );
             
            S.Sigma = cov( obj.fretTraces.Ch1.int_clean( condition ), ...
                obj.fretTraces.Ch2.int_clean( condition ) );
            
            ch1_ch2_prebleach_F(or(ch1_ch2_prebleach_F<=0,ch1_ch2_prebleach_F>600))=nan; % Eliminate rows with illogical values
            %fret_by_cell_gm = fitgmdist( [ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2)], 1, 'Start', S );
            
            fret_by_cell_gm = fitgmdist( [ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2)], 1 );
            
            % Temporary %
            figure; ax_ = axes('nextplot','add');
            [x_,y_,fret_] = deal( arrayfun( @(i) nanmean(obj.fretTraces.Ch2.x(i,:)), [1:obj.Ntracks] ),...
                arrayfun( @(i) nanmean(obj.fretTraces.Ch2.y(i,:)), [1:obj.Ntracks] ),...
                nanmean(fret,2) );
            scatter( x_, y_, 60, fret_, 'filled' )
            % Temporary %
            
            ll_data = pdf( fret_by_cell_gm, [ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2)] );
            
            init_ll = 10^-5;
            potential_fret = sum(sum(~isnan(fret)));
            result_ = 0; while numel( find(result_>0) )./sum(and(~isnan(ch1_ch2_prebleach_F(:,1)),~isnan(ch1_ch2_prebleach_F(:,2)))) < .5; 
                result_ = ll_data > init_ll;
                init_ll = init_ll-0.0001;
            end
            starttime = arrayfun( @(x) x.endOfTrace, obj.fretTraces.Ch1.traceMetadata )';
            npoints = cellfun( @(x) numel(x), fret_per_particle );
            
            if vis; 
            figure('color','k'); 
            %scatter( ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2), 50, result_, 'filled', 'MarkerFaceAlpha', 0.5 );
        
            ax_ = arrayfun(@(x) subplot(2,2,x), [1:4] );
            subplot(2,2,1);scatter( ch1_ch2_prebleach_F(:,1), ch1_ch2_prebleach_F(:,2), 50, result_, 'filled', 'MarkerFaceAlpha', 0.5 );
            % Colored by starttime
            subplot(2,2,2);scatter( ch1_ch2_prebleach_F(find(result_==1),1), ch1_ch2_prebleach_F(find(result_==1),2), 50, starttime(find(result_==1)), 'filled', 'MarkerFaceAlpha', 0.5 );
            % Colored by length of track
            subplot(2,2,3);scatter( ch1_ch2_prebleach_F(find(result_==1),1), ch1_ch2_prebleach_F(find(result_==1),2), 50, npoints(find(result_==1)), 'filled', 'MarkerFaceAlpha', 0.25 );
            ch1_ch2_tmp = ch1_ch2_prebleach_F( and(starttime<1000,result_==1),: );
            subplot(2,2,4);scatter( ch1_ch2_tmp(:,1), ch1_ch2_tmp(:,2), 50, 'filled' );
            
            myst = suptitle( sprintf('Cell %s (%s)', obj.cellnum{1}, movType ) );
            myst.Color=[1,1,1]; myst.FontWeight='normal'; myst.FontSize= 24;
            set(gcf,'Colormap',jet(40)); 
            arrayfun(@(x) set(x,'color',[0,0,0],'XLabel',text(0,0,'Acceptor F'),'YLim',[0,600],'YLabel',text(0,0,'Donor F'),'xcolor',[1,1,1],'ycolor',[1,1,1],'linewidth',4,'fontsize',18,'TickDir','out'), ax_); 
            
            end
            
            % Filter by start time %
            %fret_idxes
            %ch1_intensity
            %ch2_intensity
            %fret_idxes_ = cell2mat( arrayfun( @(x) cell2mat(fret_idxes(x))',find( result_ == 1), 'UniformOutput', false  ) );
            
            fret = cell2mat( arrayfun( @(x) cell2mat(fret_per_particle(x))', intersect( find(starttime<4000),find( result_ == 1) ), 'UniformOutput', false ) );
            
            tmp = arrayfun( @(x) cell2mat(fret_per_particle(x))', intersect( find(starttime<4000),find( result_ == 1) ), 'UniformOutput', false );
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
            fprintf('Saved to %s\n',filename);
            end
            
            obj.fretTraces.Fret.(sprintf('Calculated_%s',movType)) = fret;
            obj.fretTraces.Fret.gmfit = mygm;
            obj.fretTraces.Fret.kept_particles_logical = result_;
            
        end
        
        function watch( obj, N, varargin )
            
                % The indices that correspond to the trace %
                % In the basis of [1,obj.Ntimes]           %
                % which agrees with the other matrices     %
                
                mydata = obj.fretTraces.Ch2.x(N,:);
                frames = []; % This is empty until an optional argument fills it %
                videoFlag = 0; % By default, no video is saved %
                videoFlag = sum( strcmp( varargin,'save' ) );
                if videoFlag; videoObj = VideoWriter( sprintf('c:\\temp\\particle_%i',N), 'MPEG-4'); open(videoObj); videoFlag = 1; end
                
                % Now check for options in the input arguments %
                removeInterferences = sum( strcmp( varargin,'removeInterferences' ) );
                
                if removeInterferences;
                    
                    try; cutoff = varargin{ find( strcmp( varargin, 'removeInterferences' ) == 1) + 1 }; catch; warning(sprintf("Enter a minimum correlation value after 'removeInterferences'.\n\nie: mymovie.watch(1,'removeInterferences',0.8)\n")); return; end
                    assert( isa(cutoff,'double'), sprintf("\nEnter a minimum correlation value after 'removeInterferences'.\n\nie: mymovie.watch(1,'removeInterferences',0.8)\n") );
                    assert( and(cutoff<1,cutoff>0), 'Minimum correlation must be a decimal less than 1!');
                    
                    frames = find( obj.donorInterference(N,:) > cutoff );
                    assert(~isempty(frames),sprintf('No frames met your stringest cut-off criteria. Suggested value is %1.2f\n',nanmean(obj.donorInterference(N,:))));
                    
                    %idx = [ obj.fretTraces.Ch2.traceMetadata(N).startOfTrace:obj.fretTraces.Ch2.traceMetadata(N).endOfTrace ];
                    %frames = idx( 1 + my_subset - obj.fretTraces.Ch2.traceMetadata(N).startOfTrace ); % In the basis of [1,obj.Nframes] 
                    
                    fprintf('\nSelecting a subset of frames based on interference coefficients!\n');
                    
                end
                
                % Where no correlation requirement,
                % select tracks over
                % their entire trajectory
     
                if isempty(frames); frames = [ obj.fretTraces.Ch2.traceMetadata(N).startOfTrace:obj.fretTraces.Ch2.traceMetadata(N).endOfTrace]; end;
                
                if sum( strcmp(varargin, 'frames') ); select_frames = idx( varargin{ find( strcmp( varargin, 'frames' ) == 1) + 1 } ); frames = frames(select_frames); end
                
                figure('color','w'); 
                %ax_ = arrayfun( @(x) subplot(2,1,x, 'NextPlot', 'add'), [1:2] ); 
                %plot( ax_(1), idx, obj.donorInterference(N,idx), 'k' );
                %yyaxis(ax_(1),'right'); plot( ax_(1), obj.fretTraces.Ch2.int(N,idx) )
                dw = obj.dw;
                firstframe = 1; % A logical to check if its the first frame (for intensity adjustment)
                for frame_ = frames;
                    
                    %plot( ax_(1), frame_, obj.fretTraces.Ch2.int(N,frame_), 'r.' ); pause(0.1);
                    [x,y]= deal( fix( obj.fretTraces.Ch2.x(N,frame_) ),fix( obj.fretTraces.Ch2.y(N,frame_) ) );
                    
                    if or( isnan(x), isnan(y) )
                        continue
                    end
                    
                    if firstframe; [max_intensity.Ch1,max_intensity.Ch2] = deal( max( max( obj.donor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ) ) ), ...
                            max( max( obj.acceptor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ) ) ) );
                        
                        firstframe=0; end
                    
                    %[ax_(2).XTick,ax_(2).YTick] = deal([],[]); [ax_(2).XColor,ax_(2).YColor] = deal('w','w');
                    %imagesc(  ax_(2), obj.donor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ), [0, 400] ); axis tight;
                    %title( sprintf('Frame: %i',frame_), 'parent', ax_(2) );
                    
                    subplot(1,2,1); imagesc(  gca, obj.acceptor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ), [ 0, max_intensity.Ch2 ] ); axis tight;
                    subplot(1,2,2); imagesc(  gca, obj.donor_images{frame_}( y+[-dw:dw], x+[-dw:dw] ), [ 0, max_intensity.Ch2 ] ); axis tight;
                    %set(gca,'XTick',[],'YTick',[],'Xcolor','w','YColor','w','Position',[0,0,1,1]);
                    
                    if videoFlag
                    frame = getframe(gcf);
                    writeVideo(videoObj,frame);
                    end
                    
                    pause(0.1);

                end
                
                if videoFlag; close(videoObj); end;
                
        end
        
        function saveFret( obj )
           
            obj.donor_images = 'Cleared';
            obj.acceptor_images = 'Cleared';
            filename = sprintf('%s_Cell%s_fretObject.mat',obj.experimentName,obj.cellnum{1});
            
            save( filename, 'obj' );
            
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
        
        function output = register( im1, im2 )
            cross_corr = xcorr2( im1, im2 );
            [a,b] = max(cross_corr(:));

            [x_mesh,y_mesh] = meshgrid([ (size(cross_corr,1)-1)/2:-1:-(size(cross_corr,1)-1)/2 ]);
            [x_mesh,y_mesh] = deal( x_mesh(:),y_mesh(:) );

            output = circshift( im2, [ -y_mesh(b), -x_mesh(b) ] );
        end
        
        function output = corrcoef( input1,input2 )
            
            t = corrcoef(input1,input2);
            output=t(2);
            
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

