classdef fretTraceObj
    
    properties
        ch1_int
        ch2_int
        
        ch1_ids
        ch2_ids
    end
    
    methods
        function obj = fretTraceObj( fretTraces )
            
            obj.ch1_int = fretTraces.Ch1.int;
            obj.ch2_int = fretTraces.Ch2.int;
            
            obj.ch1_ids = cell2mat( arrayfun( @(x) str2double(regexp(x.ids, '(?<=_).*', 'match')), fretTraces.Ch1.traceMetadata , 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN ) );
            obj.ch2_ids = cell2mat( arrayfun( @(x) str2double(regexp(x.ids, '(?<=_).*', 'match')), fretTraces.Ch2.traceMetadata , 'UniformOutput', false, 'ErrorHandler', @(x,y) NaN ) );
            
        end
        
        function outputArg = method1(obj,inputArg)
            
            outputArg = obj.Property1 + inputArg;
            
        end
        
        
        %Overloaded method%
        function plot( obj, traceNum )
            
            figure;
            ax = arrayfun( @(x) subplot(1,3,x,'NextPlot','add'), [1:3] );
            plot( ax(1), obj.ch1_int( traceNum, : ) );
            plot( ax(2), obj.ch2_int( traceNum, : ) );
            plot( ax(3), obj.ch1_int( traceNum, : ) ); plot( ax(3), obj.ch2_int( traceNum, : ) );
            mykids = get(gcf,'Children');
            ylim_ = [ min( arrayfun( @(x) x.YLim(1), mykids , 'UniformOutput', true) ), max( arrayfun( @(x) x.YLim(2), mykids , 'UniformOutput', true) )];
            xlim_ = [ min( arrayfun( @(x) x.XLim(1), mykids , 'UniformOutput', true) ), max( arrayfun( @(x) x.XLim(2), mykids , 'UniformOutput', true) )];
            arrayfun( @(this_ax) set(this_ax,'YLim',ylim_,'XLim',xlim_), ax );
            
        end
        
        
        
    end
end

