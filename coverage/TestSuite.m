classdef TestSuite < handle
    properties
        printFormat
        problemFilter
        debug = 0
        nCores = +Inf % +Inf means default number of cores
        randomSeed = 0 % 0 means do not fix randomSeed
        testFunc % function for testing 
    end
    
    properties(Access = private)
        problemName
    end
    
    methods
        function o = TestSuite()
            %% Prepare the default printout info
            o.printFormat = struct;
            o.printFormat.id = '5i';
            o.printFormat.name = '50s';
            o.printFormat.success = struct('label', 'ok', 'format', '4i');
            o.printFormat.time = '10.2f';
            
            %% Setup the problem lists
            o.problemFilter = struct;
            o.problemFilter.ignoreProblems = ...
                {'netlib/pilot87$', 'netlib/fit1p$', 'netlib/fit2p$', 'netlib/qap15$', ...
                'netlib/pds_20$', 'netlib/qap12$', 'netlib/osa_60$', 'netlib/cre_b$', ...
                'netlib/cre_d$','netlib/ken_18$','netlib/dfl001$','netlib/pds_10$', ...
                'basic/random_dense@\d\d\d\d', 'basic/random_sparse@\d\d\d', 'basic/birkhoff@\d\d\d\d'};
            
            o.problemFilter.fileSizeLimit = [0 120000];
        end
        
        function test(o)
            output = TableDisplay(o.printFormat);
            output.tag = 'TestSuite';
            output.logFunc = @(tag, msg, varargin) TestSuite.output(tag, msg, varargin{:});
            output.header();
            % TODO
            l = problemList(o.problemFilter);
            success = 0; total_time = 0;
            
            o.problemName = l;
            if (o.debug || o.nCores == 1)
                for k = 1:length(l)
                    ret = o.runStep(k);
                    output.print(ret);
                    
                    total_time = total_time + ret.time;
                    success = success + ret.success;
                end
            else
                if (o.nCores ~= +Inf)
                    delete(gcp('nocreate'));
                    parpool('local', o.nCores);
                end
                
                parfor k = 1:length(l)
                    ret = o.runStep(k);
                    output.print(ret);
                    
                    total_time = total_time + ret.time;
                    success = success + ret.success;
                end
            end
            
            fprintf('%i/%i success\n', success, length(l))
            fprintf('Total time: %f\n', total_time)
        end
    end
    
    methods(Access = private)
        function ret = runStep(o, id)
            name = o.problemName{id};
            
            if o.randomSeed ~= 0, rng(o.randomSeed); end
            warning('off', 'MATLAB:nearlySingularMatrix');
            warning('off', 'MATLAB:singularMatrix');
            
            t = tic;
            if (o.debug)
                ret = o.testFunc(name);
            else
                try
                    ret = o.testFunc(name);
                catch s
                    ret = struct;
                    ret.success = 0;
                    warning('problem = %s\n%s\n\n\n', name, getReport(s,'extended'));
                end
            end
            ret.time = toc(t);
            ret.id = id;
            ret.name = name;
        end
    end
    
    methods(Static)
        function output(tag, msg, varargin)
            fprintf(msg, varargin{:})
        end
    end
end

% TODO: remove these after testing
% warning('off', 'MATLAB:nearlySingularMatrix');
% warning('off', 'MATLAB:singularMatrix');
% warning('off', 'stats:adtest:OutOfRangePLow');
% warning('off', 'stats:adtest:OutOfRangePHigh');
% warning('off', 'stats:adtest:SmallSampleSize');
% warning('off', 'stats:adtest:NotEnoughData');
% warning('off', 'unifScaleTest:size');
% warning('off', 'unifScaleTest:nonzero_grad');