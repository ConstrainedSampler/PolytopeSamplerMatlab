classdef TestSuite < handle
    properties
        printFormat
        problemFilter
        problems = [];
        debug = 0
        nCores = +Inf % +Inf means default number of cores
        randomSeed = 0 % 0 means do not fix randomSeed
        testFunc % function for testing 
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
                {'extra/fit1p$', 'extra/fit2p$', ... %matlab can't solve it
                'basic/random_dense@\d\d\d\d', 'basic/random_sparse@\d\d\d', 'basic/birkhoff@\d\d\d\d', ... % problem too large
                };
            o.problemFilter.fileSizeLimit = [0 120000];
        end
        
        function test(o)
            output = TableDisplay(o.printFormat);
            fprintf(output.header());
            if isempty(o.problems)
                o.problems = problemList(o.problemFilter);
            end
            success = 0; total_time = 0;
            
            if (o.debug || o.nCores == 1)
                for k = 1:length(o.problems)
                    ret = o.runStep(k);
                    fprintf(output.print(ret));
                    
                    total_time = total_time + ret.time;
                    success = success + ret.success;
                end
            else
                if (o.nCores ~= +Inf)
                    delete(gcp('nocreate'));
                    parpool('local', o.nCores);
                end
                
                parfor k = 1:length(o.problems)
                    ret = o.runStep(k);
                    fprintf(output.print(ret));
                    
                    total_time = total_time + ret.time;
                    success = success + ret.success;
                end
            end
            
            fprintf('%i/%i success\n', success, length(o.problems))
            fprintf('Total time: %f\n', total_time)
        end
    end
    
    methods(Access = private)
        function ret = runStep(o, id)
            name = o.problems{id};
            
            if o.randomSeed ~= 0, rng(o.randomSeed); end
            warning('off', 'MATLAB:nearlySingularMatrix');
            warning('off', 'MATLAB:singularMatrix');
            warning('off', 'MATLAB:rankDeficientMatrix');
            warning('off', 'uniformtest:size');
            warning('off', 'stats:adtest:OutOfRangePLow');
            warning('off', 'stats:adtest:OutOfRangePHigh');
            warning('off', 'stats:adtest:SmallSampleSize');
            warning('off', 'stats:adtest:NotEnoughData');
            
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
end