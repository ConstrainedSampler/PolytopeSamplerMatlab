classdef MemoryStorage < handle
    % Module for maintaining the chain in memory
    
    properties
        opts
    end
    
    methods
        function o = MemoryStorage(s)
            o.opts = s.opts.MemoryStorage;
            s.chains = zeros(s.opts.simdLen, s.ham.n, 0);
        end
        
        function o = step(o, s)
            if mod(s.i, s.iterPerRecord) == 0
                s.chains(:,:,end+1) = s.x;
                
                % Thin the chain if
                %   it exceeds the memory limit or more samples than opts.maxRecordsPerIndependentSample required
                if (isnan(s.mixingTime))
                    mixingTime = s.i;
                    recordMax = 1;
                else
                    mixingTime = s.mixingTime;
                    recordMax = mixingTime / o.opts.maxRecordsPerIndependentSample;

                    if (size(s.chains, 3) * s.iterPerRecord < 20 * mixingTime)
                        recordMax = recordMax / 10; % if ess < 20, records more
                    end
                end
                
                mem = numel(s.chains) * 8;
                if (mem > o.opts.memoryLimit || s.iterPerRecord < recordMax)
                   if (2 * s.iterPerRecord < mixingTime)
                       s.chains = s.chains(:, :, 2:2:end);
                       s.iterPerRecord = s.iterPerRecord * 2;
                   else
                       s.log('warning', 'Algorithm stops due to out of memory (opts.MemoryStorage.memoryLimit = %f byte).\n', o.opts.memoryLimit);
                       s.terminate = 2;
                   end
                end
            end
        end
        
        function o = finalize(o, s)
            if s.opts.rawOutput
                s.output.chains = s.chains;
            else
                N = size(s.chains,3);
                if o.opts.thinOutput
                    ess = min(effective_sample_size(s.chains), [], 1);
                else
                    ess = N * ones(size(s.chains,1),1);
                end
                
                out = [];
                for i = 1:numel(ess)
                    gap = ceil(N/ess(i));
                    out_i = s.chains(i, :, 1:gap:N);
                    out_i = reshape(out_i, [size(out_i,2) size(out_i,3)]);
                    out_i = s.problem.T * out_i + s.problem.y;
                    out = [out out_i];
                end
                s.output.chains = out;
            end
        end
    end
end