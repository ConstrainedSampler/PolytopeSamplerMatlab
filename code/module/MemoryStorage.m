classdef MemoryStorage < handle
    % Module for maintaining the chain in memory
    
    properties
        opts
    end
    
    methods
        function o = MemoryStorage(s)
            o.opts = s.opts.SampleStorage;
            s.chains = zeros(s.opts.simdLen, size(s.x, 2), 0);
        end
        
        function o = step(o, s)
            if mod(s.i, s.iterPerRecord) == 0
                s.chains(:,:,end+1) = s.x;
                
                % Thin the chain if
                %   it exceeds the opts.minNumRecords and
                %   we store more samples than opts.recordsPerIndependentSample required
                len = size(s.chains, 3);
                if len >= o.opts.minNumRecords && 2 * s.iterPerRecord < s.mixingTime / o.opts.recordsPerIndependentSample
                    h = max(s.mixingTime / (o.opts.recordsPerIndependentSample * s.iterPerRecord), 1);
                    h = min(h, 2 * len/o.opts.minNumRecords);
                    h = ceil(h);
                    idx = ceil(h*(1:floor(len / h))) + (len - ceil(h * floor(len / h)));
                    s.chains = s.chains(:, :, idx);
                    s.iterPerRecord = s.iterPerRecord * h;
                end
            end
        end
        
        function o = finalize(o, s)
            if s.opts.rawOutput
                s.output.chains = s.chains;
            else
                s = size(s.chains);
                n = size(s.polytope.y, 1);
                out = permute(s.chains, [2 3 1]);
                out = reshape(out, [s(2) s(3)*s(1)]);
                out = s.problem.T * out + s.problem.y;
                out = reshape(out, [n s(3) s(1)]);
                out = squeeze(num2cell(out, [1 2]))';
                s.output.chains = out;
            end
        end
    end
end