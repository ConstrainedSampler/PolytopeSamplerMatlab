classdef SampleStorage < handle
    properties
        sampler
        opts
    end
    
    methods
        function o = SampleStorage(sampler)
            o.sampler = sampler;
            o.opts = sampler.opts.SampleStorage;
        end
        
        function o = initialize(o)
            
        end
        
        function o = propose(o)
            
        end
        
        function o = step(o)
            s = o.sampler;
            if mod(s.i, s.iterPerRecord) == 0
                s.samples(:,end+1) = s.x;
                
                len = size(s.samples, 2);
                if len >= o.opts.minNumRecords && 2 * s.iterPerRecord < s.mixingTime / o.opts.recordsPerIndependentSample
                    h = max(s.mixingTime / (o.opts.recordsPerIndependentSample * s.iterPerRecord), 1);
                    h = min(h, 2 * len/o.opts.minNumRecords);
                    h = ceil(h);
                    idx = ceil(h*(1:floor(len / h))) + (len - ceil(h * floor(len / h)));
                    s.samples = s.samples(:, :, idx);
                    s.iterPerRecord = s.iterPerRecord * h;
                end
            end
        end
        
        function o = finalize(o)
            if o.opts.rawOutput
                o.sampler.output.samples = o.sampler.samples;
            else
                o.sampler.output.samples = o.sampler.polytope.T * o.sampler.samples + o.sampler.polytope.y;
            end
        end
    end
end