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
            s = o.sampler;
            s.samples = zeros(s.opts.simdLen, size(s.x, 2), 0);
        end
        
        function o = propose(o)
            
        end
        
        function o = step(o)
            s = o.sampler;
            if mod(s.i, s.iterPerRecord) == 0
                s.samples(:,:,end+1) = s.x;
                
                len = size(s.samples, 3);
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
        
        function o = sync(o)
            
        end
        
        function o = finalize(o)
            if o.sampler.opts.rawOutput
                o.sampler.output.samples = o.sampler.samples;
            else
                s = size(o.sampler.samples);
                n = size(o.sampler.polytope.y, 1);
                out = permute(o.sampler.samples, [2 3 1]);
                out = reshape(out, [s(2) s(3)*s(1)]);
                out = o.sampler.polytope.T * out + o.sampler.polytope.y;
                out = reshape(out, [n s(3) s(1)]);
                out = squeeze(num2cell(out, [1 2]))';
                o.sampler.output.samples = out;
            end
        end
    end
end