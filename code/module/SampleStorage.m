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
            switch o.sampler.opts.outputFormat
                case 'raw'
                    o.sampler.output.samples = o.sampler.samples;
                case 'separate'
                    s = size(o.sampler.samples);
                    out = permute(o.sampler.samples, [2 3 1]);
                    out = reshape(out, [s(2) s(3)*s(1)]);
                    out = o.sampler.polytope.T * out + o.sampler.polytope.y;
                    n = size(o.sampler.polytope.y, 1);
                    out = reshape(out, [n s(3) s(1)]);
                    out = permute(out, [3 1 2]);
                    out = num2cell(out, [2 3])';
                    for i = 1:numel(out)
                        out{i} = reshape(out{i}, [n s(3)]);
                    end
                    o.sampler.output.samples = out;
                case 'combine'
                    s = size(o.sampler.samples);
                    out = permute(o.sampler.samples, [2 3 1]);
                    out = reshape(out, [s(2) s(3)*s(1)]);
                    o.sampler.output.samples = o.sampler.polytope.T * out + o.sampler.polytope.y;
                    o.sampler.output.summary = summary(o.sampler.output.samples);
                otherwise
                    disp('Unknown o.opts.outputFormat.')
                    o.sampler.output.samples = o.sampler.samples;
            end
        end
    end
end