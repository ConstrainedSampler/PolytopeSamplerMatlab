function [y] = combineSamples(x)
%y = combineSamples(x)
%combine a cell of outputs generated by multiple calls to sample into 1
%
%Input:
% x - a cell of outputs
%
%Output:
% y - one outputs

if ~iscell(x)
    y = x;
    return;
end

y = struct;
y.samplesFullDim = [];
y.samples = [];
for i = 1:length(x)
    y.samplesFullDim = [y.samplesFullDim x{i}.samplesFullDim];
    y.samples = [y.samples x{i}.samples];
end
end

