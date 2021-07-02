% obj holds the default value
function obj = setfield(obj, source)
    key = fieldnames(obj);
    for i = 1:length(key)
        if isfield(source, key{i}) || isprop(source, key{i})
            obj.(key{i}) = source.(key{i});
        end
    end
    
    if isstruct(obj)
        key = fieldnames(source);
        for i = 1:length(key)
            obj.(key{i}) = source.(key{i});
        end
    end
end