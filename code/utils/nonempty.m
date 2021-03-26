function b = nonempty(obj, field)
    b = isfield(obj, field) && ~isempty(obj.(field));