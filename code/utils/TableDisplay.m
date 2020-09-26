% TODO: comments
classdef TableDisplay < handle
    properties
        format
        outputFunc
        tag
    end
    
    methods
        % format is struct
        % each field is a field in the table
        % If that field contains a string,
        %    it represents the format (according to printf)
        % otherwise
        %    it is a structure with fields
        %        format 
        %        default
        %        length
        %        label
        %        type % currently, we only support double or string
        function o = TableDisplay(format)
            fields = fieldnames(format);
            for i = 1:length(fields)
                name = fields{i};
                field = format.(fields{i});
                
                % if the field contains only a string,
                % convert it into the structure format.
                if ischar(field)
                    formattmp = field;
                    field = struct;
                    field.format = formattmp;
                end
                
                % read off the type if not specified
                if ~(isfield(field, 'type'))
                    if endsWith(field.format, 's')
                        field.type = 'string';
                    else
                        field.type = 'double';
                    end
                end
                
                % set the default if not specified
                if ~(isfield(field, 'default'))
                    if strcmp(field.type, 'string')
                        field.default = '';
                    else
                        field.default = NaN;
                    end
                end
                
                % read off the length from the format if not specified
                if ~(isfield(field, 'length'))
                    matchStr = regexp(field.format, '[0-9]*', 'match', 'once');
                    if ismissing(matchStr)
                        field.length = +Inf;
                    else
                        field.length = str2double(matchStr);
                    end
                end
                
                % use the field id as label if not specified
                if ~(isfield(field, 'label'))
                    field.label = name;
                end
                
                format.(name) = field;
            end
            o.format = format;
        end
        
        function print(o, item)
            s = "";
            fields = fieldnames(o.format);
            for i = 1:length(fields)
                name = fields{i};
                if (isfield(item, name))
                    item_i = item.(name);
                else
                    item_i = o.format.(fields{i}).default;
                end
                field = o.format.(name);
                if  strcmp(field.type, 'string') && ...
                    strlength(item_i) > field.length-1
                    item_i = extractBetween(item_i, 1, field.length-1);
                    item_i = item_i{1};
                end
                s = strcat(s, sprintf(strcat('%', field.format), item_i), ' ');
            end
            
            o.outputFunc(o.tag, '%s\n', s);
        end
        

        function header(o)
            s = '';
            fields = fieldnames(o.format);
            total_length = 0;
            for i = 1:length(fields)
                field = o.format.(fields{i});
                if field.length == +Inf
                    f = '%s';
                    total_length = total_length + strlength(field.label);
                else
                    f = strcat('%', num2str(field.length), 's');
                    total_length = total_length + field.length + 1;
                end
                s = strcat(s, sprintf(f, field.label), ' ');
            end
            total_length = total_length - 1;
            
            o.outputFunc(o.tag, '%s\n', s);
            o.outputFunc(o.tag, '%s\n', repmat('-', 1, total_length));
        end
    end
end
