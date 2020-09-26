function l = problemList(options)
if ~exist('options', 'var'), options = []; end
problems = {};

% Read files with size within the range (bytes)
default.fileSizeLimit = [0 +Inf];

%default.folders = {'basic', 'metabolic', 'metabolic2', 'netlib'};
default.folders = {'basic', 'metabolic', 'netlib'};

default.ignoreProblems = {};

default.generateDimensions = [10 100 1000 10000];

o = setfield(default, options);

curFolder = fileparts(mfilename('fullpath'));

for j = 1:length(o.folders)
    files = dir(fullfile(curFolder, o.folders{j}, '*.m*'));
    for k = 1:length(files)
        file = fullfile(files(k).folder, files(k).name);
        [~,name,ext] = fileparts(file);
        name = strcat(o.folders{j}, '/', name);
        
        if strcmp(ext, '.mat') == 1
            s = dir(file);
            if (s.bytes < o.fileSizeLimit(1) || s.bytes > o.fileSizeLimit(2))
                continue;
            end
            problems{end+1} = name;
        elseif strcmp(ext, '.m') == 1
            for l = 1:length(o.generateDimensions)
                problems{end+1} = [name '@' num2str(o.generateDimensions(l))];
            end
        end
    end
end

% removed ignored problems
l = {};
for j = 1:length(problems)
    name = problems{j};
    ignored = 0;
    for k = 1:length(o.ignoreProblems)
        if ~isempty(regexp(name, o.ignoreProblems{k}))
            ignored = 1;
            break;
        end
    end
    if ignored == 0
        l{end+1} = name;
    end
end
end
