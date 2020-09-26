function P = loadProblem(name)
path_size = split(name,'@');
path_folders = split(path_size{1},'/');
curFolder = fileparts(mfilename('fullpath'));
path = fullfile(curFolder, path_folders{:});

% check if the file exists as a mat
if exist([path '.mat'], 'file') && length(path_size) == 1
    load([path '.mat'], 'problem');
    P = problem;
elseif exist([path '.m'], 'file') && length(path_size) == 2
    [folder,name] = fileparts(path);
    prevDir = pwd;
    cd(folder);
    h = str2func(['@' name]);
    cd(prevDir);
    
    scurr = rng;
    rng(123456); % fix the seed for random generator
    P = h(str2double(path_size{2}));
    rng(scurr);
else
    error(['Problem ' name ' does not exists']);
end

P = Polytope.standardize(P);
