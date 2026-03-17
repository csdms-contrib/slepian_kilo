function setup()
    addpath(genpath('src'));

    if isfolder('external')
        addpath(genpath('external'));
    end
end
