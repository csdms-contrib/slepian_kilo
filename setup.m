function setup()
    % Add main source code
    addpath(genpath('src'));

    % Add external dependencies (if present)
    if isfolder('externals')
        addpath(genpath('externals'));
    end

    % Check for required toolbox
    assert(license('test', 'Optimization_Toolbox') == 1, ...
        'Optimization Toolbox is required.');
end
