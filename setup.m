function setup()
    % Add main source code
    addpath(genpath('src'));

    % Add external dependencies (if present)
    addpath(genpath('externals'));

    % Check for required toolboxes
    assert(license('test', 'Optimization_Toolbox') == 1, ...
        'Optimization Toolbox is required.');
end
