function format_staged(files)
% Add MBeautifier path if needed
% addpath('externals/MBeautifier'); % Already done in shell script

% Load default options
options = MBeautify.defaultOptions();

% Enable “beautiful” spacing fixes
options.insertSpaceAroundOperators = true;  % Adds spaces around +, -, *, /, etc.
options.alignComments = true;                % Aligns comments
options.indentSpaces = 2;                    % Use 2-space indentation
options.useOperatorPaddingRulesFromXML = true; % Respect XML rules

for k = 1:length(files)
    file = files{k};
    try
        before = fileread(file);
        % Format file in-place
        MBeautify.formatFile(file, file, options);
        after = fileread(file);

        if ~strcmp(before, after)
            fprintf('Formatted: %s\n', file);
        else
            fprintf('Skipped (already formatted): %s\n', file);
        end
    catch ME
        fprintf('Error formatting %s: %s\n', file, ME.message);
        exit(1);
    end
end

end
