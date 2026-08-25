function [subjectIds, indexByList, unmatchedIds, normalizedIds] = groupSIFT_matchSubjects(fileNameLists)
% groupSIFT_matchSubjects() - Match and order subjects using fileNameList.
%
% Usage:
%   [subjectIds, indexByList, unmatchedIds, normalizedIds] = ...
%       groupSIFT_matchSubjects({fileNameList1, fileNameList2, ...});
%
% Subject identity is the filename without its path and final extension.
% Matching is case-insensitive and exact; condition labels are never
% guessed or stripped. The order of subjectIds follows the first list.

if ~iscell(fileNameLists) || isempty(fileNameLists)
    error('groupSIFT:InvalidFileNameLists', ...
        'Input must be a non-empty cell array of fileNameList variables.');
end

numberOfLists = numel(fileNameLists);
normalizedIds = cell(1, numberOfLists);
displayIds = cell(1, numberOfLists);

for listIdx = 1:numberOfLists
    currentList = fileNameLists{listIdx};
    if ischar(currentList)
        currentList = cellstr(currentList);
    elseif isstring(currentList)
        currentList = cellstr(currentList(:));
    end
    if ~iscell(currentList) || isempty(currentList)
        error('groupSIFT:InvalidFileNameList', ...
            'fileNameList %d must be a non-empty cell array of filenames.', listIdx);
    end

    currentList = currentList(:);
    normalizedIds{listIdx} = cell(size(currentList));
    displayIds{listIdx} = cell(size(currentList));
    for subjectIdx = 1:numel(currentList)
        currentName = currentList{subjectIdx};
        if isstring(currentName) && isscalar(currentName)
            currentName = char(currentName);
        end
        if ~ischar(currentName) || isempty(strtrim(currentName))
            error('groupSIFT:InvalidSubjectFilename', ...
                'Every fileNameList entry must be a non-empty character vector or string scalar.');
        end
        [~, baseName, ~] = fileparts(strtrim(currentName));
        if isempty(baseName)
            error('groupSIFT:InvalidSubjectFilename', ...
                'Could not extract a subject ID from "%s".', currentName);
        end
        displayIds{listIdx}{subjectIdx} = baseName;
        normalizedIds{listIdx}{subjectIdx} = lower(baseName);
    end

    if numel(unique(normalizedIds{listIdx})) ~= numel(normalizedIds{listIdx})
        error('groupSIFT:DuplicateSubjectId', ...
            'fileNameList %d contains duplicate subject IDs after normalization.', listIdx);
    end
end

commonMask = true(size(normalizedIds{1}));
for listIdx = 2:numberOfLists
    commonMask = commonMask & ismember(normalizedIds{1}, normalizedIds{listIdx});
end
commonNormalizedIds = normalizedIds{1}(commonMask);
subjectIds = displayIds{1}(commonMask);

indexByList = cell(1, numberOfLists);
indexByList{1} = find(commonMask);
for listIdx = 2:numberOfLists
    [isPresent, location] = ismember(commonNormalizedIds, normalizedIds{listIdx});
    if ~all(isPresent)
        error('groupSIFT:InternalSubjectMatchFailure', ...
            'Subject matching failed unexpectedly for fileNameList %d.', listIdx);
    end
    indexByList{listIdx} = location;
end

unmatchedIds = cell(1, numberOfLists);
for listIdx = 1:numberOfLists
    unmatchedMask = ~ismember(normalizedIds{listIdx}, commonNormalizedIds);
    unmatchedIds{listIdx} = displayIds{listIdx}(unmatchedMask);
end
end
