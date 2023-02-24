function IDFrameSet = ImportRawFiles(filename, pathname)

inputFormat = 'raw_cmbu';
outputFormat = 'raw10';

if iscell(filename)
    FN = filename;
else
    FN{1} = filename;
end
numFiles = length(FN);

%% calculation
for j = 1:numFiles
    %Results{j}.fullFileName = [pathname,FN{j}];
    %Results{j}.fileName = FN{j};
    % read in rawq data
    tmpID = getFile([pathname,FN{j}], inputFormat, outputFormat);
%     if (size(tmpID, 1) == 1314)
%         tmpID = tmpID(3:end, :);
%     end
    IDFrameSet(:,:,j) = double(tmpID);
    %IDFrameSet(:,:,j) = double(IDFrameSet(:,:,j));
end
