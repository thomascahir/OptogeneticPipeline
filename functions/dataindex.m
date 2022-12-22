%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Data Indexing Function                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DataIndex, FolderCount, experiment_list, StimList] = dataindex()

% Get list of all subfolders.
TopFolder = uigetdir('..\OptoPipeline\data\',"Select Data Folder.");
allSubFolders = genpath(TopFolder);
% Parse into a cell array.
remain = allSubFolders;
DataIndex = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	DataIndex = [DataIndex singleSubFolder];
end
DataIndex = transpose(DataIndex);
FolderCount = length(DataIndex);

%%%%%%%%%% BUILD EXPERIMENT LIST
% Get a list of all files and folders in this folder.
files = dir(TopFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subDirs = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
experiment_list = {subDirs(3:end).name};

%%%%%%%%%% BUILD SUBEXPERIMENT / STIMULUS TYPE LIST
StimFolders = string(DataIndex(1))+'\'+string((experiment_list(1)));
% Get a list of all files and folders in this folder.
files2 = dir(StimFolders);
% Get a logical vector that tells which is a directory.
dirFlags2 = [files2.isdir];
% Extract only those that are directories.
subDirs2 = files2(dirFlags2); % A structure with extra info.
% Get only the folder names into a cell array.
StimList = {subDirs2(3:end).name};

end

