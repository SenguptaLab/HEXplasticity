%% Script to sequentially run basic analysis scripts:

function wormTrackers
useDefaultBinSett = true; %make false to adjust bin settings for analyses

setts = getComputerSettings;

folder_name = uigetdir(setts.searchPath,'Select Folder for Analysis');


pathall = [';' genpath(folder_name)];

foldersSep = strfind(pathall,';');
movFiles = {};
mCount = 0;
for i = 1:length(foldersSep)-1
    p = pathall(foldersSep(i)+1:foldersSep(i+1)-1);
    movF = dir([p filesep '*_worms.avi']);
    if ~isempty(movF)
        for j = 1:length(movF)
            mCount = mCount +1;
            movFiles{mCount} = [p filesep movF(j).name]; %#ok<*AGROW>
        end
    end
end
%% now run for each movie
defaultBinSettings = [0.5 10];

if useDefaultBinSett
    binSettings = defaultBinSettings;
else
    
    answer = dagetnum({'density plot: Time Bin (min):','density plot: Area Bin (pixels):'}, ...
        defaultBinSettings);
    binSettings(1) = answer(1).num;  % minutes per bin
    binSettings(2) = answer(2).num; % in pixels
    
end

for iMov = 1:mCount
    close all
    try
        matFile = trackWorms(movFiles{iMov});
        %     matFile = findMatFile(movFiles{iMov});
        segmenttracks_fast(matFile);
        wormDensity_fast(matFile,binSettings);
        residencyPlots_fast(matFile,binSettings);
        ethogram_fast(matFile);
      %  dirkstates_fast(matFile);
        
    catch me
        warning(['This file could not be analyzed: ' movFiles{iMov} ': ' me.message]);
        continue
    end
    display([datestr(now) '  success: ' movFiles{iMov}])
end
