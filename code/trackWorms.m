function outputFile = trackWorms(file)

%% Defaults for worm tracking! From Dirks defaults
trackSett.MinWormArea = 0.4;          % Min area for object to be a valid worm (0.5)
trackSett.MaxWormArea = 1.7;          % Max area for single worm
trackSett.MinDistance = 15;           % Min Distance for connecting a new worm to an existing track (in pixels)
trackSett.SizeChangeThreshold = 100;  % Max size change between frames (in pixels)
trackSett.MinTrackLength = 10;        % Min Length of valid track (in frames)
trackSett.PlotRate = 200;         % Display tracking results every 'PlotFrameRate' frames
trackSett.DefaultSize  = 200;        % Default worm size (added by Till)
trackSett.LevelEstimatePerFrame  = 200;        % Estimate levels for every 200 frames and choose best value (added by Till)
trackSett.StoreImage = true;  % Will store an image of each work for each frame. Only for debugging.

setts = getComputerSettings;

if nargin <1
    [movieName, pathName] = uigetfile([setts.searchPath '*worms.avi'], 'Please choose worm movie');
else
    [pathName,movieName,extension] = fileparts(file);
    movieName = [movieName extension];
end
%%
shortNameTmp = movieName;
if ~isempty(strfind(shortNameTmp,'.avi'))
    shortNameTmp(strfind(shortNameTmp,'.avi'):end) = '';
    if ~isempty(strfind(shortNameTmp,'_worms'))
        shortNameTmp(strfind(shortNameTmp,'_worms'):end) = '';
    end
    if strcmp(pathName(end),filesep), pathName = pathName(1:end-1);end
    filsepFound = strfind(pathName,filesep);
    arenaName = pathName(filsepFound(end)+1:end);
%     if ~isempty(strfind(shortName,[arenaName '_']))
%         shortName(1:length(arenaName)+1) = '';
%     end


end
shortName = arenaName;
%%
if exist([pathName filesep shortName '.mat'],'file')
    dataLoad = load([pathName filesep shortName '.mat']);
    outputFile = [pathName filesep shortName '.mat'];
else
    error('run:trackWorms','run preprocess first');
    %     matFilesInDir = dir([pathName filesep '*.mat']);
    %     if length(matFilesInDir)==1
    %         outputFile = [pathName filesep matFilesInDir(1).name];
    %         dataLoad = load(outputFile);
    %     end
end
if isfield(dataLoad,'data') % old style
    data = dataLoad.data;
    if isfield(data,'data') % old style
        data.preprocess = data.data;
    elseif isfield(data,'pixelSize')
       data.preprocess = data;
    end
elseif isfield(dataLoad,'preprocess') % new style
    data  = dataLoad ;
end

%% ANALYZE FRAGMENTS
%----------------------
WTFigH = findobj('Tag', ['WTFIG_' arenaName]);
if isempty(WTFigH)
    WTFigH = figure('NumberTitle', 'off', 'Tag', ['WTFIG_' arenaName]);colormap('gray')
else
    figure(WTFigH);clf;colormap('gray')
end
%FragmentFrames = 100;
% for Fragment = 1:File(MovieNum).Fragments
%     FragmentTrackerTH(Mov,File,FragmentSaveNames,FragmentFrames,Fragment,Mask,TrackDye,WTFigH);
% end
movH = VideoReader([pathName filesep movieName]);

trackSett.nrFrames = movH.NumberOfFrames;
trackSett.endFrame = trackSett.nrFrames-1 ;
trackSett.startFrame = 1;
mov = read(movH,[trackSett.startFrame trackSett.endFrame]);
mov = single(squeeze(mov(:,:,1,:)));
mov = mov-median(mov(:));
%mov = single(mov);

narrowFilter = single(fspecial('average',3)); %default narrow filter is 3

mov= convn(mov,narrowFilter,'same');
%%
if strcmpi(data.preprocess.userInput.concMode,'pulse');
    trackSett.pulseOn = true;
    pulseMovieName = strrep(movieName,'worms','pulse');
    movDyeH = VideoReader([ pathName filesep pulseMovieName]);
    movDye = read(movDyeH,[trackSett.startFrame trackSett.endFrame]);
    movDye = logical(squeeze(movDye(:,:,1,:)));
    movSize = [size(movDye,1) size(movDye,2)];
else
    trackSett.pulseOn = false;
end
%%
framesForLevel =  trackSett.startFrame: trackSett.LevelEstimatePerFrame : trackSett.endFrame;
levelNr = length(framesForLevel);
levelForFrame = NaN(1,levelNr);
animalPixPerFrame = NaN(1,levelNr);
for i = 1:levelNr
    [levelForFrame(i),animalPixPerFrame(i)] = autoThreshold_filt(mov(:,:,framesForLevel(i)),trackSett);
end

trackSett.levelEstimation.allLevels = levelForFrame;
trackSett.levelEstimation.animalPixPerFrame = animalPixPerFrame;
trackSett.levelEstimation.framesForLevel = framesForLevel;

trackSett.Level = median(levelForFrame);
trackSett.AnimalPix = 350;%median(animalPixPerFrame);
warning('temporarly switched auto AnimalPix off and choose 350')
 Settings = [NaN, trackSett.Level, trackSett.AnimalPix];
%%
% background = getbackground(movH);
% imagesc(background);colormap('gray');colorbar
% title('Almost all white(between 240 and 255)')
% suptitle('will not substract background')
% trackSett = getLeveLAndAnimalsize(mov(:,:,1:trackSett.LevelEstimatePerFrame:end),trackSett);
% imagesc(background);colormap('gray');colorbar
% title('Almost all white(between 240 and 255)')
% suptitle('will not substract background')

%%
mov = (mov <= trackSett.Level);


MinWormArea = trackSett.MinWormArea;          % Min area for object to be a valid worm
MaxWormArea =trackSett.MaxWormArea;          % Max area for single worm
MinDistance  = trackSett.MinDistance ;           % Min Distance for connecting a new worm to an existing track (in pixels)
SizeChangeThreshold = trackSett.SizeChangeThreshold ;  % Max size change between frames (in pixels)
MinTrackLength = trackSett.MinTrackLength;        % Min Length of valid track (in frames)
PlotFrameRate = trackSett.PlotRate ;         % Display tracking results every 'PlotFrameRate' frames
AnimalPix = trackSett.AnimalPix;

% Initialize variables
tracks = [];
inactiveTracks = [];

%%
nrFrames = size(mov,3);
frames = 1:nrFrames;

hProgress = waitbar(0,'message');

tic
allData = cell(1,nrFrames); %cell arrays can be filled much easier than a matrix array

for Frame = frames
    %%
    if ~mod(Frame,10)
        waitbar(Frame/nrFrames,hProgress,['Analysing frame ' num2str(Frame) ' (total: ' num2str(nrFrames) ')']);
    end
    
    % Segment objects
    L = bwlabel(mov(:,:,Frame));
    STATS = regionprops  (L, {'Area', 'Centroid', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Image', 'BoundingBox'});
    %%
    % Identify all worms by size, get their centroid coordinates
    WormIndices = find([STATS.Area] > MinWormArea*AnimalPix & [STATS.Area] < MaxWormArea*AnimalPix);
    
    NumWorms = length(WormIndices);
    if isempty(NumWorms)
        break
    end
    
    WormCentroids = [STATS(WormIndices).Centroid];
    WormCoordinates = reshape(WormCentroids,2,NumWorms)'; 
    WormSizes = [STATS(WormIndices).Area];
    WormEccentricities = [STATS(WormIndices).Eccentricity];
    WormMajorAxes = [STATS(WormIndices).MajorAxisLength];
    WormMinorAxes = [STATS(WormIndices).MinorAxisLength];
    WormOrientation = [STATS(WormIndices).Orientation];
    WormBox = [STATS(WormIndices).BoundingBox];
    WormBoundingBox = [WormBox(1:4:4*NumWorms)', WormBox(2:4:4*NumWorms)'];
    if trackSett.StoreImage
        for wi = 1:length(WormIndices)
            WormImage(wi).Image = STATS(WormIndices(wi)).Image;
        end
    end
    if trackSett.pulseOn
         wormLocationIndex = sub2ind(movSize,round(WormCoordinates(:,2)),round(WormCoordinates(:,1)));
         dyeMap = movDye(:,:,Frame);
        WormPulse = dyeMap(wormLocationIndex);
     else
         WormPulse = false(NumWorms,1);
     end
    allData{Frame} = cat(2,zeros(length(STATS),1)+Frame,cat(1,STATS(:).Area),cat(1,STATS(:).Centroid), cat(1,STATS(:).Eccentricity));
    % Track worms
    % -----------
    if ~isempty(tracks)
        ActiveTracks = find([tracks.Active]);
    else
        ActiveTracks = [];
    end
    
    % Update active tracks with new coordinates
    for i = 1:length(ActiveTracks)
        DistanceX = WormCoordinates(:,1) - tracks(ActiveTracks(i)).LastCoordinates(1);
        DistanceY = WormCoordinates(:,2) - tracks(ActiveTracks(i)).LastCoordinates(2);
        Distance = sqrt(DistanceX.^2 + DistanceY.^2);
        [MinVal, MinIndex] = min(Distance);
        if (MinVal <= MinDistance) && (abs(WormSizes(MinIndex) - tracks(ActiveTracks(i)).LastSize) < SizeChangeThreshold)
            tracks(ActiveTracks(i)).Path = [tracks(ActiveTracks(i)).Path; WormCoordinates(MinIndex, :)];
            tracks(ActiveTracks(i)).LastCoordinates = WormCoordinates(MinIndex, :);
            tracks(ActiveTracks(i)).Frames = [tracks(ActiveTracks(i)).Frames, Frame];
            tracks(ActiveTracks(i)).Size = [tracks(ActiveTracks(i)).Size, WormSizes(MinIndex)];
            tracks(ActiveTracks(i)).LastSize = WormSizes(MinIndex);
            tracks(ActiveTracks(i)).Eccentricity = [tracks(ActiveTracks(i)).Eccentricity, WormEccentricities(MinIndex)];
            tracks(ActiveTracks(i)).MajorAxes = [tracks(ActiveTracks(i)).MajorAxes, WormMajorAxes(MinIndex)];
            tracks(ActiveTracks(i)).MinorAxes = [tracks(ActiveTracks(i)).MinorAxes, WormMinorAxes(MinIndex)];
            tracks(ActiveTracks(i)).Orientation = [tracks(ActiveTracks(i)).Orientation, WormOrientation(MinIndex)];
            tracks(ActiveTracks(i)).Box = [tracks(ActiveTracks(i)).Box; WormBoundingBox(MinIndex,:)];
            TrackFrameNum = length(tracks(ActiveTracks(i)).Size);
            tracks(ActiveTracks(i)).Frame(TrackFrameNum).Image = WormImage(MinIndex).Image;
            tracks(ActiveTracks(i)).Pulse = [tracks(ActiveTracks(i)).Pulse,  WormPulse(MinIndex)];
             %tracks(ActiveTracks(i)).Image = {tracks(ActiveTracks(i)).Image,WormPulse(MinIndex);
            WormCoordinates(MinIndex,:) = NaN;
        else
            tracks(ActiveTracks(i)).Active = 0;
            if length(tracks(ActiveTracks(i)).Frames) < MinTrackLength
                tracks(ActiveTracks(i)) = [];
                ActiveTracks = ActiveTracks - 1;
            end
        end
    end
    
    % Start new tracks for coordinates not assigned to existing tracks
    NumTracks = length(tracks);
    for i = 1:length(WormCoordinates(:,1))
        Index = NumTracks + i;
        tracks(Index).Active = 1;
        tracks(Index).Path = WormCoordinates(i,:);
        tracks(Index).LastCoordinates = WormCoordinates(i,:);
        tracks(Index).Frames = Frame;
        tracks(Index).Size = WormSizes(i);
        tracks(Index).LastSize = WormSizes(i);
        tracks(Index).Eccentricity = WormEccentricities(i);
        tracks(Index).MajorAxes = WormMajorAxes(i);
        tracks(Index).MinorAxes = WormMinorAxes(i);
        tracks(Index).Orientation = WormOrientation(i);
        tracks(Index).Box = WormBoundingBox(i,:);
        tracks(Index).Frame(1).Image = WormImage(i).Image;
        tracks(Index).Pulse = WormPulse(i);
    end
    
    
    %% Display every PlotFrameRate'th frame
    if ~mod(Frame, PlotFrameRate)
        t1 = toc;
        PlotFrame(WTFigH, ~mov(:,:,Frame), tracks);

        t2 = toc; tic;
        fps = PlotFrameRate/t1;
        fprintf('\nFrame: %5d - Time: %1.3f fps (%1.2f s) Level: (%1.3f/%3d)',Frame,fps,t2-t1,trackSett.Level,trackSett.AnimalPix)
    end
    
    
end
allData = cat(1,allData{:}); 
close (hProgress);
%
% Get rid of invalid tracks
tracks = cat(2,tracks,inactiveTracks);

DeleteTracks = [];
for i = 1:length(tracks)
    if length(tracks(i).Frames) < MinTrackLength
        DeleteTracks = [DeleteTracks, i];
    end
end
tracks(DeleteTracks) = [];
%%

expData.PixelSize = data.preprocess.pixelSize;
expData.ArenaSize = [size(mov,2) size(mov,1)];
expData.FrameRate = data.preprocess.userInput.framerate;
expData.TrackTime = datestr(now);
expData.TrackedFrames = nrFrames;
expData.TrackStats.MinWormArea = MinWormArea;
expData.TrackStats.MaxWormArea = MaxWormArea;
expData.TrackStats.AnimalPix = AnimalPix;
expData.TrackStats.Level = trackSett.Level;
expData.TrackStats.MinDistance = MinDistance;
expData.TrackStats.SizeChangeThreshold = SizeChangeThreshold;
expData.TrackStats.MinTrackLength = MinTrackLength;
expData.TrackSettings = Settings;

% Save Fragment File
data.movieName = movieName;
data.arenaName = arenaName;
data.wormtrack.tracks = tracks;
data.wormtrack.allData = allData;
data.wormtrack.expData = expData;
data.wormtrack.trackSett = trackSett; %#ok<STRNU>


saveViaAppendData(outputFile,data); % resaves all data, slow! but works with new format
disp([datestr(now),' *** Save complete *** ']);

end