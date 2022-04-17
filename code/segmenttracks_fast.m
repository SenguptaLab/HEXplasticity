%
% Identifies instantaneous behavioral state from worm tracks and
% morphological data from ArenaTracker.
% Saves behavioral data to a file (*_seg.mat)
%
% USAGE:
%   SaveList = SegmentTracks(FileName)
%
%   FileName: single filename of tracked data or cell array of multiple
%               filenames. Select with user input if none given.
%   SaveList: cell array of output filename(s).

%----------------------------
% Dirk Albrecht
% Version 1.0
% 30-Mar-2011 14:17:48
%----------------------------

function SaveList = segmenttracks_fast(fileName)
setts = getComputerSettings;

if nargin < 1
    [fileName, pathName] = uigetfile([setts.searchPath '*.mat'],'Select Track File for Analysis','MultiSelect','off');
    fullName = [pathName filesep fileName];
else
    fullName = fileName;
end

SaveList = {};

% Initialize Segment Analysis Settings
SegmentationSettings;



%%
[pathname,filename,ext] = fileparts(fullName);  
SaveName = [filename ext];
FullSaveName = fullfile(pathname, SaveName);
%behSaveName = [FullSaveName(1:end-8),'_beh.mat'];

%% ------------
% Load Data
%------------

%folderdate = []; founddate = strfind(fullName,'200'); if length(founddate)>0 folderdate = fullName(founddate+(0:7)); end
%disp([datestr(now),': Loading File #',num2str(fnum),': [',folderdate,'] ',filename]);

data = load(fullName,'wormtrack');
% data = data.data;
tracks = data.wormtrack.tracks;
trackSett = data.wormtrack.trackSett;
expData = data.wormtrack.expData;
allData = data.wormtrack.allData;

%%

disp(['...loaded ',datestr(now)]);

%allData = single(allData);

if exist('expData','var') && ~isempty(expData)
    framerate = expData.FrameRate;
    pixelsize = expData.PixelSize;
else
    framerate = 2;
    pixelsize = 1;
    disp('WARNING: NO SCALING DATA');
end
% now overwrite the defaults!!!
Settings.FrameRate = framerate;
Settings.PixelSize = pixelsize;
Settings.XBorderPadding = 0.5*pixelsize;
Settings.StallDistance = 0.022;

Settings.MaxPathAngleDev = 25;
Settings.MaxBodyPathAngleDev = 25;
Settings.MinSpeedFor = 1.2;
Settings.MinFwdRunFr = 2;
Settings.CollEffectTime = 2;

if isfield(tracks,'Segment')
    ButtonName=questdlg('Segments have already been analyzed.','','Reanalyze', 'Stop', 'Reanalyze');
    switch ButtonName
        case 'Stop'
            disp('ending.'),
            return;
    end
end

TrackArena = ones(1, length(tracks));

%---------------------------
numtracks = length(tracks);
trackstats = [];
totalfr = 0;
sprg = -1;
    gSigma = 0.5;
    g = fspecial('gauss',round((gSigma*2)*2+1),gSigma);
    g = g(round((gSigma*2)*2+1),:)./sum(g(round((gSigma*2)*2+1),:));

%%
tic;

for tr = 1:length(tracks)
    %% --------------------------------------------------------------------------
    %  LOOP FOR EACH TRACK
    %--------------------------------------------------------------------------
    newtracks = tracks(tr);
    newtracks.Path(:,1) = imfilter(newtracks.Path(:,1),g','replicate','conv','same')';
    newtracks.Path(:,2) = imfilter(newtracks.Path(:,2),g','replicate','conv','same')';

    [Segment, SegAnalysis] = SegmentTrack(newtracks,Settings);
    trackbox = [1 1 expData.ArenaSize];   
    TrackAnalysis = AnalyzeTrack(Segment,SegAnalysis,trackbox,Settings);
    %%
    reversal = (TrackAnalysis.FRP == 2);
    
    %------------------------------
    % reorganize data for plotting
    %------------------------------
    % segment codes:
    %   1 - fwd
    %   2 - fwd w/ revs
    %   3 - lane chg
    %   4 - small curve (60deg)
    %   5 - short reversal
    %   6 - long curve / loopy movement
    %   7 - pirouette (sharp curve following reversal)
    %   8 - unknown/erratic movement
    %   9 - pause
    %  10 - out of bounds / collision
    
    distance = []; pathang = []; pathangvel = [];
    for seg = 1:length(Segment)
        distance(Segment(seg).RealFrames) = Segment(seg).Distance';
        pathang(Segment(seg).RealFrames) = Segment(seg).PathAngle';
        pathangvel(Segment(seg).RealFrames) = Segment(seg).PathAngVel';
    end
    
    distance(SegAnalysis.StallFr) = 0;
    pathang(SegAnalysis.StallFr) = NaN;
    pathangvel(SegAnalysis.StallFr) = NaN;
    
    spd = distance / expData.PixelSize * expData.FrameRate .* (1-2*reversal);
    
    t = toc;
    trlength = length(TrackAnalysis.Code);
    numsegments = length(Segment);
    totalfr = totalfr + trlength;
    
    status = sprintf('Track: %d/%d [%d fr] %d @ %d fps]',tr,numtracks,trlength,totalfr,round(totalfr/t));
    sprg = showprog(status,sprg);
    
    trackstats = [trackstats; trlength, numsegments];
    
    tracks(tr).Code = TrackAnalysis.Code;
    tracks(tr).Distance = distance;
    tracks(tr).Speed = spd;
    tracks(tr).PathAngle = pathang;
    tracks(tr).PathAngVel = pathangvel;
    tracks(tr).Reverse = reversal;
    tracks(tr).OriginalDistance = SegAnalysis.OriginalDistance;
    tracks(tr).Segment = Segment;
    tracks(tr).Stall = SegAnalysis.StallFr;
    tracks(tr).NoStall = SegAnalysis.NoStallFr;
    tracks(tr).Omega = SegAnalysis.Omega;
    tracks(tr).FTurnCCW = SegAnalysis.FTurnCCW;
    if size(tracks(tr).FTurnCCW,1) > 0 tracks(tr).FTurnCCW = tracks(tr).FTurnCCW(find(TrackAnalysis.Code(tracks(tr).FTurnCCW(:,1)) ~= 10),:); end
    tracks(tr).FTurnCW = SegAnalysis.FTurnCW;
    if size(tracks(tr).FTurnCW,1) > 0 tracks(tr).FTurnCW = tracks(tr).FTurnCW(find(TrackAnalysis.Code(tracks(tr).FTurnCW(:,1)) ~= 10),:); end
    
    tracks(tr).TurnOmega = TrackAnalysis.TurnOmega;
    tracks(tr).Beh = TrackAnalysis.Beh;
end

if ~isfield(expData,'FrameRate') expData.FrameRate = 2; end

Arena = 1;
expData.Arena = Arena;

%-----------------------
% Link Tracks together
%-----------------------
linkoutput = LinkTracks(tracks,0.5*expData.PixelSize/expData.FrameRate,Inf);  %Max avg velocity = 0.5mm/s
[tracks.OriginalTrack] = deal(linkoutput.OriginalTrack);
expData.Animals = max(struct2mat(1,linkoutput,[],{'OriginalTrack'}));


segmentation.settings = Settings;
segmentation.tracks = tracks;
segmentation.expData = expData;
segmentation.allData = allData;
segmentation.trackstats = trackstats;
segmentation.link = linkoutput;

save(FullSaveName,'segmentation', '-append');

%Save condensed data
disp('Condensing behavior data...');
% 
% %behSaveName = [FullSaveName(1:end-8),'_beh.mat'];
% tfields = fieldnames(tracks);
% rfields = tfields(find(~strcmp(tfields,'Beh') & ...
%     ~strcmp(tfields,'Frames') & ...
%     ~strcmp(tfields,'OriginalTrack') & ...
%     ~strcmp(tfields,'PathAngle') & ...
%     ~strcmp(tfields,'Speed')));
% for tr=1:length(tracks)
%     tracks(tr).X = tracks(tr).Path(:,1)';
%     tracks(tr).Y = tracks(tr).Path(:,2)';
% end
% tracks = rmfield(tracks,rfields);
% tracks = singleStruct(tracks);

%save(behSaveName,'Tracks');
disp(['done. ',datestr(now)]);

%SaveList = cat(1,SaveList,{FullSaveName});

end


