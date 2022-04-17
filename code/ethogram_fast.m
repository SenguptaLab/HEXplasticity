%
% Displays ethogram of behavior and summarizes behavioral state probability
% and speed over time.  Adjusts timing according to flow properties, either
% via user input or automatically from dye experiments. Saves data to a
% file (*_ethogram.mat) and prints PDF summary pages.
%
% USAGE:
%   Ethogram(FileName)
%
%   FileName: single filename of segmented data or cell array of multiple
%               filenames. Select with user input if none given.

%----------------------------
% Dirk Albrecht
% Version 1.0
% 31-Mar-2011 12:09:13
%----------------------------

function ethogram_fast(fileName)
setts = getComputerSettings;
%%

if nargin < 1
    % Get track data for analysis
    % --------------------------
    [fileName, pathName] = uigetfile('*.mat','Select segmented track File(s) For Analysis','MultiSelect','on');
    fullName = [pathName filesep fileName];
else
    fullName = fileName;
end
%%
[pathname,filename,ext] = fileparts(fullName);
SaveName = [filename ext];
FullSaveName = fullfile(pathname, SaveName);


data = load(fullName,'segmentation');


Settings = data.segmentation.settings;
Tracks = data.segmentation.tracks;
ExpData = data.segmentation.expData;
allData = data.segmentation.allData;
trackstats = data.segmentation.trackstats;
linkoutput = data.segmentation.link;
%%
    stimulusdata = [];
    ontime = NaN;
    flowinfo.StartDelayFr = NaN;
    flowinfo.TimeDelayFr = NaN;
    flowinfo.UpDownDelayFr = NaN;
    flowinfo.VelocityPixPerFr = NaN;
    flowinfo.XUpDownPix = [1 ExpData.ArenaSize(1)];

%flowinfo = DyeAnalysis(DyeData,ontime*60*ExpData.FrameRate);
ExpData.Flow = flowinfo(1);
ExpData.Flow.VelocityMMperSec = ExpData.Flow.VelocityPixPerFr * Settings.FrameRate / Settings.PixelSize;
ExpData.Flow.StartDelaySec = ExpData.Flow.StartDelayFr / Settings.FrameRate;

% use velocity and start delay data from experiment info file
if isfield(stimulusdata,'FlowVel')
    ExpData.Flow.VelocityMMperSec = str2num(char(stimulusdata.FlowVel));
    ExpData.Flow.VelocityPixPerFr = ExpData.Flow.VelocityMMperSec * Settings.PixelSize / Settings.FrameRate;
end
if isfield(stimulusdata,'Delay')
    ExpData.Flow.StartDelaySec = str2num(char(stimulusdata.Delay));
    ExpData.Flow.StartDelayFr = ExpData.Flow.StartDelaySec * Settings.FrameRate;
end

%--------------------------------------
% Generate ethogram and data matrices
%--------------------------------------

ExpData.Flow.DelayFr = ExpData.Flow.StartDelayFr(1) - ExpData.Flow.XUpDownPix(1) / ExpData.Flow.VelocityPixPerFr(1);

ethfig = findobj(get(0,'Children'),'Tag','Ethogram');
if isempty(ethfig) ethfig = figure; set(ethfig,'Tag','Ethogram'); end

Data = Tracks2Matrix(Tracks,ExpData.Flow,~isnan(ExpData.Flow.DelayFr));

%------------------------------
% Save Data and summary plots
%------------------------------
%%
ethogram = Data;
ethogram.ExpData = ExpData;
save(FullSaveName,'ethogram','-append');

% page 1: Raw matrices
orient(ethfig,'landscape');
saveas(ethfig,strrepl(FullSaveName,'.mat','_raw.pdf'),'pdf');
% Behavior code: 1-F, 2-LF, 3-R, 4-P, 5-OmR, 6-Omf, 7-?, 8-OOB
%cmap = [1 1 1;.7 .7 .7; .7 .7 .7; 0 0 0; .3 .3 .3; .6 0 0; 1 .2 .2; 1 1 1; .9 .9 .9];
cmap = [.7 .7 .7; .7 .7 .7; 0 0 0;.3 .3 .3;.6 0 0; 1 .2 .2; 1 1 1; 1 1 1; .9 .9 .9];
   
clf;
t = (1:size(Data.behmat,2)) / 60 / ExpData.FrameRate; % min
subplot(5,1,1); image(t,1:size(Data.behmat,1),ind2rgb(Data.behmat,cmap)); ylabel('Animal #');
if isnan(ExpData.Flow.VelocityMMperSec)
    FlowLabel = '(real time)';
else
    FlowLabel = sprintf('(flow: %1.3f mm/s, delay: %1.3f s)',ExpData.Flow.VelocityMMperSec(1),ExpData.Flow.StartDelaySec(1));
end

%%
title(sprintf('%s %s',SaveName,FlowLabel),'interpreter','none');
subplot(5,1,2); stateplot(Data.behprob,[],t,0,0,0); ylabel('State probability');
subplot(5,1,[3:4]); stateplot(Data.behprob,[],t,.05,0,0); hilite(ontime,[],[1 1 .5]); ylabel('State probability'); xlim([0 t(end)]);
subplot(5,1,5); plot(t,Data.speed.fwdpause); ylim([0 0.4]); hilite(ontime,[],[1 1 .5]); ylabel('Speed (mm/s)'); xlim([0 t(end)]);
xlabel('Time (min)');
orient(ethfig,'tall');
saveas(ethfig,strrepl(FullSaveName,'.mat','_avg.pdf'),'pdf');

% page 3: Population state probability and speed, averaged over
% repeated cycles
if exist('numcycles') && numcycles > 1
    clf;
    cycbehmat = split(Data.behmat,4,2,1);
    t = (1:size(cycbehmat,2)) / 60 / ExpData.FrameRate; % min
    subplot(5,1,1); image(t,1:size(cycbehmat,1),ind2rgb(cycbehmat,cmap)); ylabel('Animal #');
    title(sprintf('%s %s',SaveName,FlowLabel),'interpreter','none');
    cycbehhist = hist(cycbehmat,1:8);
    cycbehnum = sum(cycbehhist(1:6,:));
    cycbehprob = cycbehhist(1:6,:) ./ repmat(cycbehnum,6,1);
    subplot(5,1,2); stateplot(cycbehprob,[],t,0,0,0); ylabel('State probability');
    subplot(5,1,[3:4]); stateplot(cycbehprob,[],t,.05,0,0); hilite(ontime,[],[1 1 .5]); ylabel('State probability');  xlim([0 t(end)]);
    subplot(5,1,5); plot(t,nanmean(reshape(Data.speed.fwdpause',[],numcycles)')); ylim([0 0.4]); hilite(ontime,[],[1 1 .5]); ylabel('Speed (mm/s)'); xlim([0 t(end)]);
    xlabel('Time (min)');
    orient(ethfig,'tall');
    saveas(ethfig,strrepl(FullSaveName,'.mat','_cycavg.pdf'),'pdf');
end
