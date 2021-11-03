%% add user input for date, genotype, condition, exp #, stimulus

function preprocessMovie
%%
global hconcTop hconcBottom hconcEdge hconcCenter hConcMode


setts = getComputerSettings;
maxY = setts.maxY;

UserPrompt = ('Select AVI File For Analysis:');
[movieName, pathName] = uigetfile([setts.searchPath '*.avi'], UserPrompt);
MovieFile = [pathName movieName];
mov = VideoReader(MovieFile);
firstFrame = read(mov, 2);
%%
prompt = {'date','genotype','condition','experiment number','stimulus','orient'};
dlg_title = 'Enter User Input';
num_lines = 1;
def = {'20150131','N2','sparse','1','hex','0'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
date = answer{1}; %yyyymmdd
genotype = answer{2};
condition = answer{3};
experiment = answer{4};
stimulus = answer{5};
grad_orient = answer{6}; %binary ans 0 = [high] at bottom (y = 16mm), 1 = [high] at top
ExpInfo.date = date;
ExpInfo.genotype = genotype;
ExpInfo.condition = condition;
ExpInfo.experimentNum = experiment;
ExpInfo.stimulus = stimulus;
ExpInfo.gradOrient = grad_orient;

%%
% defaults:
nrArenas = 2;
arenaSize = 16.1;
defaultStartFrame = 2;

% show first frame
close all
WTFigH = findobj('Tag', ['WTFIG_' movieName]);
if isempty(WTFigH)
    WTFigH = figure('NumberTitle', 'off', 'Tag', ['WTFIG_' movieName]);clf
else
    figure(WTFigH);clf
end
set(gca,'position',[0 0 1 1],'units','normalized')

movSize = size(firstFrame);
imshow(firstFrame);

set(WTFigH,'MenuBar','none');
set(WTFigH,'Toolbar','none');
imageScale = 1;
if movSize(2) >maxY+100
    imageScale = (maxY+100)./movSize(2);
end

set(WTFigH,'Position',[100 60 movSize(2)*imageScale movSize(1)*imageScale]);
axis off
%
% --------get information----------
buttony = 10; buttonyWidth = 17;
nrFrames = mov.NumberOfFrames;
FrameRate = mov.FrameRate;

uicontrol('Style','text','Position',[10 buttony 40 buttonyWidth],'String','Frames');
hStartFrame = uicontrol('style','edit','Position',[60 buttony 50 buttonyWidth],'String',num2str(defaultStartFrame));
hStopFrame = uicontrol('style','edit','Position',[110 buttony 50 buttonyWidth],'String',num2str(nrFrames));
uicontrol('Style','text','Position',[170 buttony 60 buttonyWidth],'String','NrArenas');
hNrArenas = uicontrol('style','edit','Position',[220 buttony 20 buttonyWidth],'String',num2str(nrArenas));

uicontrol('Style','text','Position',[260 buttony 60 buttonyWidth],'String','Arena(mm)');
hArena = uicontrol('style','edit','Position',[320 buttony 40 20],'String',num2str(arenaSize));

uicontrol('Style','text','Position',[370 buttony 60 buttonyWidth],'String','FrameRate');
hFrameRate = uicontrol('style','edit','Position',[430 buttony 30 20],'String',num2str(FrameRate));

uicontrol('Style','text','Position',[470 buttony 80 buttonyWidth],'String','background Fr.');
hBG = uicontrol('style','edit','Position',[550 buttony 30 buttonyWidth],'String',num2str(600));
%%
uicontrol('Style','text','Position',[590 buttony 60 buttonyWidth],'String','Conc Mode');
hConcMode = uicontrol('Style','popupmenu','Position',[660 buttony 70 buttonyWidth],'String',{'linear','quadratic'},'Callback', @setconcmode_helper);
setconcmode_helper


%%

hStartArenaDef = uicontrol('Style','togglebutton','Position',[800 buttony+40 50 buttonyWidth+60],'String','start');
%%
arenaNames = struct;
arenaNames.hstat = uicontrol('Style','text','Position',[980 buttony 60 buttonyWidth],'String','arena 1 name:');
arenaNames.hedit = uicontrol('style','edit','Position',[1040 buttony 80 20],'String','arena_1');
arenaNames.end   = 1150;
arenaNames.active = 1;
%
while ~get(hStartArenaDef,'Value')
    
    if arenaNames(1).active < str2double(get(hNrArenas,'String'))
        arenaNames(1).active = arenaNames(1).active+1;
        arenaNames(arenaNames(1).active).hstat = uicontrol('Style','text','Position',[arenaNames(arenaNames(1).active-1).end+10 buttony 60 buttonyWidth],'String',['arena ' num2str(arenaNames(1).active) ' name:']);
        arenaNames(arenaNames(1).active).hedit = uicontrol('style','edit','Position',[arenaNames(arenaNames(1).active-1).end+10+70 buttony 80 20],'String',['arena_' num2str(arenaNames(1).active)] );
        arenaNames(arenaNames(1).active).end = arenaNames(arenaNames(1).active-1).end+160;
    elseif  arenaNames(1).active > str2double(get(hNrArenas,'String'))
        set(arenaNames(arenaNames(1).active).hstat,'Visible','off');
        set(arenaNames(arenaNames(1).active).hedit,'Visible','off');
        arenaNames(1).active = arenaNames(1).active-1;
    end
    drawnow;
end
%
preprocess.userInput.nrArenas = str2double(get(hNrArenas,'String'));
preprocess.userInput.startFrame = str2double(get(hStartFrame,'String'));
preprocess.userInput.stopFrame = str2double(get(hStopFrame,'String'));
preprocess.userInput.framerate  = str2double(get(hFrameRate,'String'));
preprocess.userInput.arenaSize  = str2double(get(hArena,'String'));
preprocess.userInput.arenaName = cell(1,arenaNames(1).active);
preprocess.userInput.bgFrames = str2double(get(hBG,'String'));
concmode = get(hConcMode,'String');
concval  = get(hConcMode,'Value');
preprocess.userInput.concMode = concmode{concval};
switch concmode{concval}
    case 'linear'
        preprocess.userInput.concTop = str2double(get(hconcTop,'String'));
        preprocess.userInput.concBottom = str2double(get(hconcBottom,'String'));
        preprocess.userInput.concMax = max([preprocess.userInput.concTop preprocess.userInput.concBottom]);
        preprocess.userInput.concMin = min([preprocess.userInput.concTop preprocess.userInput.concBottom]);
    case 'quadratic'
        preprocess.userInput.concEdge = str2double(get(hconcEdge,'String'));
        preprocess.userInput.concCenter = str2double(get(hconcCenter,'String'));
        preprocess.userInput.concMax = max([preprocess.userInput.concEdge  preprocess.userInput.concCenter]);
        preprocess.userInput.concMin = min([preprocess.userInput.concEdge  preprocess.userInput.concCenter]);
end


for ar = 1:arenaNames(1).active
    preprocess.userInput.arenaName{ar} = get(arenaNames(ar).hedit,'String');
end

% get scale info
txt = ['GET SCALING: pick points on top and bottom arena edges '  ' = ' num2str(preprocess.userInput.arenaSize) 'mm'];
label = text(movSize(2)/2+20,movSize(1)*0.4+40,txt,'FontSize',18,'HorizontalAlignment','center');

title(txt); set(label,'String',txt);
[~,Y] = ginput(2);
scalepix = abs(Y(2)-Y(1));
set(label,'String','');

txt = ['SELECT TRACKING AREA (start with #1): click inside to confirm, outside to redo.'];
label = text(movSize(2)/2,movSize(1)*0.4,txt,'FontSize',18,'HorizontalAlignment','center');
title(txt); set(label,'String',txt);
for ar = 1:preprocess.userInput.nrArenas
    success = 0;
    while success == 0
        box = getrect(1); box = box+(box==0); h = rectangle('Position',box); set(h,'EdgeColor',[1,0,0]);
        [X,Y,button] = ginput(1);
        if button == 1 && (X-box(1) >= 0 && X-box(1) <= box(3) && Y-box(2) >= 0 && Y-box(2) <= box(4))
            preprocess.arena(ar).trackBox = box+(box==0);
            set(h,'EdgeColor',[0,0,1]);
            success = true;
        else
            set(h,'Visible','off');
        end
    end
end
set(label,'String','');
%%
txt = ['SELECT BACKGROUND AREA (start with #1): click inside to confirm, outside to redo.'];
label = text(movSize(2)/2,movSize(1)*0.4,txt,'FontSize',18,'HorizontalAlignment','center');
title(txt); set(label,'String',txt);
for ar = 1:preprocess.userInput.nrArenas
    success = 0;
    while success == 0
        box = getrect(1); box = box+(box==0); h = rectangle('Position',box); set(h,'EdgeColor',[0,1,0]);
        [X,Y,button] = ginput(1);
        if button == 1 && (X-box(1) >= 0 && X-box(1) <= box(3) && Y-box(2) >= 0 && Y-box(2) <= box(4))
            backGround(ar).trackBox = box+(box==0);
            set(h,'EdgeColor',[0,0.5,1]);
            success = true;
        else
            set(h,'Visible','off');
        end
    end
    preprocess.background(ar).trackBox = backGround(ar).trackBox;
    preprocess.background(ar).bgRecalculated = true;
end
set(label,'String','');
drawnow
%%
preprocess.NumArenas = ar;
%%
preprocess.pixelSize = scalepix / preprocess.userInput.arenaSize;
preprocess.movSize = movSize;
preprocess.originalFilename = movieName;
preprocess.originalPath = pathName;


%% create movie handles and prepare to write into them!
[~,newName] = fileparts(MovieFile);
ifcomp = strfind(newName,'_comp');
if ~isempty(ifcomp)
    newName(ifcomp:ifcomp+4) = '';
end
%%
for ar = 1:preprocess.userInput.nrArenas
    
    preprocess.outputPaths{ar} =  [pathName   preprocess.userInput.arenaName{ar}];
    if ~exist(preprocess.outputPaths{ar},'dir')
        mkdir(preprocess.outputPaths{ar});
    end
    
    preprocess.wormMovieName{ar} = preprocess.userInput.arenaName{ar};
    
    if exist([preprocess.outputPaths{ar} filesep preprocess.wormMovieName{ar} '_worms.avi'],'file')
        button = questdlg(['Video Exists: ' preprocess.wormMovieName{ar} '_worms.avi, overwrite?'],'Video Writing','yes','no','no');
        switch button
            case 'yes'
                display('video will be overwritten')
            case 'no'
                error('user cancelled storing of video')
        end
    end
    
    
    movWriteWorm{ar} = VideoWriter([preprocess.outputPaths{ar} filesep preprocess.wormMovieName{ar} '_worms'],'Motion JPEG AVI'); %#ok<*AGROW,TNMLP>
    set(movWriteWorm{ar},'FrameRate',preprocess.userInput.framerate);
    set(movWriteWorm{ar},'Quality',85);
    open(movWriteWorm{ar});
    
    
    preprocess.bgMovieName{ar} = [preprocess.userInput.arenaName{ar}];
    
    if exist([preprocess.outputPaths{ar} filesep preprocess.bgMovieName{ar} '_bg.avi'],'file')
        
        button = questdlg(['Video Exists: ' preprocess.bgMovieName{ar} '_bg.avi, overwrite?'],'Video Writing','yes','no','no');
        switch button
            case 'yes'
                display('video will be overwritten')
            case 'no'
                error('user cancelled storing of video')
        end
    end
    
    movWritebg{ar} = VideoWriter([preprocess.outputPaths{ar} filesep preprocess.bgMovieName{ar} '_bg'],'Motion JPEG AVI'); %#ok<TNMLP>
    set(movWritebg{ar},'FrameRate',preprocess.userInput.framerate);
    set(movWritebg{ar},'Quality',75);
    open(movWritebg{ar});
    
    
end




%%
trackBox = round(  cat(1, preprocess.arena(:).trackBox));
if any(trackBox(:,3)> movSize(2))
    decreaseSize = find(trackBox(:,3)> movSize(2));
    for ar = 1:length(decreaseSize)
        trackBox(ar,3) =  movSize(2)-trackBox(ar,1);
    end
end
if any(trackBox(:,4)> movSize(1))
    decreaseSize = find(trackBox(:,4)> movSize(1));
    for ar = 1:length(decreaseSize)
        trackBox(ar,4) =  movSize(1)-trackBox(ar,2);
    end
end

framesToRead = preprocess.userInput.startFrame:preprocess.userInput.stopFrame;
totalFramesToRead = length(framesToRead);
for ar = 1:preprocess.userInput.nrArenas % for plotting later
    wormDataYT{ar} = zeros([ trackBox(ar,4)+1 totalFramesToRead],'single');
    bgDataYT{ar}   = zeros([ trackBox(ar,4)+1 totalFramesToRead],'single');
    bgUseX{ar} = false( 1,trackBox(ar,3)+1);
    
    bgBox =     preprocess.background(ar).trackBox(1,:);
    bgStartPix = bgBox(1)-trackBox(ar,1);
    bgEndPix   = bgBox(1)-trackBox(ar,1)+bgBox(3);
    bgUseX{ar}(round(bgStartPix :bgEndPix))  = true;
end

%%
for frags = 1:ceil(totalFramesToRead./preprocess.userInput.bgFrames)
    
    movieSection = read(mov,[framesToRead((frags-1)*preprocess.userInput.bgFrames+1) framesToRead(min([(frags)*preprocess.userInput.bgFrames totalFramesToRead]))]);
    for ar = 1:preprocess.userInput.nrArenas
        sectionData = squeeze(movieSection(trackBox(ar,2):trackBox(ar,2)+trackBox(ar,4),trackBox(ar,1):trackBox(ar,1)+trackBox(ar,3),1,:));
        newd = single(sectionData);
        meannewd = mean(newd,3);
        newd = bsxfun(@minus,newd,meannewd);
        newd(newd>0) = 0;
        newd = -newd;
        
        if frags == 1
            scalingFactor(ar) = 220./max(newd(:)); %#ok<AGROW>
        end
        writeData = sectionData + uint8(newd);
        writeData  = reshape( writeData ,[size( writeData,1) size( writeData ,2) 1 size(writeData,3)]);
        writeVideo(movWritebg{ar}, writeData);
        
        bgDataYT{ar}(:,(frags-1)*preprocess.userInput.bgFrames+1:min([(frags)*preprocess.userInput.bgFrames totalFramesToRead])) =  mean(single(writeData(:, bgUseX{ar},:,:)),2);
        
        
        writeData = uint8(255-(newd.*scalingFactor(ar)));
        writeData  = reshape( writeData ,[size( writeData,1) size( writeData ,2) 1 size(writeData,3)]);
        writeVideo(movWriteWorm{ar}, writeData );
        wormDataYT{ar}(:,(frags-1)*preprocess.userInput.bgFrames+1:min([(frags)*preprocess.userInput.bgFrames totalFramesToRead])) =  mean(single(writeData),2);
        
    end
    clear movieSection
end
preprocess.scalingFactor = scalingFactor;
clear movieSection newd meannewd movieSection writeData
%%

for ar = 1:preprocess.userInput.nrArenas
    close(movWriteWorm{ar});
    close(movWritebg{ar});
end

%% now the plotting
for ar = 1:preprocess.userInput.nrArenas
    newFig{ar} = figure;clf
    set(newFig{ar},'PaperOrientation','portrait','PaperSize',[8.5 11],'PaperPositionMode','manual','PaperPosition',[0 0 8.5 11]);
    subplot(311)
    meanData = wormDataYT{ar};
    meanData = (meanData-min(meanData(:)));
    meanData = meanData ./max(meanData(:));
    imagesc([0 totalFramesToRead./preprocess.userInput.framerate./60],[0 trackBox(ar,3)]./preprocess.pixelSize, meanData,[0.5 1] );
    preprocess.wormLuminaceDistribution.data{ar} = meanData;
    preprocess.wormLuminaceDistribution.vertPix{ar}  = (0:trackBox(ar,4))./preprocess.pixelSize;
    preprocess.wormLuminaceDistribution.minutes{ar}  = [0 totalFramesToRead./preprocess.userInput.framerate./60];
    
    
    xlabel('time (min)')
    ylabel('vertical position (mm)')
    colormap('gray')
    title('worm body distribution (luminace)')
    subplot(312)
    meanData = bgDataYT{ar};
    meanData = (meanData-min(meanData(:)));
    meanData = meanData ./max(meanData(:));
    imagesc([0 totalFramesToRead./preprocess.userInput.framerate./60],[0 trackBox(ar,3)]./preprocess.pixelSize, meanData,[0 1] );
    xlabel('time (min)')
    ylabel('vertical position (mm)')
    colormap('gray')
    title('background (luminace)')
    preprocess.bgLuminaceDistribution.data{ar} = meanData;
    preprocess.bgLuminaceDistribution.vertPix{ar}  = (0:trackBox(ar,4))./preprocess.pixelSize;
    preprocess.bgLuminaceDistribution.minutes{ar}  = [0 totalFramesToRead./preprocess.userInput.framerate./60];
    
    subplot(313)
    plot((0:trackBox(ar,4))./preprocess.pixelSize, mean(meanData,2))
    xlabel('vertical position (mm)')
    ylabel('luminance')
    axis tight;
    
    suptitle([date ':    ' newName '    -    ' preprocess.userInput.arenaName{ar}]);%,'interpreter','none')
    
    saveas(newFig{ar},[preprocess.outputPaths{ar} filesep newName '_preprocessing.pdf']);
    
    data.ExpInfo = ExpInfo;
    data.preprocess = preprocess;
    
    outputFile = [ preprocess.outputPaths{ar} filesep preprocess.userInput.arenaName{ar} '.mat']; %#ok<NASGU>
    saveViaAppendData(outputFile,data); % resaves all data, slow! but works with new format
    disp([datestr(now),' *** Save complete *** ']);
end
