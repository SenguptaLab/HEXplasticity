%%Script acquires video at 2fps

vid = videoinput('pixelinkimaq', 1, 'MONO8_1920x1200');

src = getselectedsource(vid);
src.Exposure = 3;
src.Gain = '0.00';
src.ActualFrameRate = 4;
src.FrameRate = 4;
vid.ReturnedColorspace = 'grayscale';

prompt = {'Enter Videolength(frames):','Enter Filename:', 'Enter Framerate (s-1)'};
dlg_title = 'Video Specifications';
num_lines = 1;
def = {'2401','','4'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
filename = answer{2};
src.ActualFrameRate = str2double(answer{3});
src.FrameRate = str2double(answer{3});


vid.FramesPerTrigger = str2double(answer{1}); %3600 = 30min at 2fps
% vid.FrameGrabInterval = 7;
vid.LoggingMode = 'disk';
folder_name = uigetdir();
diskLogger = VideoWriter([folder_name filesep filename], 'Motion JPEG AVI');


diskLogger.FrameRate = src.ActualFrameRate; 
vid.DiskLogger = diskLogger;
preview(vid);

button = questdlg('Click OK to start video!','Video Start Box','OK','cancel','cancel');
switch button
    case 'OK'
        run = true;
        display('Recording started.')
        start(vid);
            case 'cancel'
        run = false;
        display('Recording canceled.')
end
%%