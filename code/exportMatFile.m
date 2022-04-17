[file, path] = uigetfile
exampleObject = matfile([path,file]);
dataset = exampleObject.ethogram;
info = exampleObject.ExpInfo;
preprocess = exampleObject.preprocess;
pixelSize = preprocess.pixelSize;
bgLuminance = preprocess.bgLuminaceDistribution.data{1,1};
dlmwrite([path,'behmat.csv'], dataset.behmat);
dlmwrite([path,'dirmat.csv'], dataset.dirmat);
dlmwrite([path,'spdmat.csv'], dataset.spdmat);
dlmwrite([path,'ymat.csv'], dataset.ymat);
dlmwrite([path,'xmat.csv'], dataset.xmat);
writetable(struct2table(info),[path,'Expinfo.csv']);
writetable(struct2table(preprocess),[path, 'preprocess.csv']);
csvwrite([path, 'luminance.csv'], bgLuminance);