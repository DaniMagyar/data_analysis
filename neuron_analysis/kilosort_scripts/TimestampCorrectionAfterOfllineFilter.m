function TimestampCorrectionAfterOfllineFilter

info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');
folder = [cd '\continuous\' info.Header.folder_name];
contFile=fullfile(folder,'continuous.dat');
file=dir(contFile);
samples=file.bytes/2/info.Header.num_channels;
oldTS = readNPY([folder 'timestamps.npy']);
oldSynchTS = readNPY([folder 'synchronized_timestamps.npy']);
timestamps = oldTS(1:samples);
synchronized_timestamps = oldSynchTS(1:samples);
writeNPY(timestamps,[folder 'timestamps.npy']);
writeNPY(synchronized_timestamps,[folder 'synchronized_timestamps.npy']);


