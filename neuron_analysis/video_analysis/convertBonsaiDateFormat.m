clear all;
Table = readtable('3camera_TS.csv');
date = Table.Var17;

time = cellfun(@(fun) fun(12:end-6), date, 'UniformOutput', false); % @(fun) fun...: way to write a function

hours = str2num(cell2mat(cellfun(@(time) time(1:2), time, 'UniformOutput', false)));
minutes = str2num(cell2mat(cellfun(@(time) time(4:5), time, 'UniformOutput', false)));
seconds = str2num(cell2mat(cellfun(@(time) time(7:end), time, 'UniformOutput', false)));

conv_time = hours*60*60 + minutes*60 + seconds;
conv_diff = diff(conv_time);
fps = numel(conv_time) / (conv_time(end) - conv_time(1));

plot(conv_diff)