function wf = waveform_analysis(neuron,varargin)

% This function calculates the features of the waveforms and/ or plots the
% waveforms of the channel where the spike is detected. User input is
% needed when using the 'features' option.
% 
% INPUTS: 
%   - group: .mat file of the group where the neuron is. ex. 'GR1.mat'
%   - neuron: .mat file with thw timestamps of the neuron. ex. 'GR1_1.mat'
%   Varargin
%   - 'plot': When you want to plot the waveform. The default plot is the
%   mean waveform. You can choose to plot only the mean waveform ('mean') or
%   all the waveforms with overlaped mean waveform ('all'). xaxis values
%   15u = 0.5ms.
%   - 'features': generates the structure wf, that contains all the
%   calculated features.
%   - 's': saves the figure in the specified file format (jpg, png, svg,
%   etc.), this format can be changed in the 'type' variable in line 44.
%   See saveas help for more info. Default NOT save. Default format 'jpg'.
%
% OUTPUTS: 
%   - wf: a structure with the calculated features. For more information
%   about these features click to see <a href="matlab:A =
%   imread('L:\Cecilia\WFparams.tif'); imshow(A);">Features image</a>. 
%   - without output arguments, the function plots the mean waveform.
% 
% Examples: 
% wf = waveform_analysis(neuron,'features');
% wf = waveform_analysis(neuron,'features','plot','mean','s','jpg');
% waveform_analysis(neuron); 
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
%
% MATLAB toolboxes: Signal Processing Toolbox.
% -------------------------------------------------------------------------

clearvars -except neuron varargin 

% Defining params ---------------------------------------------------------
if ~isempty(varargin) && max(strcmp(varargin, 's')) == 1
    s = 1;
    type = 'jpg';
else
    s = 0;
end
if max(strcmp(varargin, 'NoType'))
    TypeID = 0;
end
group = neuron(1:end-6);

% Load variables from gourp and neuron ------------------------------------
load(group,'TimeStamps');
load(group,'WaveForms');
load(neuron,'TS');

A = TS/10000;
C = zeros(length(TimeStamps),1);
time = linspace(0, 0.9, size(WaveForms,3)); 
% NOTE: time in ms, 0 to 0.9 defined by the oedisc_probe waveform window
% Waveform length (27 units) = 900us; 15 units = 500us.

% Timestamps relative to the neuron ---------------------------------------
for kk = 1:length(A)
    B = TimeStamps==A(kk,1); 
    C(B) = 1;    
end

% Calculate the mean waveform of the neuron -------------------------------
for n = 1:size(WaveForms,2)
D = WaveForms(:,n,:);
E = D(logical(C),:,:);
F = mean(E(:,:));
F2(n,:) = F; 
end

H = [min(F2,[],2),max(F2,[],2)];
H2 = H(:,2)-H(:,1);
H3 = max(H2);
G = H2 == H3;
D = WaveForms(:,G,:); % Select the wave from G channel
E = D(logical(C),:,:);
allUP = E(:,:); % All wf of the neuron
allDOWN = -allUP; % Inverted values
mUP = mean(allUP); % Mean wf
mDOWN = -mUP; % Inverted values

% Select the type of neuron -----------------------------------------------
if TypeID == 0
else
    neuronType = imread('L:\Cecilia\WFparams_Types.tif');
    figure('name','WaveformType','Position',[100 200 [838, 600]]);
    subplot(2,1,1);
    plot(time,mUP,'b','LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
    title('Your waveform');
    subplot(2,1,2,'Position',[.0 -.25 [1, 1]]);
    imshow(neuronType);
    disp('Please select the waveform type,');
    TypeID = input('[1/ 2]:');
    close WaveformType
end

% Calculate features using the mean waveform, depends on waveform type ----
if max(strcmp(varargin, 'features')) == 1 && max(strcmp(varargin, 'plot') == 1)
    
    % Plot 
    if max(strcmp(varargin, 'mean')) == 1 % Plot mean waveform
        plot(mUP);
        xticks([0 15 30]);
        xticklabels([0 0.5 1]);
        xlabel('msec')
        ylabel('Amplitude')
        title(['Waveform - ', neuron(1,1:end-4),' - CH',...
            num2str((find (G == 1)))],'Interpreter','none');
        set(gca,'fontname','arial');
        
    elseif max(strcmp(varargin, 'all')) == 1 % Plot all waveforms
        for kk = 1:length(E)
        plot(allUP(kk,:),'k');
        hold on;
        end
        hold on;
        plot(mUP,'r');
        hold off;
        xticks([0 15 30]);
        xticklabels([0 0.5 1]);
        xlabel('msec')
        ylabel('Amplitude')
        title(['Waveform - ', neuron(1,1:end-4),' - CH', ...
            num2str((find (G == 1)))],'Interpreter','none');
        set(gca,'fontname','arial');
    end

    % Features 
    pk  = zeros(1,2); % pk or p refer to the height
    loc = zeros(1,2); % loc or l refer to the time
    pr  = zeros(2,2); % pr or r refer to the halfwidth
    switch TypeID
        case 0
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mDOWN,time,'MinPeakHeight',0.01);
        case 1
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mDOWN,time,'MinPeakHeight',0.01);
        case 2
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mDOWN,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mUP,time,'MinPeakHeight',0.01);
            pk2(1,1)  = -pk2(1,1);
            p = -p;
    end
    
    % We want these variables being a 1 row vector.
    if size(p,1) > size(p,2)
        p = transpose(p);
    end
    if size(l,1) > size(l,2)
        l = transpose(l);
    end
    if size(r,1) > size(r,2)
        r = transpose(r);
    end
    
    if size(p,2) == 1 || size(l,2) == 1 || size(r,2) == 1
        pk(:,1)  = p;
        loc(:,1) = l;
        pr(:,1)  = r;
    elseif size(p,2) >= 2
        A = find(l < pk2(2,1));
        B = find(l > pk2(2,1));
        try
            pk(1,1) = p(A(end));
        catch
            pk(1,1) = 0;
        end
        pk(1,2)  = p(B(1));
        try
            loc(1,1) = l(A(end));
        catch
            loc(1,1) = 0;
        end
        loc(1,2) = l(B(1));
        try
            pr(1,1) = r(A(end));
        catch
            pr(1,1) = 0;
        end        
        pr(1,2)  = r(B(1));
    end
    
    % Calculate values
    wf.a = pr(1,1);
    wf.b = pr(1,2);
    wf.pk1(1,1) = -pk(1,1);
    wf.pk1(1,2) = loc(1,1);
    wf.pk2 = transpose(pk2);
    wf.pk3(1,1) = -pk(1,2);
    wf.pk3(1,2) = loc(1,2);
    wf.c = wf.pk2(1,2) - wf.pk1(1,2);
    wf.d = wf.pk3(1,2) - wf.pk2(1,2);
    wf.e = abs(wf.pk2(1,1)) + abs(wf.pk1(1,1));
    wf.f = abs(wf.pk3(1,1)) + abs(wf.pk2(1,1));
    wf.g = w2;
    switch TypeID
        case 0
            F  = mDOWN>0;
            F2 = (time > wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'linear');
            F4 = time(find(F3,2,'first'));
            I(1,2) = F4(2);
            F  = mDOWN>0;
            F2 = (time < wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'variance');
            if max(F3) ~= 1
                F4 = 0;
            else
                F4 = time(find(F3,2,'last'));
            end
            I(1,1) = F4(1);
        case 1
            F  = mDOWN>0;
            F2 = (time > wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'linear');
            F4 = time(find(F3,2,'first'));
            I(1,2) = F4(2);
            F  = mDOWN>0;
            F2 = (time < wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'variance');
            if max(F3) ~= 1
                F4 = 0;
            else
                F4 = time(find(F3,2,'last'));
            end
            I(1,1) = F4(1);
        case 2
            F  = mUP>0;
            F2 = (time > wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'linear');
            F4 = time(find(F3,2,'first'));
            I(1,2) = F4(2);
            F = mUP>0;
            F2 = (time < wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'variance');
            F4 = time(find(F3,2,'last'));
            I(1,1) = F4(1);
    end
    wf.first = wf.pk2(1,2) - I(1,1);
    wf.last  = I(1,2) - wf.pk2(1,2);
    wf.mWF   = mUP;
    wf.time  = time;

% Plots the waveform ------------------------------------------------------
elseif max(strcmp(varargin, 'plot')) == 1
    if max(strcmp(varargin, 'mean')) == 1 % Plot mean waveform
        plot(mUP);
        xticks([0 15 30]);
        xticklabels([0 0.5 1]);
        xlabel('msec')
        ylabel('Amplitude')
        title(['Waveform - ', neuron(1,1:end-4),' - CH',...
            num2str((find (G == 1)))],'Interpreter','none');
        set(gca,'fontname','arial');
    elseif max(strcmp(varargin, 'all')) == 1 % Plot all waveforms
        for kk = 1:length(E)
        plot(allUP(kk,:),'k');
        hold on;
        end
        hold on;
        plot(mUP,'r');
        hold off;
        xticks([0 15 30]);
        xticklabels([0 0.5 1]);
        xlabel('msec')
        ylabel('Amplitude')
        title(['Waveform - ', neuron(1,1:end-4),' - CH',...
            num2str((find (G == 1)))],'Interpreter','none');
        set(gca,'fontname','arial');
    end

% Calculate features using the mean waveform, depends on waveform type ----
elseif max(strcmp(varargin, 'features')) == 1
    
    pk  = zeros(1,2); % pk or p refer to the height
    loc = zeros(1,2); % loc or l refer to the time
    pr  = zeros(2,2); % pr or r refer to the halfwidth
    switch TypeID
        case 0
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mDOWN,time,'MinPeakHeight',0.01);
        case 1
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mUP,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mDOWN,time,'MinPeakHeight',0.01);
        case 2
            [pk2(1,1),pk2(2,1),w2,~] = findpeaks(mDOWN,time,'SortStr',...
                'descend','NPeaks',1);
            [p,l,~,r] = findpeaks(mUP,time,'MinPeakHeight',0.01);
            pk2(1,1)  = -pk2(1,1);
            p = -p;
    end
    
    % We want these variables being a 1 row vector.
    if size(p,1) > size(p,2)
        p = transpose(p);
    end
    if size(l,1) > size(l,2)
        l = transpose(l);
    end
    if size(r,1) > size(r,2)
        r = transpose(r);
    end
    
    if size(p,2) == 1 || size(l,2) == 1 || size(r,2) == 1
        pk(:,1)  = p;
        loc(:,1) = l;
        pr(:,1)  = r;
    elseif size(p,2) >= 2
        A = find(l < pk2(2,1));
        B = find(l > pk2(2,1));
        try
            pk(1,1) = p(A(end));
        catch
            pk(1,1) = 0;
        end
        pk(1,2)  = p(B(1));
        try
            loc(1,1) = l(A(end));
        catch
            loc(1,1) = 0;
        end
        loc(1,2) = l(B(1));
        try
            pr(1,1) = r(A(end));
        catch
            pr(1,1) = 0;
        end        
        pr(1,2)  = r(B(1));
    end
    
    % Calculate values
    wf.a = pr(1,1);
    wf.b = pr(1,2);
    wf.pk1(1,1) = -pk(1,1);
    wf.pk1(1,2) = loc(1,1);
    wf.pk2 = transpose(pk2);
    wf.pk3(1,1) = -pk(1,2);
    wf.pk3(1,2) = loc(1,2);
    wf.c = wf.pk2(1,2) - wf.pk1(1,2);
    wf.d = wf.pk3(1,2) - wf.pk2(1,2);
    wf.e = abs(wf.pk2(1,1)) + abs(wf.pk1(1,1));
    wf.f = abs(wf.pk3(1,1)) + abs(wf.pk2(1,1));
    wf.g = w2;
    switch TypeID
        case 0
            F  = mDOWN>0;
            F2 = (time > wf.pk2(2));
            F3 = ischange(double(F & F2),'linear');
            F4 = time(find(F3,2,'first'));
            try 
                I(1,2) = F4(2);
            catch
                I(1,2) = 0;
            end
            F  = mDOWN>0;
            F2 = (time < wf.pk2(2));
            F3 = ischange(double(F & F2),'variance');
            if max(F3) ~= 1
                F4 = 0;
            else
                F4 = time(find(F3,2,'last'));
            end
            I(1,1) = F4(1);
        case 1
            F  = mDOWN>0;
            F2 = (time > wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'linear');
            F4 = time(find(F3,2,'first'));
            I(1,2) = F4(2);
            F  = mDOWN>0;
            F2 = (time < wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'variance');
            if max(F3) ~= 1
                F4 = 0;
            else
                F4 = time(find(F3,2,'last'));
            end
            I(1,1) = F4(1);
        case 2
            F  = mUP>0;
            F2 = (time > wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'linear');
            F4 = time(find(F3,2,'first'));
            I(1,2) = F4(2);
            F = mUP>0;
            F2 = (time < wf.pk2(2));
            F3 = ischange(double(F & transpose(F2)),'variance');
            F4 = time(find(F3,2,'last'));
            I(1,1) = F4(1);
    end
    wf.first = wf.pk2(1,2) - I(1,1);
    wf.last  = I(1,2) - wf.pk2(1,2);
    wf.mWF   = mUP;
    wf.time  = time;

% Default only mean plot --------------------------------------------------
else % If there are no varargins plot mean waveform
    plot(mUP);
    xticks([0 15 30]);
    xticklabels([0 0.5 1]);
    xlabel('msec')
    ylabel('Amplitude')
    title(['Waveform - ', neuron(1,1:end-4),' - CH',...
        num2str((find (G == 1)))],'Interpreter','none');
    set(gca,'fontname','arial');
end

% Saving ------------------------------------------------------------------
if s == 1
    saveas(gcf, genvarname([neuron(1:end-4),'_',num2str(G),'_wf']), type);
elseif s ~= 1
end    

if s == 1 && max(strcmp(varargin, 'features')) == 1 % Save structure if features calculated
    disp(['Save in ' neuron '?']);
    saveID = input('[Y/ N]:','s');
    switch saveID
        case {'Y', 'y',''}
            save(neuron(1:end-4), 'wf', '-append');
        case {'N', 'n'}
            % Nothing happens
    end
end

clearvars TimeStamps WaveForms TS

end
