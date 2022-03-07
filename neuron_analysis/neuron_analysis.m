function neuron_analysis(neuron, varargin)

% This function calculates several data from a neuron: features of the
% waveforms, plots the waveforms, ISI histogram, Cross correlation and
% Autocorrelation. In the VARARGIN you can coose whitch things to
% calculate, the default mode calculates everything.
% 
% INPUTS: 
%   - neuron: neuron to analyse.
%   Varargin
%   - 'plot': when you want to plot the waveform. The default plot is the
%   mean waveform. You can choose to plot only the mean waveform ('mean') or
%   all the waveforms with overlaped mean waveform ('all').
%   - 'wf': generates the structure wf, that contains all the
%   calculated features.
%   - 'ISI': calculates the interspike interval.
%   - 'XCorr': calculates the cross correlation.
%   - 'ACorr': calculates the auto-correlation.
%   - 's': saves the figure in the specified file format (jpg, png, svg,
%   etc.), this format can be changed in the 'type' variable in line 53.
%   See saveas help for more info. Default NOT save. Default format 'jpg'.
%
% OUTPUTS: 
%   - the calculated variables are automatically saved inside the neuron
%   variable.
%   
% Examples: 
% neuron_analysis(neuron);
% neuron_analysis(neuron,'wf','plot','mean','ISI');
% 
% -------------------------------------------------------------------------
% Cecília Pardo-Bellver, 2019
% Laboratory of Network Neurophysiology
% Institute of Experimental Medicine, Hungary.
%
% Uses code from MClust.
%
% MATLAB toolboxes: Signal Processing Toolbox.
% -------------------------------------------------------------------------

% Default params
ToCalculate = {'wf', 'plot', 'mean', 'ISI', 'XCorr', 'ACorr'};

% Calculations defined in varargin
if ~isempty(varargin)
    ToCalculate = varargin;
else
end

s = 0;
if max(strcmp(ToCalculate, 's')) == 1 % Save the plots
    s = 1;
    type = 'jpg';
elseif max(strcmp(ToCalculate, 's')) ~= 1 % Not save the plots
    s = 0;
end

group = [neuron(1:end-6) '.mat'];
load(group,'TimeStamps','WaveForms'); % Load 
load(neuron,'TS');

% Params, read from varargin
for ii = 1:size(ToCalculate,2)
    switch ToCalculate{ii}
      case 'wf' % Waveform calculations -----------------------------------
        wf = waveform_analysis(group,neuron,'features');
                
      case 'plot' % Plot waveforms ----------------------------------------
        if max(strcmp(ToCalculate, 'mean')) == 1 % Plot mean waveform
            waveform_analysis(group,neuron,'plot','mean');
        elseif max(strcmp(ToCalculate, 'all')) == 1 % Plot all waveforms
            waveform_analysis(group,neuron,'plot','all');
        end
        
        % Extract the waveforms corresponding to the TS (neuron)
        A = TS/10000;
        C = zeros(length(TimeStamps),1);

            for kk = 1:length(A)
                B = TimeStamps==A(kk,1);
                C(B) = 1;    
            end

            for jj = 1:size(WaveForms,2)
            D = WaveForms(:,jj,:);
            E = D(logical(C),:,:);
            F = mean(E(:,:));
            F2(jj,:) = F; 
            end
          
        % From the 4 possible wf selects the one with highest amplitude
        H  = [min(F2,[],2),max(F2,[],2)];
        H2 = H(:,2)-H(:,1);
        H3 = max(H2);
        G  = H2 == H3;
        
        if s == 1
            title(['Waveform - ', neuron(1,1:end-4),' - CH', num2str(G)],'Interpreter','none');
            saveas(gcf, genvarname([neuron(1,1:end-4),...
                 '_',num2str(G), '_wf']), type);
        elseif s ~= 1
        end 

      case 'ISI' % ISI hist -----------------------------------------------
        [ISIh, ISIbins] = ISI_hist(neuron);
                
        if s == 1
            title(['ISI - ', neuron(1,1:end-4)],'Interpreter','none');
            saveas(gcf, genvarname([neuron(1,1:end-4),...
                  '_ISIhist']), type);
        elseif s ~= 1
        end
       
      case 'XCorr' % XCorrelogram -----------------------------------------
        [XCorrVals,XCorrX] = XCorrelogram(neuron, neuron, 'bins', 2);
        
        if s == 1
            title(['XCorr - ', neuron(1,1:end-4)],'Interpreter','none');
            saveas(gcf, genvarname([neuron(1,1:end-4),...
                 '_XCorr']), type);
        elseif s ~= 1
        end 
      
      case 'ACorr' % ACorrelogram -----------------------------------------
        [ACorrVals,ACorrX] = ACorrelogram(neuron);
        
        if s == 1
            title(['ACorr - ', neuron(1,1:end-4)],'Interpreter','none');
            ssaveas(gcf, genvarname([neuron(1,1:end-4),...
                 '_ACorr']), type);
        elseif s ~= 1
        end
  
      otherwise 
        
    end
end

% Save the generated variables  
for ii = 1:size(ToCalculate,2)
    switch ToCalculate{ii}
        case 'wf' % Save waveform calculations ----------------------------
            save(neuron,'wf','-append');
        case 'ISI' % Save ISI hist ----------------------------------------
            save(neuron,'ISIh','ISIbins','-append');
        case 'XCorr' % Save XCorrelogram ----------------------------------
            save(neuron,'XCorrVals','XCorrX','-append');
        case 'ACorr' % Save ACorrelogram ----------------------------------
            save(neuron,'ACorrVals','ACorrX','-append');
        otherwise
    end
end

end