function [incoming_data] = DM_detrend_data(varargin)
% detrending data
prs = inputParser;
addRequired(prs,'g',@isstruct)
addRequired(prs,'incoming_data',@isnumeric)
addParameter(prs,'extracut', 10, @isnumeric) % samplepoints before stim onset. (10 ~ 0.3ms) for appropriate linspace. (NOT IN MS!!)
addParameter(prs,'breakpoints', [4 10 50], @isnumeric) % in ms, default 4 10 50 ms. 
parse(prs,varargin{:})

g = prs.Results.g;
incoming_data = prs.Results.incoming_data;
artefacts_sn = double(DM_load_artefacts(g));

% for ch = 1:g.numChannels
%     for ii = 1:size(artefacts_sn,1)
%         idx1 = artefacts_sn(ii);
%         idx2 = round(artefacts_sn(ii)+(g.detrend(1)*g.fs));
%         idx3 = round(artefacts_sn(ii)+(g.detrend(2)*g.fs));
%         incoming_data(ch,idx2:idx3) = int16(detrend(double(incoming_data(ch,idx2:idx3)),3)); % sample number at 30kHz
%         temp2 = linspace(double(incoming_data(ch,idx1)), double(incoming_data(ch, idx2+1)),idx2-idx1+1); % linspacing between artefact start and detrend begining
%         incoming_data(ch,idx1:idx2) = int16(temp2);
%     end    
%     disp(['Channel ' num2str(ch) ' detrended'])
% end  

loadBar = waitbar(0,'Detrending data...');
for ch = 1:g.numChannels
    for ii = 1:size(artefacts_sn,1)
        idx1 = artefacts_sn(ii);
        idx2 = round(artefacts_sn(ii)+(g.detrend(1)*g.fs));
        idx3 = round(artefacts_sn(ii)+(g.detrend(2)*g.fs));
        detrended_data = detrend(double(incoming_data(ch,idx2:idx3)),3, prs.Results.breakpoints*(g.fs/1000)); % sample number at 30kHz
        detrended_data = detrended_data + (double(incoming_data(ch,idx3+1)) - detrended_data(end));
        incoming_data(ch,idx2:idx3) = int16(detrended_data);
        temp2 = linspace(double(incoming_data(ch,idx1-prs.Results.extracut)), detrended_data(1),idx2-idx1+1); % linspacing between artefact start and detrend begining
        incoming_data(ch,idx1:idx2) = int16(temp2);
    end    
    waitbar(ch/g.numChannels,loadBar);
end
close(loadBar)