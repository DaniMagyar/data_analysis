function [incoming_data] = DM_detrend_data(varargin)
% detrending data
prs = inputParser;
addRequired(prs,'g',@isstruct)
addRequired(prs,'incoming_data',@isnumeric)
parse(prs,varargin{:})

g = prs.Results.g;
incoming_data = prs.Results.incoming_data;
artefacts_sn = double(DM_load_artefacts(g));

for ch = 1:g.numChannels
    for ii = 1:size(artefacts_sn,1)
        idx1 = artefacts_sn(ii);
        idx2 = round(artefacts_sn(ii)+(g.detrend(1)*g.fs));
        idx3 = round(artefacts_sn(ii)+(g.detrend(2)*g.fs));
        incoming_data(ch,idx2:idx3) = int16(detrend(double(incoming_data(ch,idx2:idx3)),3)); % sample number at 30kHz
        temp2 = linspace(double(incoming_data(ch,idx1)), double(incoming_data(ch, idx2+1)),idx2-idx1+1);
        incoming_data(ch,idx1:idx2) = int16(temp2);
    end    
end  