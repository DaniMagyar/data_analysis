rawDiff = diff(TAG_all);
roundDiff = round(rawDiff,3);
[GC,GR] = groupcounts(roundDiff);
pulseIndex = NaN(max(GC), numel(GR(GR<1)));
pulseMatrix = NaN(max(GC), numel(GR(GR<1)));

freqs = GR(GR<1); % assuming optotag when pulsediff <1
freqStr=string(freqs);
freqCell = cellstr(freqStr);
freqNamesForSave = matlab.lang.makeValidName(freqCell);

for ii = 1:numel(freqs) 
    pulseNum = numel(roundDiff(roundDiff==GR(ii)));
    pulseIndex(1:pulseNum,ii)= find(roundDiff==GR(ii));
    pulseMatrix(1:pulseNum,ii) = TAG_all(pulseIndex(1:pulseNum,ii));
    saveStruct.(freqNamesForSave{ii}) = TAG_all(pulseIndex(1:pulseNum,ii));
end

save('optoTTL.mat','-struct','saveStruct')

% The last pulses of each train are missing here. Can do it later
