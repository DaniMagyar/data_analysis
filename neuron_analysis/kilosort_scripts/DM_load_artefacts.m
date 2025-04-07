function [artefacts_sn] = DM_load_artefacts(g)

channel_states = readNPY([cd '\events\' g.OriginalHeader.folder_name '\TTL\states.npy']);
event_sample_numbers = readNPY([cd '\events\' g.OriginalHeader.folder_name '\TTL\sample_numbers.npy']);
data_sample_numbers = readNPY([cd '\continuous\' g.OriginalHeader.folder_name '\sample_numbers.npy']);

for jj = 1:max(channel_states)
    TTL_ON = find(channel_states==int16(jj));
    TTL_channels{1,jj*2-1} = ['CH' num2str(jj) 'ON_sample_num'];
    TTL_channels{2,jj*2-1} = event_sample_numbers(TTL_ON);
    
    TTL_OFF = find(channel_states==-int16(jj));
    TTL_channels{1,jj*2} = ['CH' num2str(jj) 'OFF_sample_num'];
    TTL_channels{2,jj*2} = event_sample_numbers(TTL_OFF);
end

switch g.experiment
    case 'PFC_shock_run_rest'
        artefacts_sn = sort([TTL_channels{2,1}]); % number of samples where stimulation starts
    case 'PFC_shock_run_rest_ONOFF'
        artefacts_sn = sort([TTL_channels{2,1}; TTL_channels{2,2}]); % number of samples where stimulation starts
    case 'PFC_BAopto_run_rest'
         artefacts_sn = sort([TTL_channels{2,7}]);
    case 'PFC_BAopto_laser_remove'
        artefacts_sn = sort([TTL_channels{2,7}; TTL_channels{2,9}]);
    case 'PFC_BAopto_laser_ON_OFF_remove'
        artefacts_sn = sort([TTL_channels{2,7}; TTL_channels{2,8}; TTL_channels{2,9}; TTL_channels{2,10};]);
    otherwise 
        error('unknown experiment')
end
artefacts_sn = (artefacts_sn - data_sample_numbers(1)) +1; % normalize sample number to 'incoming_data'; start from 1 instead of 0

if g.firstSec > 0
    artefacts_sn = artefacts_sn(artefacts_sn>g.firstSec*30000);
end
if g.lastSec > 0
    artefacts_sn = artefacts_sn(artefacts_sn<g.lastSec*30000);
end
if g.firstSec > 0
    artefacts_sn = artefacts_sn - firstEL; % normalize pulses_sn to fristEL
end