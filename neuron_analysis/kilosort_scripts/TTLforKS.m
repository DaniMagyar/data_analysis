function TTLforKS(experiment)

% Navigate to the folder containing structure.oebin

info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');
zeroTime = info.Timestamps(1);

channel_states = readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\channel_states.npy']);
channels = readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\channels.npy']);
timestamps = double(readNPY([cd '\events\Rhythm_FPGA-100.0\TTL_1\timestamps.npy']));

for ii = 1:max(channels)
    TTL_ON = find(channel_states==ii);
    TTL_channels{1,ii} = ['CH' num2str(ii)];
    TTL_channels{2,ii} = timestamps(TTL_ON);
end
switch experiment
    case 'BA_50_5Hz'
        pulses = (TTL_channels{2,3}-double(zeroTime))/30000; % start of each pulse, in seconds
        pulses_5Hz = pulses(1:500);
        BA_50_5Hz = pulses_5Hz(1:10:end);
        save('TTLsKS.mat','BA_50_5Hz');
    case 'PFC_ChETA'
        pulses = (TTL_channels{2,3}-double(zeroTime))/30000; % start of each pulse, in seconds
        ChETA_50_20Hz = pulses(1:40:end);
        save('TTLsKS.mat','ChETA_50_20Hz');
    case 'PFC_shock_laser'
        shocks = (TTL_channels{2,2}-double(zeroTime))/30000;
        shock_only = shocks(1:2:end);
        shock_inh = shocks(2:2:end);
        save('TTLsKS.mat','shock_only','shock_inh');
    case 'FFI_5Hz_10Hz'
        pulses = (TTL_channels{2,3}-double(zeroTime))/30000; % start of each pulse, in seconds
        pulses_5Hz = pulses(1:500);
        for i = 1:10
            BA_250_5Hz(:,i) = pulses_5Hz(i:20:end);
        end
        BA_250_5Hz = BA_250_5Hz(:);
        BA_250_5Hz = sort(BA_250_5Hz);
        BA_25_5Hz = BA_250_5Hz(1:10:end);
        
        for i = 11:20
            TO_250_5Hz(:,i-10) = pulses_5Hz(i:20:end);
        end
        TO_250_5Hz = TO_250_5Hz(:);
        TO_250_5Hz = sort(TO_250_5Hz);
        TO_25_5Hz = TO_250_5Hz(1:10:end);
        
        %10Hz pulse trains
        pulses_10Hz = pulses(501:end);
        for i = 1:10
            BA_250_10Hz(:,i) = pulses_10Hz(i:20:end);
        end
        BA_250_10Hz = BA_250_10Hz(:);
        BA_250_10Hz = sort(BA_250_10Hz);
        BA_25_10Hz = BA_250_10Hz(1:10:end);
        
        for i = 11:20
            TO_250_10Hz(:,i-10) = pulses_10Hz(i:20:end);
        end
        TO_250_10Hz = TO_250_10Hz(:);
        TO_250_10Hz = sort(TO_250_10Hz);
        TO_25_10Hz = TO_250_10Hz(1:10:end);
        
        save('TTLsKS.mat','BA_250_5Hz','BA_25_5Hz','TO_250_5Hz','TO_25_5Hz',...
            'BA_250_10Hz','BA_25_10Hz','TO_250_10Hz','TO_25_10Hz');

    case 'footshock'
        shocks = (TTL_channels{2,2}-double(zeroTime))/30000; % start of each shock, in seconds
        shock_only = shocks(1:2:end);
        shock_inh = shocks(2:2:end);
        save('TTLsKS.mat','shock_only','shock_inh','inhibitions');

    case 'M2_footshock_only'
        shocks = (TTL_channels{2,1}-double(zeroTime))/30000; % start of each shock, in seconds
        save('TTLsKS.mat','shocks');
end