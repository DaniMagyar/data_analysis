function load_digital_ttl(experiment)
    
% Input experiment:
% 'FFI'
% 'optotagging'


    channel_states = readNPY('channel_states.npy');
    channels = readNPY('channels.npy');
    timestamps = double(readNPY('timestamps.npy'));
    for ii = 1:max(channels)
        TTL_ON = find(channel_states==ii);
        TTL_channels{1,ii} = ['CH' num2str(ii)];
        TTL_channels{2,ii} = timestamps(TTL_ON);
    end
    switch experiment
        case 'FFI'
            pulses = (TTL_channels{2,3})/30000; % start of each pulse, in seconds
            for i = 1:10
                BA_250(:,i) = pulses(i:20:end);
            end
            BA_250 = BA_250(:);
            BA_250 = sort(BA_250);
            BA_25 = BA_250(1:10:end);
    
            for i = 11:20
                TO_250(:,i) = pulses(i:20:end);
            end
            TO_250 = TO_250(:);
            TO_250 = sort(TO_250);
            TO_250(1:250) = [];
            TO_25 = TO_250(1:10:end);
    
            save('TTLs.mat','BA_250','BA_25','TO_250','TO_25');

        case 'FFI_5Hz_10Hz'
            % 5Hz pulse trains
            pulses = (TTL_channels{2,3})/30000; % start of each pulse, in seconds
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

            save('TTLs.mat','BA_250_5Hz','BA_25_5Hz','TO_250_5Hz','TO_25_5Hz',...
                'BA_250_10Hz','BA_25_10Hz','TO_250_10Hz','TO_25_10Hz');

            
        case 'optotagging'
            BA_500 = (TTL_channels{2,3})/30000; % start of each pulse, in seconds
            BA_50 = BA_500(1:10:end);
            TAG_all = (TTL_channels{2,2})/30000;
            save('TTLs.mat','BA_500','BA_50','TAG_all');
    end