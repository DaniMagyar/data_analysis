function TTLforKS_OEPv063(experiment)
    
info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');
zeroTime = info.Timestamps(1);

states = readNPY([cd '\events\Acquisition_Board-100.Rhythm Data\TTL\states.npy']);
timestamps = double(readNPY([cd '\events\Acquisition_Board-100.Rhythm Data\TTL\timestamps.npy']));

for jj = 1:max(states)
    TTL_ON = find(states==int16(jj));
    TTL_channels{1,jj*2-1} = ['CH' num2str(jj) 'ON_sample_num'];
    TTL_channels{2,jj*2-1} = timestamps(TTL_ON);
    
    TTL_OFF = find(states==-int16(jj));
    TTL_channels{1,jj*2} = ['CH' num2str(jj) 'OFF_sample_num'];
    TTL_channels{2,jj*2} = timestamps(TTL_OFF);
end

switch experiment
    case 'PFC_shock_run_rest'
        shocks = TTL_channels{2,1}-double(zeroTime); % start of each shock, in seconds
    
        for i = 1:10
            shock_motor_rest(:,i) = shocks(i:20:end);
        end
        shock_motor_rest = shock_motor_rest(:);
        shock_motor_rest = sort(shock_motor_rest);
        
        for i = 11:20
            shock_motor_run(:,i-10) = shocks(i:20:end);
        end
        shock_motor_run = shock_motor_run(:);
        shock_motor_run = sort(shock_motor_run);
    
            save('TTLsKS.mat','shocks', 'shock_motor_rest', 'shock_motor_run');

    case 'PFC_shock_run_rest_compare'
        shocks = TTL_channels{2,1}-double(zeroTime); % start of each shock, in seconds
    
        for i = 1:10
            shock_motor_rest(:,i) = shocks(i:30:end);
        end
        shock_motor_rest = shock_motor_rest(:);
        shock_motor_rest = sort(shock_motor_rest);
        
        for i = 11:20
            shock_motor_run_small(:,i-10) = shocks(i:30:end);
        end
        shock_motor_run_small = shock_motor_run_small(:);
        shock_motor_run_small = sort(shock_motor_run_small);

        for i = 21:30
            shock_motor_run_high(:,i-20) = shocks(i:30:end);
        end
        shock_motor_run_high = shock_motor_run_high(:);
        shock_motor_run_high = sort(shock_motor_run_high);

            save('TTLsKS.mat','shocks', 'shock_motor_rest', 'shock_motor_run_small', 'shock_motor_run_high');

%     case 'PFC_BAopto_run_rest'
%         stim = TTL_channels{2,1}-double(zeroTime);
% 
%         for ii = 1:5
%             BA_5Hz_motor_rest_all(:,ii) = stim(ii:20:end);
%         end
%         BA_5Hz_motor_rest_all = BA_5Hz_motor_rest_all(:);
%         BA_5Hz_motor_rest_all = sort(BA_5Hz_motor_rest_all);
% 
%         for jj = 6:10
%             BA_20Hz_motor_rest_all(:,jj-5) = stim(jj:20:end);
%         end
%         BA_20Hz_motor_rest_all = BA_20Hz_motor_rest_all(:);
%         BA_20Hz_motor_rest_all = sort(BA_20Hz_motor_rest_all);
% 
%         for kk = 11:15
%             BA_5Hz_motor_run_all(:,kk-10) = stim(kk:20:end);
%         end
%         BA_5Hz_motor_run_all = BA_5Hz_motor_run_all(:);
%         BA_5Hz_motor_run_all = sort(BA_5Hz_motor_run_all);
% 
%         for ll = 16:20
%             BA_20Hz_motor_run_all(:,ll-15) = stim(ll:20:end);
%         end
%         BA_20Hz_motor_run_all = BA_20Hz_motor_run_all(:);
%         BA_20Hz_motor_run_all = sort(BA_20Hz_motor_run_all);
% 
%         BA_5Hz_motor_rest_first = BA_5Hz_motor_rest_all(1:5:end);
%         BA_20Hz_motor_rest_first = BA_20Hz_motor_rest_all(1:5:end);
%         BA_5Hz_motor_run_first = BA_5Hz_motor_run_all(1:5:end);
%         BA_20Hz_motor_run_first = BA_20Hz_motor_run_all(1:5:end);
% 
%         
%         save('TTLsKS.mat','stim', 'BA_5Hz_motor_rest_all', 'BA_20Hz_motor_rest_all', 'BA_5Hz_motor_run_all', 'BA_20Hz_motor_run_all', ...
%             'BA_5Hz_motor_rest_first', 'BA_20Hz_motor_rest_first', 'BA_5Hz_motor_run_first', 'BA_20Hz_motor_run_first');

    case 'PFC_BAopto_run_rest'
         stim = TTL_channels{2,1}-double(zeroTime);

         rounds = numel(stim)/100; % 1 run+rest cycle contains 2 rounds (10x5+5 rest, then 10x5+5 run)
         rounds_rest = 1:2:rounds;
         rounds_run = 2:2:rounds;

         for ii = rounds_rest
             stim_rest(:,ii) = stim(ii*100-99:ii*100);
         end
         for jj = rounds_run
             stim_run(:,jj) = stim(jj*100-99:jj*100);
         end
         stim_rest = stim_rest(:);
         stim_run = stim_run(:);
         stim_rest(stim_rest==0) = [];
         stim_run(stim_run==0) = [];

         for kk = 1:5
             BA_5Hz_motor_rest_all(:,kk) = stim_rest(kk:10:end);
         end
         BA_5Hz_motor_rest_all = BA_5Hz_motor_rest_all(:);
         BA_5Hz_motor_rest_all = sort(BA_5Hz_motor_rest_all);

         for ll = 6:10
             BA_20Hz_motor_rest_all(:,ll-5) = stim_rest(ll:10:end);
         end
         BA_20Hz_motor_rest_all = BA_20Hz_motor_rest_all(:);
         BA_20Hz_motor_rest_all = sort(BA_20Hz_motor_rest_all);


         for mm = 1:5
             BA_5Hz_motor_run_all(:,mm) = stim_run(mm:10:end);
         end
         BA_5Hz_motor_run_all = BA_5Hz_motor_run_all(:);
         BA_5Hz_motor_run_all = sort(BA_5Hz_motor_run_all);

         for nn = 6:10
             BA_20Hz_motor_run_all(:,nn-5) = stim_run(nn:10:end);
         end
         BA_20Hz_motor_run_all = BA_20Hz_motor_run_all(:);
         BA_20Hz_motor_run_all = sort(BA_20Hz_motor_run_all);

         BA_5Hz_motor_rest_first = BA_5Hz_motor_rest_all(1:5:end);
         BA_20Hz_motor_rest_first = BA_20Hz_motor_rest_all(1:5:end);
         BA_5Hz_motor_run_first = BA_5Hz_motor_run_all(1:5:end);
         BA_20Hz_motor_run_first = BA_20Hz_motor_run_all(1:5:end);

         save('TTLsKS.mat','stim', 'BA_5Hz_motor_rest_all', 'BA_20Hz_motor_rest_all', 'BA_5Hz_motor_run_all', 'BA_20Hz_motor_run_all', ...
             'BA_5Hz_motor_rest_first', 'BA_20Hz_motor_rest_first', 'BA_5Hz_motor_run_first', 'BA_20Hz_motor_run_first');

% plot(BA_5Hz_motor_rest_all, ones(numel(BA_5Hz_motor_rest_all)),'.',BA_5Hz_motor_run_all, ones(numel(BA_5Hz_motor_run_all)),'.', BA_20Hz_motor_rest_all, ones(numel(BA_20Hz_motor_run_all)),'.',BA_20Hz_motor_run_all, ones(numel(BA_5Hz_motor_rest_all)), '.')

end