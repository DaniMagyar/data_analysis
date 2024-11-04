function TTLforKS_OEPv063(experiment)

firstSec = 0; % if Preprocessing had a "firstSec", then use that value to synchronize TTLs. 
              % Remove negative times manually!!

    
info = load_open_ephys_binary([cd '\structure.oebin'],'continuous',1, 'mmap');
zeroTime = info.Timestamps(1) + firstSec;

states = readNPY([cd '\events\' info.Header.folder_name '\TTL\states.npy']);
timestamps = double(readNPY([cd '\events\' info.Header.folder_name '\TTL\timestamps.npy']));

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

    case 'PFC_BAopto_shock'
        shocks = TTL_channels{2,9}-double(zeroTime); % start of each shock, in seconds
    
        for i = 1:5
            BA_20Hz_only(:,i) = shocks(i:10:end);
        end
        BA_20Hz_only = BA_20Hz_only(:);
        BA_20Hz_only = sort(BA_20Hz_only);
        
        for i = 6:10
            BA_20Hz_shock(:,i-5) = shocks(i:10:end);
        end
        BA_20Hz_shock = BA_20Hz_shock(:);
        BA_20Hz_shock = sort(BA_20Hz_shock);

        BA_20Hz_only_first = BA_20Hz_only(1:5:end);
        BA_20Hz_shock_first = BA_20Hz_shock(1:5:end);
    
            save('TTLsKS.mat','shocks', 'BA_20Hz_only', 'BA_20Hz_shock', 'BA_20Hz_only_first', 'BA_20Hz_shock_first');

    case 'M2_footshock_only'
         shocks = TTL_channels{2,1}-double(zeroTime);
         save('TTLsKS.mat','shocks');

    case 'motor'
        motorON = TTL_channels{2,11}-double(zeroTime);
        motorDiff = diff(motorON);
        startIdx = [1; (find(motorDiff > std(motorDiff)*3)+1)];
        motorStart = motorON(startIdx);

        
        save('TTLsKS.mat', 'motorStart',  "-append")
        



%MD199 custom: 
% RED_ON = TTL_channels{2,7}-double(zeroTime);
% RED_ON_singles = RED_ON(1:400);
% RED_ON_singles(401:600) = RED_ON(501:700);
% plot(diff(RED_ON_singles))
% BA_singles_all = RED_ON_singles(1:2:end);
% TO_singles_all = RED_ON_singles(2:2:end);
% save('TTLsKS.mat','BA_singles_all','TO_singles_all')
% RED_ON_trains = RED_ON(401:500);
% RED_ON_trains(101:800) = RED_ON(701:1400);
% BA_trains_all = vertcat(RED_ON_trains(1:10:end),RED_ON_trains(2:10:end),RED_ON_trains(3:10:end),RED_ON_trains(4:10:end),RED_ON_trains(5:10:end));
% TO_trains_all = vertcat(RED_ON_trains(6:10:end),RED_ON_trains(7:10:end),RED_ON_trains(8:10:end),RED_ON_trains(9:10:end),RED_ON_trains(10:10:end));
% BA_trains_all = sort(BA_trains_all);
% TO_trains_all = sort(TO_trains_all);
% BA_trains_first = BA_trains_all(1:5:end);
% TO_trains_first = TO_trains_all(1:5:end);

% VALIDATE: 
% plot(BA_singles_all,ones(300), 'b.')
% hold on
% plot(TO_singles_all,ones(300),'r.')
% plot(BA_trains_all,ones(400), 'g.')
% plot(TO_trains_all,ones(400), 'm.')

% MD200 custom: 
% RED_ON = TTL_channels{2,7}-double(zeroTime);
% RED_ON_singles = RED_ON(1:400);
% plot(diff(RED_ON_singles))
% BA_singles_all = RED_ON_singles(1:2:end);
% TO_singles_all = RED_ON_singles(2:2:end);
% RED_ON_trains = RED_ON(401:end);
% BA_only_trains_all = vertcat(RED_ON_trains(1:15:end),RED_ON_trains(2:15:end),RED_ON_trains(3:15:end),RED_ON_trains(4:15:end),RED_ON_trains(5:15:end));
% BA_shock_trains_all = vertcat(RED_ON_trains(6:15:end),RED_ON_trains(7:15:end),RED_ON_trains(8:15:end),RED_ON_trains(9:15:end),RED_ON_trains(10:15:end));
% TO_shock_trains_all = vertcat(RED_ON_trains(11:15:end),RED_ON_trains(12:15:end),RED_ON_trains(13:15:end),RED_ON_trains(14:15:end),RED_ON_trains(15:15:end));
% BA_only_trains_all = sort(BA_only_trains_all);
% BA_shock_trains_all = sort(BA_shock_trains_all);
% TO_shock_trains_all = sort(TO_shock_trains_all);
% BA_only_trains_first = BA_only_trains_all(1:5:end);
% BA_only_trains_first = BA_only_trains_all(1:5:end);
% BA_shock_trains_first = BA_shock_trains_all(1:5:end);
% TO_shock_trains_first = TO_shock_trains_all(1:5:end);
% shocks_ON = TTL_channels{2,1}-double(zeroTime);
% shock_only = shocks_ON(1:3:end);
% save('TTLsKS.mat','BA_singles_all','TO_singles_all', 'BA_only_trains_all', 'BA_only_trains_first', 'BA_shock_trains_all', 'BA_shock_trains_first', 'TO_shock_trains_all', 'TO_shock_trains_first', 'RED_ON', 'RED_ON_singles', 'RED_ON_trains', 'shock_only')


% MD202 
% RED_ON = TTL_channels{2,7}-double(zeroTime);
% RED_ON_trains = RED_ON(1:2000);
% RED_ON_trains(2001:2500) = RED_ON(2451:2950);
% BA_trains_all = vertcat(RED_ON_trains(1:10:end),RED_ON_trains(2:10:end),RED_ON_trains(3:10:end),RED_ON_trains(4:10:end),RED_ON_trains(5:10:end));
% TO_trains_all = vertcat(RED_ON_trains(6:10:end),RED_ON_trains(7:10:end),RED_ON_trains(8:10:end),RED_ON_trains(9:10:end),RED_ON_trains(10:10:end));
% BA_trains_all = sort(BA_trains_all);
% TO_trains_all = sort(TO_trains_all);
% BA_trains_first = BA_trains_all(1:5:end);
% TO_trains_first = TO_trains_all(1:5:end);

% MD204
% RED_ON = TTL_channels{2,7}-double(zeroTime);
% RED_ON_trains = RED_ON(1:2500);
% BA_trains_all = vertcat(RED_ON_trains(1:20:end),RED_ON_trains(2:20:end),RED_ON_trains(3:20:end),RED_ON_trains(4:20:end),RED_ON_trains(5:20:end),RED_ON_trains(6:20:end),RED_ON_trains(7:20:end),RED_ON_trains(8:20:end),RED_ON_trains(9:20:end),RED_ON_trains(10:20:end));
% TO_trains_all = vertcat(RED_ON_trains(11:20:end),RED_ON_trains(12:20:end),RED_ON_trains(13:20:end),RED_ON_trains(14:20:end),RED_ON_trains(15:20:end),RED_ON_trains(16:20:end),RED_ON_trains(17:20:end),RED_ON_trains(18:20:end),RED_ON_trains(19:20:end),RED_ON_trains(20:20:end));
% BA_trains_all = sort(BA_trains_all);
% TO_trains_all = sort(TO_trains_all);
% BA_trains_first = BA_trains_all(1:10:end);
% TO_trains_first = TO_trains_all(1:10:end);
% save('TTLsKS.mat', 'BA_trains_all', 'TO_trains_all', 'BA_trains_first', 'TO_trains_first')

% MD205
% RED_ON = TTL_channels{2,7}-double(zeroTime);
% RED_ON_trains = RED_ON(1:2000);
% BA_trains_all = vertcat(RED_ON_trains(1:20:end),RED_ON_trains(2:20:end),RED_ON_trains(3:20:end),RED_ON_trains(4:20:end),RED_ON_trains(5:20:end),RED_ON_trains(6:20:end),RED_ON_trains(7:20:end),RED_ON_trains(8:20:end),RED_ON_trains(9:20:end),RED_ON_trains(10:20:end));
% TO_trains_all = vertcat(RED_ON_trains(11:20:end),RED_ON_trains(12:20:end),RED_ON_trains(13:20:end),RED_ON_trains(14:20:end),RED_ON_trains(15:20:end),RED_ON_trains(16:20:end),RED_ON_trains(17:20:end),RED_ON_trains(18:20:end),RED_ON_trains(19:20:end),RED_ON_trains(20:20:end));
% BA_trains_all = sort(BA_trains_all);
% TO_trains_all = sort(TO_trains_all);
% BA_trains_first = BA_trains_all(1:10:end);
% TO_trains_first = TO_trains_all(1:10:end);
% save('TTLsKS.mat', 'BA_trains_all', 'TO_trains_all', 'BA_trains_first', 'TO_trains_first')

% MD233
% shocks = TTL_channels{2,1}-double(zeroTime);
% shockdiff = diff(shocks);
% shockdiff1 = [1; shockdiff];
% shocks1 = shocks(shockdiff1>=1);
% shocks2 = reshape(shocks1, [10,6])
% shock_motor_rest = [shocks2(:,2);shocks2(:,3);shocks2(:,5)]  % olvasd el a jegyzokonyvet, fel volt cserelve az elejen a run-rest
% shock_motor_run = [shocks2(:,1);shocks2(:,4);shocks2(:,6)]
% save('TTLsKS.mat','shocks', 'shock_motor_rest', 'shock_motor_run');'

% MD236_rec001
% RED_ON = TTL_channels{2,7}-double(zeroTime);
% BA_trains_all = RED_ON(1:200);
% BA_trains_first = RED_ON(1:5:200);
% RED_ON2 = RED_ON(201:600);
% RED_ON3 = reshape(RED_ON2, [5,80]);
% BAstimON_BeforeShock = reshape(RED_ON3(:, 1:2:end), [1,200])';
% BAstimON_AfterShock = reshape(RED_ON3(:, 2:2:end), [1,200])';
% BAstimON_BeforeShock_first = BAstimON_BeforeShock(1:5:end);
% BAstimON_AfterShock_first = BAstimON_AfterShock(1:5:end);
% save('TTLsKS.mat', 'BA_trains_all', 'BA_trains_first', 'BAstimON_BeforeShock', 'BAstimON_BeforeShock_first', 'BAstimON_AfterShock', 'BAstimON_AfterShock_first')
% 


% MD236_rec002 
% RED_ON = TTL_channels{2,7}-double(zeroTime);
% RED_ON2 = reshape(RED_ON, [5,160]);
% BA_trains_first = RED_ON2(1,:)';
% BA_trains_all = RED_ON;
% save('TTLsKS.mat', 'BA_trains_all', 'TO_trains_all', 'BA_trains_first', 'TO_trains_first')

% MD237
% RED_ON = TTL_channels{2,7}-double(zeroTime);
% BA_trains_all = RED_ON(1:400);
% BA_trains_first = RED_ON(1:5:400);
% shocks = TTL_channels{2,1}-double(zeroTime);
% shocks2 = reshape(shocks, [10,8]);
% shock_BAstim = vertcat(shocks2(:, 1:2:end));
% shock_only = reshape(vertcat(shocks2(:, 2:2:end)), [numel(shocks2(:, 2:2:end)),1]);
% shock_BAstim = reshape(vertcat(shocks2(:, 1:2:end)), [numel(shocks2(:, 1:2:end)),1]);
% RED_ON2 = RED_ON(401:end);
% RED_ON3 = reshape(RED_ON2, [100,8]);
% BAstimON_BeforeAfterShock = reshape(RED_ON3(:,1:2:end), [numel(RED_ON3(:,1:2:end)),1]);
% BAstimOFF_BeforeAfterShock = reshape(RED_ON3(:,2:2:end), [numel(RED_ON3(:,2:2:end)),1]);
% BAstimON_BeforeAfterShock2 = reshape(BAstimON_BeforeAfterShock, [5,80]);
% BAstimOFF_BeforeAfterShock2 = reshape(BAstimOFF_BeforeAfterShock, [5,80]);
% BAstimON_BeforeShock = reshape(BAstimON_BeforeAfterShock2(:, 1:2:end), [1,200])';
% BAstimON_AfterShock = reshape(BAstimON_BeforeAfterShock2(:, 2:2:end), [1,200])';
% BAstimOFF_BeforeShock = reshape(BAstimOFF_BeforeAfterShock2(:, 1:2:end), [1,200])';
% BAstimOFF_AfterShock = reshape(BAstimOFF_BeforeAfterShock2(:, 2:2:end), [1,200])';
% save('TTLsKS.mat', 'BA_trains_all', 'BA_trains_first', 'shock_only', 'shock_BAstim', 'BAstimON_BeforeShock','BAstimON_AfterShock', 'BAstimOFF_BeforeShock', 'BAstimOFF_AfterShock')
% csvwrite('events.csv',BA_trains_first)


% MD238
% shocks = TTL_channels{2,1}-double(zeroTime);
% sound_all = TTL_channels{2,5}-double(zeroTime);
% sound_all_diff = diff(sound_all);
% sound_all_diff1 = [0.5; sound_all_diff];
% sound_all1 = sound_all(sound_all_diff1>=0.5);
% sound_habit_all = sound_all1(1:240);
% sound_habit_resh = reshape(sound_habit_all, [30,8]);
% tone_habit_all = reshape(sound_habit_resh(:,1:2:end), [120,1]);
% noise_habit_all = reshape(sound_habit_resh(:,2:2:end), [120,1]);
% sound_cond_all = sound_all1(241:600);
% sound_cond_resh = reshape(sound_cond_all, [30,12]);
% tone_cond_all = reshape(sound_cond_resh(:,1:2:end), [180,1]);
% noise_cond_all = reshape(sound_cond_resh(:,2:2:end), [180,1]);
% sound_recall_all = sound_all1(601:840);
% sound_recall_resh = reshape(sound_recall_all, [30,8]);
% tone_recall_all = reshape(sound_recall_resh(:,1:2:end), [120,1]);
% noise_recall_all = reshape(sound_recall_resh(:,2:2:end), [120,1]);
% tone_habit_first = tone_habit_all(1:30:end);
% noise_habit_first = noise_habit_all(1:30:end);
% tone_cond_first = tone_cond_all(1:30:end);
% noise_cond_first = noise_cond_all(1:30:end);
% tone_recall_first = tone_recall_all(1:30:end);
% noise_recall_first = noise_recall_all(1:30:end);
% save('TTLsKS.mat', 'shocks', 'sound_all1', 'tone_habit_all', 'tone_cond_all', 'tone_recall_all', 'noise_habit_all', 'noise_cond_all', 'noise_recall_all', 'tone_habit_first', 'tone_cond_first', 'tone_recall_first', 'noise_habit_first', 'noise_cond_first', 'noise_recall_first')

% MD243, MD250, MD251, MD254
LED_synch = TTL_channels{2,3}-double(zeroTime);
shocks = TTL_channels{2,1}-double(zeroTime);
sound_all = TTL_channels{2,5}-double(zeroTime);
sound_all_diff = diff(sound_all);
sound_all_diff1 = [0.5; sound_all_diff];
sound_all1 = sound_all(sound_all_diff1>=0.5);
sound_habit_all = sound_all1(1:100);
sound_habit_resh = reshape(sound_habit_all, [5,20]);
tone_habit_all = reshape(sound_habit_resh(:,1:2:end), [50,1]);
noise_habit_all = reshape(sound_habit_resh(:,2:2:end), [50,1]);
sound_cond_all = sound_all1(101:300);
sound_cond_resh = reshape(sound_cond_all, [5,40]);
tone_cond_all = reshape(sound_cond_resh(:,1:2:end), [100,1]);
noise_cond_all = reshape(sound_cond_resh(:,2:2:end), [100,1]);
sound_recall_all = sound_all1(301:400);
sound_recall_resh = reshape(sound_recall_all, [5,20]);
tone_recall_all = reshape(sound_recall_resh(:,1:2:end), [50,1]);
noise_recall_all = reshape(sound_recall_resh(:,2:2:end), [50,1]);
tone_habit_first = tone_habit_all(1:5:end);
noise_habit_first = noise_habit_all(1:5:end);
tone_cond_first = tone_cond_all(1:5:end);
noise_cond_first = noise_cond_all(1:5:end);
tone_recall_first = tone_recall_all(1:5:end);
noise_recall_first = noise_recall_all(1:5:end);
save('TTLsKS.mat', 'shocks', 'sound_all1', 'tone_habit_all', 'tone_cond_all', 'tone_recall_all', 'noise_habit_all', 'noise_cond_all', 'noise_recall_all', 'tone_habit_first', 'tone_cond_first', 'tone_recall_first', 'noise_habit_first', 'noise_cond_first', 'noise_recall_first', 'LED_synch')
plot(tone_habit_all, ones(numel(tone_habit_all)),'.',noise_habit_all, ones(numel(noise_habit_all)),'.',tone_cond_all, ones(numel(tone_cond_all)),'.',noise_cond_all, ones(numel(noise_cond_all)),'.',tone_recall_all, ones(numel(tone_recall_all)),'.',noise_recall_all, ones(numel(noise_recall_all)),'.')

% interchanged CS, MD252, MD253, MD268, MD277, MD278
LED_synch = TTL_channels{2,3}-double(zeroTime);
shocks = TTL_channels{2,1}-double(zeroTime);
sound_all = TTL_channels{2,5}-double(zeroTime);
sound_all_diff = diff(sound_all);
sound_all_diff1 = [0.5; sound_all_diff];
sound_all1 = sound_all(sound_all_diff1>=0.5);
sound_habit_all = sound_all1(1:100);
sound_habit_resh = reshape(sound_habit_all, [5,20]);
tone_habit_all = reshape(sound_habit_resh(:,2:2:end), [50,1]);
noise_habit_all = reshape(sound_habit_resh(:,1:2:end), [50,1]);
sound_cond_all = sound_all1(101:300);
sound_cond_resh = reshape(sound_cond_all, [5,40]);
tone_cond_all = reshape(sound_cond_resh(:,2:2:end), [100,1]);
noise_cond_all = reshape(sound_cond_resh(:,1:2:end), [100,1]);
sound_recall_all = sound_all1(301:400);
sound_recall_resh = reshape(sound_recall_all, [5,20]);
tone_recall_all = reshape(sound_recall_resh(:,2:2:end), [50,1]);
noise_recall_all = reshape(sound_recall_resh(:,1:2:end), [50,1]);
tone_habit_first = tone_habit_all(1:5:end);
noise_habit_first = noise_habit_all(1:5:end);
tone_cond_first = tone_cond_all(1:5:end);
noise_cond_first = noise_cond_all(1:5:end);
tone_recall_first = tone_recall_all(1:5:end);
noise_recall_first = noise_recall_all(1:5:end);
save('TTLsKS.mat', 'shocks', 'sound_all1', 'tone_habit_all', 'tone_cond_all', 'tone_recall_all', 'noise_habit_all', 'noise_cond_all', 'noise_recall_all', 'tone_habit_first', 'tone_cond_first', 'tone_recall_first', 'noise_habit_first', 'noise_cond_first', 'noise_recall_first', 'LED_synch')
plot(tone_habit_all, ones(numel(tone_habit_all)),'.',noise_habit_all, ones(numel(noise_habit_all)),'.',tone_cond_all, ones(numel(tone_cond_all)),'.',noise_cond_all, ones(numel(noise_cond_all)),'.',tone_recall_all, ones(numel(tone_recall_all)),'.',noise_recall_all, ones(numel(noise_recall_all)),'.')






% plot(BA_5Hz_motor_rest_all, ones(numel(BA_5Hz_motor_rest_all)),'.',BA_5Hz_motor_run_all, ones(numel(BA_5Hz_motor_run_all)),'.', BA_20Hz_motor_rest_all, ones(numel(BA_20Hz_motor_run_all)),'.',BA_20Hz_motor_run_all, ones(numel(BA_5Hz_motor_rest_all)), '.')



% MD246
pulses = TTL_channels{2,7}-double(zeroTime);
pulses_firsts = pulses(1:5:end);
BA_trains_first = pulses_firsts(1:2:end);
TO_trains_first = pulses_firsts(2:2:end);

pulses_resh = reshape(pulses, [5,200]);
BA_trains_all = reshape(pulses_resh(:,1:2:end), [500,1]);
TO_trains_all = reshape(pulses_resh(:,2:2:end), [500,1]);
save('TTLsKS.mat', 'BA_trains_all', 'TO_trains_all', 'BA_trains_first', 'TO_trains_first')

pulses = TTL_channels{2,9}-double(zeroTime);
BLUE_only = pulses(1:2:end);
save('TTLsKS.mat', 'BLUE_only', '-append')

% MD259
shocks = TTL_channels{2,1}-double(zeroTime);
pulses = TTL_channels{2,7}-double(zeroTime);
save('TTLsKS.mat', 'shocks', 'pulses')
csvwrite('events.csv',pulses)

% MD271
shocks = TTL_channels{2,1}-double(zeroTime);
shocks_only  = shocks(1:2:end);
shocks_light = shocks(2:2:end);
pulses = TTL_channels{2,7}-double(zeroTime);
save('TTLsKS.mat', 'shocks', 'shocks_only', 'shocks_light', 'pulses')
csvwrite('events.csv',shocks_light)


%MD243 SYNCH LED TS
TTL_synch = TTL_channels{2,3}-double(zeroTime);
save('MD243_TTL_synch.mat', 'TTL_synch')

end