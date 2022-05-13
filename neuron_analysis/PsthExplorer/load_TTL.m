function [ttl] = load_TTL(Stim)
load([cd '\TTLsKS.mat']); % #ok<LOAD>
    switch Stim
        case 'BA_500'
            ttl = BA_500;
        case 'BA_500_5Hz'
            ttl = BA_500_5Hz;
        case 'BA_50'
            ttl = BA_50;
        case 'TTL_500'
            ttl = TTL;
        case 'TTL_50'
            ttl = TTL(1:10:end);
        case 'BA_25'
            ttl = BA_25;
        case 'BA_25_5Hz'
            ttl = BA_25_5Hz;
    %            ttl = BA_25_5Hz(11:25);
        case 'BA_25_10Hz'
            ttl = BA_25_10Hz;
        case 'TO_25'
            ttl = TO_25;
        case 'TO_25_5Hz'
            ttl = TO_25_5Hz;
    %            ttl = TO_25_5Hz(1:10);
        case 'TO_25_10Hz'
            ttl = TO_25_10Hz;
        case 'BA_250'
            ttl = BA_250;
        case 'BA_250_5Hz'
            ttl = BA_250_5Hz;
    %             ttl = BA_250_5Hz(1:125);
    %             ttl = [ttl; BA_250_5Hz(1:10:end)];
    %             ttl = [ttl; BA_250_5Hz(2:10:end)];
    %             ttl = [ttl; BA_250_5Hz(3:10:end)];
    %             ttl = [ttl; BA_250_5Hz(4:10:end)];
    %             ttl = sort(ttl);
        case 'BA_250_10Hz'
            ttl = BA_250_10Hz;
        case 'TO_250'
            ttl = TO_250;
        case 'TO_250_5Hz'
            ttl = TO_250_5Hz;
    %            ttl = TO_250_5Hz(1:125);
        case 'TO_250_10Hz'
            ttl = TO_250_10Hz;
        case 'BA_500_5Hz_10Hz'
            ttl = BA_250_5Hz;
            ttl = [ttl; BA_250_10Hz];
            ttl = sort(ttl);
        case 'SK'
            ttl = TTL; % The SK TTL is already called TTL.
        case 'shock_only'
            ttl = shock_only;
        case 'shock_inh'
            ttl = shock_inh;
    end