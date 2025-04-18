function [channel_map] = DM_load_channel_map(channel_map)

switch channel_map
    case 'none'
        channel_map = 1:64;
        disp('No channel map selected')
    case 'A1x32'
        channel_map = [1 17 16 32 3 19 14 30 9 25 10 20 8 24 2 29 7 26 15 21 11 23 12 28 6 18 13 22 5 27 4 31];
        disp('Selected channel map: A1x32')
    case 'Buzs32'
        channel_map = [31 26 27 21 22 23 18 28 29 17 24 32 20 19 25 30 1 2 16 8 3 10 14 9 7 4 15 5 11 13 12 6];
        disp('Selected channel map: Buzs32')
    case 'Cambridge64_H2_Adpt_NNx'
        channel_map = [3 1 31 28 30 23 7 5 25 6 8 29 27 32 2 4 9 11 14 16 18 20 22 24 26 21 19 17 15 13 12 10 56 54 51 49 47 45 43 41 39 44 46 48 50 52 53 55 62 64 34 37 35 42 58 60 40 59 57 36 38 33 63 61];
        disp('Selected channel map: Cambridge64_H2_Adpt_NNx')
    case 'Cambridge64_H2_Adpt_Cambridge'
        channel_map = [14 13 12 7 6 57 1 2 9 15 16 11 10 5 4 3 49 50 51 52 53 54 55 56 8 58 59 60 61 62 63 64 48 47 46 45 44 43 42 41 25 39 38 37 36 35 34 33 19 20 21 26 27 40 32 31 24 18 17 22 23 28 29 30];
        disp('Selected channel map: Cambridge64_H2_Adpt_Cambridge')
    case 'Cambridge64_H7_Adpt_NNx'
        channel_map = [3 10 1 12 31 13 28 15 30 17 23 19 7 21 5 26 25 24 6 22 8 20 29 18 27 16 32 14 2 11 4 9 56 61 54 63 51 33 49 38 47 36 45 57 43 59 41 40 39 60 44 58 46 42 48 35 50 37 52 34 53 64 55 62];
        disp('Selected channel map: Cambridge64_H7_Adpt_NNx')
    case 'Cambridge64_H7_Adpt_Cambridge'
        channel_map = [14 64 13 63 12 62 7 61 6 60 57 59 1 58 2 8 9 56 15 55 16 54 11 53 10 52 5 51 4 50 3 49 48 30 47 29 46 28 45 23 44 22 43 17 42 18 41 24 25 31 39 32 38 40 37 27 36 26 35 21 34 20 33 19];
        disp('Selected channel map: Cambridge64_H7_Adpt_Cambridge')
    case 'Cambridge64_P1'
        channel_map = [23 3 7 1 5 31 25 28 6 30 8 29 27 32 2 4 20 9 22 11 24 14 26 16 21 18 19 17 15 13 12 10 45 56 43 54 41 51 39 49 44 47 46 48 50 52 53 55 42 62 58 64 60 34 40 37 59 35 57 36 38 33 63 61];
        disp('Selected channel map: Cambridge64_P1')
    case 'Cambridge64_P1_Adpt_Cambridge'
        channel_map = [57 14 1 13 2 12 9 7 15 6 16 11 10 5 4 3 54 49 55 50 56 51 8 52 58 53 59 60 61 62 63 64 43 48 42 47 41 46 25 45 39 44 38 37 36 35 34 33 40 19 32 20 31 21 24 26 18 27 17 22 23 28 29 30];
        disp('Selected channel map: Cambridge64_P1_Adpt_Cambridge')
    case 'Cambridge64_H5_Adpt_Cambridge'
        channel_map = [50 14 51 13 52 12 53 7 54 6 55 57 56 1 8 2 58 9 59 15 60 16 61 11 62 10 63 5 64 4 48 3 47 49 46 45 44 43 42 41 25 39 38 37 36 35 34 33 19 20 21 26 27 40 32 31 24 18 17 22 23 28 29 30];
        disp('Selected channel map: Cambridge64_H5_Adpt_Cambridge')
    case 'Cambridge64_H10_Adpt_Cambridge'
        channel_map = [49 12 62 3 13 63 50 14 64 51 7 61 52 6 60 53 57 59 54 1 58 4 2 8 5 9 56 10 15 55 11 16 19 46 28 33 47 29 20 48 30 21 45 23 26 44 22 27 43 17 40 42 18 34 41 24 35 25 31 36 39 32 37 38];
        disp('Selected channel map: Cambridge64_H10_Adpt_Cambridge')
    otherwise
        error('read_openephys: Unknown channel map.')
end  