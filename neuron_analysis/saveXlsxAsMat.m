opts = detectImportOptions('PFC_Chrimson_stGtACR.xlsx');
opts = setvartype(opts, 'Group', 'char');
PFC_Chrimson_stGtACR = readtable('PFC_Chrimson_stGtACR.xlsx', opts);
save ('PFC_Chrimson_stGtACR.mat', 'PFC_Chrimson_stGtACR')