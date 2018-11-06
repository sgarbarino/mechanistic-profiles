GP_csv = array2table(-zscoresROIs');
GP_csv.Properties.VariableNames=ROIs;
GP_csv.RID = demo.RID;
GP_csv.age = demo.AGE;
GP_csv.DX = demo.DX_bl;

writetable(GP_csv, 'GP.csv')