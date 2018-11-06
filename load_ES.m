%% create end-stage data only

[RID_unique,last_row]=unique(demo.RID,'last');
[RID_unique_HC,last_row_HC]=unique(demoHC.RID,'last');
[RID_unique_PAT,last_row_PAT]=unique(demoPAT.RID,'last');

demoES = demo(last_row,:);
demoHC_ES = demoHC(last_row_HC,:);
demoPAT_ES = demoPAT(last_row_PAT,:);

zscoresROIs_ES = zscoresROIs(:,last_row);

YY = (mean(zscoresROIs_ES,2) - min(mean(zscoresROIs_ES,2)))./(max(mean(zscoresROIs_ES,2))-min(mean(zscoresROIs_ES,2)));
Y_STD = (std(zscoresROIs_ES')'- min(std(zscoresROIs_ES')'))./(max(std(zscoresROIs_ES')')-min(std(zscoresROIs_ES')'));

table_OUT_es = table(YY,Y_STD);
lunchruns;
for ind = 1:length(runs)
    struct_OUT_ES{ind}.subj = (zscoresROIs_ES(:,ind) - min(zscoresROIs_ES(:,ind)))./(max(zscoresROIs_ES(:,ind))-min(zscoresROIs_ES(:,ind)));; 
end