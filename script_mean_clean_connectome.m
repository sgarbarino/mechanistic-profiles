%% script to clean connectome and keep just the ROIs selected by ROIs selection

connectome_mean = connectome_mean(find(h),find(h));
connectome_std = connectome_std(find(h),find(h));
connectome_label = connectome_label(find(h));
connect = connect(find(h),find(h),:,:);