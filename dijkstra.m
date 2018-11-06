SALVA_VOL = zeros(size(SALVA_BETW));

if strcmpi(ROIs(epi_mean),'Hippocampus') || strcmpi(ROIs(epi_mean),'Accumbens')
    new_vol_accazzo = [epi_mean find(strcmpi(ROIs,'Thalamus')) ...
        find(strcmpi(ROIs,'Temporal')) find(strcmpi(ROIs,'Putamen')) ...
        find(strcmpi(ROIs,'Caudate')) find(strcmpi(ROIs,'Accumbens')) ...
        find(strcmpi(ROIs,'Pallidum')) find(strcmpi(ROIs,'Insula')) ...
        find(strcmpi(ROIs,'Amygdala')) find(strcmpi(ROIs,'Cingulate'))...
        find(strcmpi(ROIs,'Parietal')) find(strcmpi(ROIs,'Frontal'))...
        find(strcmpi(ROIs,'Occipital')) ];
    
    SALVA_VOL(:,1) = epi_mean;
    
    SALVA_VOL(1:(size(SALVA_BETW,1)/3),2) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL((size(SALVA_BETW,1)/3)+1:2*(size(SALVA_BETW,1)/3),2) = find(strcmp(ROIs,'Putamen'));
    SALVA_VOL(2*(size(SALVA_BETW,1)/3):end,2) = find(strcmp(ROIs,'Temporal'));
    SALVA_VOL(1:(size(SALVA_BETW,1)/3),3) = find(strcmp(ROIs,'Putamen'));
    SALVA_VOL((size(SALVA_BETW,1)/3)+1:2*(size(SALVA_BETW,1)/3),3) = find(strcmp(ROIs,'Temporal'));
    SALVA_VOL(2*(size(SALVA_BETW,1)/3):end,3) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(1:(size(SALVA_BETW,1)/3),4) = find(strcmp(ROIs,'Temporal'));
    SALVA_VOL((size(SALVA_BETW,1)/3)+1:2*(size(SALVA_BETW,1)/3),4) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(2*(size(SALVA_BETW,1)/3):end,4) = find(strcmp(ROIs,'Putamen'));
    
    SALVA_VOL(1:5,5) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(6:10,5) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(11:15,5) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(16:20,5) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(21:end,5) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(1:5,6) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(6:10,6) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(11:15,6) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(16:20,6) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(21:end,6) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(1:5,7) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(6:10,7) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(11:15,7) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(16:20,7) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(21:end,7) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(1:5,8) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(6:10,8) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(11:15,8) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(16:20,8) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(21:end,8) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(1:5,9) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(6:10,9) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(11:15,9) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(16:20,9) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(21:end,9) = find(strcmp(ROIs,'Caudate'));
    
    SALVA_VOL(1:(size(SALVA_BETW,1)/2),10) = find(strcmp(ROIs,'Cingulate'));
    SALVA_VOL((size(SALVA_BETW,1)/2)+1:end,10) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL(1:(size(SALVA_BETW,1)/2),11) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL((size(SALVA_BETW,1)/2)+1:end,11) = find(strcmp(ROIs,'Cingulate'));
    
    SALVA_VOL(1:(size(SALVA_BETW,1)/2),12) = find(strcmp(ROIs,'Frontal'));
    SALVA_VOL((size(SALVA_BETW,1)/2)+1:end,12) = find(strcmp(ROIs,'Occipital'));
    SALVA_VOL(1:(size(SALVA_BETW,1)/2),13) = find(strcmp(ROIs,'Occipital'));
    SALVA_VOL((size(SALVA_BETW,1)/2)+1:end,13) = find(strcmp(ROIs,'Frontal'));
    
elseif strcmpi(ROIs(epi_mean),'Insula')
    new_vol_accazzo = [epi_mean find(strcmp(ROIs,'Temporal')) ...
        find(strcmp(ROIs,'Putamen')) find(strcmp(ROIs,'Caudate'))...
        find(strcmp(ROIs,'Accumbens')) find(strcmp(ROIs,'Hippocampus')) ...
        find(strcmp(ROIs,'Amygdala')) find(strcmp(ROIs,'Thalamus')) ...
        find(strcmp(ROIs,'Pallidum')) find(strcmp(ROIs,'Frontal')) ...
        find(strcmp(ROIs,'Parietal'))  find(strcmp(ROIs,'Occipital'))  ...
        find(strcmp(ROIs,'Cingulate')) ];
    
    SALVA_VOL(:,1) = epi_mean;
    
    SALVA_VOL(1:(size(SALVA_BETW,1)/2),2) = find(strcmp(ROIs,'Temporal'));
    SALVA_VOL((size(SALVA_BETW,1)/2)+1:end,2) = find(strcmp(ROIs,'Putamen'));
    SALVA_VOL(1:(size(SALVA_BETW,1)/2),3) = find(strcmp(ROIs,'Putamen'));
    SALVA_VOL((size(SALVA_BETW,1)/2)+1:end,3) = find(strcmp(ROIs,'Temporal'));
    
    SALVA_VOL(1:4,4) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(5:8,4) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(9:12,4) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(13:16,4) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(17:20,4) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(21:end,4) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(1:4,5) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(5:8,5) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(9:12,5) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(13:16,5) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(17:20,5) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(21:end,5) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(1:4,6) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(5:8,6) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(9:12,6) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(13:16,6) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(17:20,6) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(21:end,6) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(1:4,7) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(5:8,7) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(9:12,7) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(13:16,7) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(17:20,7) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(21:end,7) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(1:4,8) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(5:8,8) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(9:12,8) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(13:16,8) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(17:20,8) = find(strcmp(ROIs,'Caudate'));
    SALVA_VOL(21:end,8) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(1:4,9) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(5:8,9) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(9:12,9) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(13:16,9) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(17:20,9) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(21:end,9) = find(strcmp(ROIs,'Caudate'));
    
    SALVA_VOL(1:(size(SALVA_BETW,1)/3),10) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL((size(SALVA_BETW,1)/3)+1:2*(size(SALVA_BETW,1)/3),10) = find(strcmp(ROIs,'Occipital'));
    SALVA_VOL(2*(size(SALVA_BETW,1)/3):end,10) = find(strcmp(ROIs,'Frontal'));
    SALVA_VOL(1:(size(SALVA_BETW,1)/3),11) = find(strcmp(ROIs,'Occipital'));
    SALVA_VOL((size(SALVA_BETW,1)/3)+1:2*(size(SALVA_BETW,1)/3),11) = find(strcmp(ROIs,'Frontal'));
    SALVA_VOL(2*(size(SALVA_BETW,1)/3):end,11) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL(1:(size(SALVA_BETW,1)/3),12) = find(strcmp(ROIs,'Frontal'));
    SALVA_VOL((size(SALVA_BETW,1)/3)+1:2*(size(SALVA_BETW,1)/3),12) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL(2*(size(SALVA_BETW,1)/3):end,12) = find(strcmp(ROIs,'Occipital'));
    
    SALVA_VOL(:,13) = find(strcmp(ROIs,'Cingulate'));
    
    
elseif strcmpi(ROIs(epi_mean),'Caudate') % for P2
    new_vol_accazzo = [epi_mean find(strcmp(ROIs,'Thalamus')) ...
        find(strcmp(ROIs,'Accumbens')) find(strcmp(ROIs,'Pallidum'))...
        find(strcmp(ROIs,'Putamen')) find(strcmp(ROIs,'Amygdala')) ...
        find(strcmp(ROIs,'Hippocampus')) find(strcmp(ROIs,'Insula')) ...
        find(strcmp(ROIs,'Temporal')) find(strcmp(ROIs,'Cingulate')) ...
        find(strcmp(ROIs,'Occipital'))  find(strcmp(ROIs,'Frontal'))  ...
        find(strcmp(ROIs,'Parietal')) ];
    
    SALVA_VOL(:,1) = epi_mean;
    
    SALVA_VOL(1:6,2) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(7:12,2) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(13:20,2) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(21:end,2) = find(strcmp(ROIs,'Putamen'));
    SALVA_VOL(1:6,3) = find(strcmp(ROIs,'Putamen'));
    SALVA_VOL(7:12,3) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(13:20,3) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(21:end,3) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(1:6,4) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(7:12,4) = find(strcmp(ROIs,'Putamen'));
    SALVA_VOL(13:20,4) = find(strcmp(ROIs,'Thalamus'));
    SALVA_VOL(21:end,4) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(1:6,5) = find(strcmp(ROIs,'Accumbens'));
    SALVA_VOL(7:12,5) = find(strcmp(ROIs,'Pallidum'));
    SALVA_VOL(13:20,5) = find(strcmp(ROIs,'Putamen'));
    SALVA_VOL(21:end,5) = find(strcmp(ROIs,'Thalamus'));
    
    SALVA_VOL(1:6,6) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(7:12,6) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(13:20,6) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(21:end,6) = find(strcmp(ROIs,'Temporal'));
    SALVA_VOL(1:6,7) = find(strcmp(ROIs,'Temporal'));
    SALVA_VOL(7:12,7) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(13:20,7) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(21:end,7) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(1:6,8) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(7:12,8) = find(strcmp(ROIs,'Temporal'));
    SALVA_VOL(13:20,8) = find(strcmp(ROIs,'Amygdala'));
    SALVA_VOL(21:end,8) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(1:6,9) = find(strcmp(ROIs,'Hippocampus'));
    SALVA_VOL(7:12,9) = find(strcmp(ROIs,'Insula'));
    SALVA_VOL(13:20,9) = find(strcmp(ROIs,'Temporal'));
    SALVA_VOL(21:end,9) = find(strcmp(ROIs,'Amygdala'));
    
    SALVA_VOL(1:6,10) = find(strcmp(ROIs,'Cingulate'));
    SALVA_VOL(7:12,10) = find(strcmp(ROIs,'Occipital'));
    SALVA_VOL(13:20,10) = find(strcmp(ROIs,'Frontal'));
    SALVA_VOL(21:end,10) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL(1:6,11) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL(7:12,11) = find(strcmp(ROIs,'Cingulate'));
    SALVA_VOL(13:20,11) = find(strcmp(ROIs,'Occipital'));
    SALVA_VOL(21:end,11) = find(strcmp(ROIs,'Frontal'));
    SALVA_VOL(1:6,12) = find(strcmp(ROIs,'Frontal'));
    SALVA_VOL(7:12,12) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL(13:20,12) = find(strcmp(ROIs,'Cingulate'));
    SALVA_VOL(21:end,12) = find(strcmp(ROIs,'Occipital'));
    SALVA_VOL(1:6,13) = find(strcmp(ROIs,'Occipital'));
    SALVA_VOL(7:12,13) = find(strcmp(ROIs,'Frontal'));
    SALVA_VOL(13:20,13) = find(strcmp(ROIs,'Parietal'));
    SALVA_VOL(21:end,13) = find(strcmp(ROIs,'Cingulate'));

    
else
    SALVA_VOL = SALVA_SHORT;
    %     new_vol_accazzo = new_short_accazzo;
end