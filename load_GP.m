
%% read files and store them
addpath('npy-matlab-master/')

% separate HC and PAT and their RID (individual ID)

if strcmp(dataset, 'ADNI') || strcmp(dataset, 'dummy')
    % for AD and dummy
    demoFINAL = [demoHC ; demoPAT];
    % add a column, similar to DX_bl, but DX_last, that has the info on the last diagnosis
    [aa,bb]= unique(demoFINAL.RID);
    demoFINAL_RID = demoFINAL(bb,:);
    demoFINAL_RID_PAT = demoFINAL_RID(demoFINAL_RID.DX_bl~='CN',:);
elseif strcmp(dataset,'PPMS')
    % for MS
    demoFINAL = [demoHC ; demoPAT];
    % % add a column, similar to DX_bl, but DX_last, that has the info on the last diagnosis
    [aa,bb]= unique(demoFINAL.ID2);
    demoFINAL_RID = demoFINAL(bb,:);
    demoFINAL_RID_PAT = demoFINAL_RID(demoFINAL_RID.subtype~='PPMS',:);
%% HA is already OK because there is no distinction between HC/Patients
end

%% load RIDs 
RID = unique(Data_tot.RID);
% if just PAT
% RID = unique(demoPAT.RID);

% if just healthy sub
% RID = unique(demoHC.RID);

load data/ROIs.mat
names_biomarkers = ROIs;

N_RUN = length(RID); % number or runs - of individual subjects
N_biom = size(names_biomarkers,2); % number of biomarkers
mydataSTD = cell(N_biom,N_RUN);
mydataX = cell(N_biom,1);
mydataY = cell(N_biom,1);

for kbiom = 0:N_biom-1
    for k = 0:N_RUN-1
        myfilenameSTD = sprintf('results-GP/xyaxis/stdyaxis_%d_%d.npy', kbiom,k);
        mydataSTD{kbiom+1,k+1} = readNPY(myfilenameSTD);
    end
    myfilenameX = sprintf('results-GP/xyaxis/xaxis%d.npy', kbiom);

    mydataX{kbiom+1} = readNPY(myfilenameX);
    myfilenameY = sprintf('results-GP/xyaxis/yaxis%d.npy', kbiom);
    mydataY{kbiom+1} = readNPY(myfilenameY);
    disp(kbiom)
end
X = (cell2mat(mydataX'))'+30;
Y = (cell2mat(mydataY'))';
N_points = length(mydataX{1}); % number of points on Xaxis

% compute std for every point across all the lines
M_Max = zeros(N_biom,N_points);
m_min = zeros(N_biom,N_points);

for kbiom = 1: N_biom
    STD_tmp = mydataSTD(kbiom,:);
    STD = cell2mat(STD_tmp);
    M_Max(kbiom,:) = max(STD');
    m_min(kbiom,:) = min(STD');
end

%% load also measurements 

mydataXDOTS = cell(length(RID),N_biom); % cells with columns = biomarkers and rows = subject measurements
mydataYDOTS = cell(length(RID),N_biom);

for kbiom = 0:N_biom-1
    for kd = 0:length(RID)-1
        myfilenameXDOTS = sprintf('results-GP/xyaxis/xaxis_dots%d_%d.npy', kd, kbiom);

        mydataXDOTS{kd+1,kbiom+1} = readNPY(myfilenameXDOTS);
        myfilenameYDOTS = sprintf('results-GP/xyaxis/yaxis_dots%d_%d.npy', kd, kbiom);

        mydataYDOTS{kd+1,kbiom+1} = readNPY(myfilenameYDOTS);
    end
    disp(kbiom)
end

%% derivative

ALL_CURVES = zeros(N_biom,N_points,N_RUN);
DER_all = zeros(N_biom,N_points-1,N_RUN); % they are one dimension shorter on x-axis

for kbiom = 1:N_biom
    ALL_CURVES(kbiom,:,:)=cell2mat(mydataSTD(kbiom,:));
end
for kbiom = 1:N_biom
    for krun = 1:N_RUN
        DER_all(kbiom,:,krun)=(diff(ALL_CURVES(kbiom,:,krun)))./diff(X(kbiom,:));
    end
end     

DER_mean = zeros(N_biom,N_points-1); % they are one dimension shorter on x-axis
for kbiom = 1:N_biom
    DER_mean(kbiom,:)=(diff(Y(kbiom,:)))./diff(X(kbiom,:));
end

%******* compute max derivative 
val = zeros(N_biom,1);
pos = zeros(N_biom,1);

for kbiom = 1:N_biom
    [val(kbiom),pos(kbiom)] = max(DER_mean(kbiom,:));
end
% and order ROIs accordingly
table_OUT = table(val, pos, [1:N_biom]',names_biomarkers', std(DER_mean')');
table_OUT.Properties.VariableNames = {'der_values','der_position','biom_position','biom_names', 'der_values_std'}';
table_OUT_sorted=sortrows(table_OUT,'der_values','descend');
% normalization
table_OUT.der_values = (table_OUT.der_values- min(table_OUT.der_values))./(max(table_OUT.der_values)-min(table_OUT.der_values));

% ******** compute max derivative also of all the subjs (repeat above steps, just bigger table)

val_all = zeros(N_biom,N_RUN);
pos_all = zeros(N_biom,N_RUN);

for kbiom = 1:N_biom
    for krun = 1:N_RUN
        [val_all(kbiom,krun),pos_all(kbiom,krun)] = max(DER_all(kbiom,:,krun));
    end
end

lunchruns; % load names of subjects (i.e. their number) - different for different dataset

for ind = 1:length(runs)
    struct_OUT_all.(runs{ind}) = table(val_all(:,ind),pos_all(:,ind),[1:N_biom]',names_biomarkers');
    struct_OUT_all.(runs{ind}).Properties.VariableNames = {'der_values','der_position','biom_position','biom_names'}';
    struct_OUT_all_sorted.(runs{ind}) = sortrows(struct_OUT_all.(runs{ind}),'der_values','descend');
end
% normalization
for ind = 1:length(runs)
    struct_OUT_all.(runs{ind}).der_values = (struct_OUT_all.(runs{ind}).der_values- min(struct_OUT_all.(runs{ind}).der_values))./(max(struct_OUT_all.(runs{ind}).der_values)-min(struct_OUT_all.(runs{ind}).der_values));

end
struct_OUT_all_sorted.('run0') = table_OUT_sorted;

%% ordering

% table_OUT_sorted_PURE_alltempt_biom_pos = zeros(N_biom,length(runs));
% for ind = 1:length(runs)
%     table_OUT_sorted_PURE_alltempt_biom_pos(:,ind) = struct_OUT_all_sorted.(runs{ind}).biom_position;
% end
% 
% for i=1:5
%     table_OUT_sorted_PURE_alltempt_biom_pos(:,end+1) = struct_OUT_all_sorted.('run0').biom_position; % aggiungo il run0
% end
% table_OUT_sorted_PURE_alltempt_biom_pos = table_OUT_sorted_PURE_alltempt_biom_pos';

%% all figures at the end

%********** figure all biom 
figure
grid on

for kbiom = [1:N_biom] 
    x = X(kbiom,:)'; y = Y(kbiom,:)';
    hold on; plot(x,y,'linewidth',2)
end

legend(names_biomarkers)
xlabel('Re-parametrized age ','Fontsize',20)
ylabel('Abnormality','Fontsize',20)

%% fig derivative

figure, grid on
for kbiom=1:N_biom
    hold on
    plot(X(kbiom,1:end-1),(DER_mean(kbiom,:)),'linewidth',2)
end

legend(names_biomarkers)
xlabel('Re-parametrized age ','Fontsize',20)
ylabel('Derivatives','Fontsize',20)

%% %*********** now again but all separate plots
figure
grid on
for kbiom = 1:N_biom
    x = X(kbiom,:)'; y = Y(kbiom,:)';hi = M_Max(kbiom,:)'/1; lo = m_min(kbiom,:)'/1;
            subplot(4,4,kbiom) 
    hold on;
    
        for kd = 1:length(RID)
            % color differently the diagnosis -- AD or dummy (+30)
            if demoFINAL_RID.DX_bl(kd)=='CN'
                hold on
                plot(mydataXDOTS{kd,kbiom}(1)+30,mydataYDOTS{kd,kbiom}(1),'*g') % adding (1) I select just BL
            elseif demoFINAL_RID.DX_bl(kd)=='EMCI'
                plot(mydataXDOTS{kd,kbiom}(1)+30,mydataYDOTS{kd,kbiom}(1),'*m')
            elseif demoFINAL_RID.DX_bl(kd)=='LMCI'
                plot(mydataXDOTS{kd,kbiom}+30,mydataYDOTS{kd,kbiom},'*m')
            elseif demoFINAL_RID.DX_bl(kd)=='AD'
                plot(mydataXDOTS{kd,kbiom}(1)+30,mydataYDOTS{kd,kbiom}(1),'*k')        
            end
           % color differently the age -- RS
%             if Data_demo.age(kd)<65
%                 hold on
%                 plot(mydataXDOTS{kd,kbiom},mydataYDOTS{kd,kbiom},'.m')
%             elseif Data_demo.age(kd)>85
%                 plot(mydataXDOTS{kd,kbiom},mydataYDOTS{kd,kbiom},'.k')  
%             else
%                      plot(mydataXDOTS{kd,kbiom},mydataYDOTS{kd,kbiom},'.g')
%             end

%            %color differently the diagnosis -- MS %%%%%% something wrong
%             if demoFINAL_RID.subtype(kd)=='HC'
%                 hold on
%                 plot(mydataXDOTS{kd,kbiom}(1),mydataYDOTS{kd,kbiom},'*m') % adding (1) I select just BL
%             elseif demoFINAL_RID.subtype(kd)=='PPMS'
%                 plot(mydataXDOTS{kd,kbiom}(1),mydataYDOTS{kd,kbiom},'*g')
%             end
% 
        end                                         
    
    plot(x,y+hi/2,'-.r'); plot(x,y-hi/2,'-.r')
    plot(x,y,'linewidth',2) % plot traiettorie
    hl = line(x,y);
    set(hl, 'color', 'r','linewidth',2);%, 'marker', 'x');
    xlabel('Re-parametrized Age','Fontsize',10)
    ylabel('Abnormality','Fontsize',10)
%     legend(hl,names_biomarkers{kbiom})
        axis([50 100 -1 1]) % ADNI


    axis square
    title(names_biomarkers{kbiom})
end