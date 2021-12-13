%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPARISON DALELIO FN PCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
clear
clc

%% DATA AGGREGATOR 
val.Nsupply = [5];%[3,10];
val.tau     = [0.2 0.4 0.6 0.8 1 1/0.8 1/0.6 1/0.4 1/0.2];%[0.2, 0.4, 0.6, 0.8, 1, 1/0.8, 1/0.6, 1/0.4, 1/0.2];
val.T       = [20];%[5, 10, 15, 20, 25, 30];
val.gamma   = [0.2:0.15:1];
val.tempfuA = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.tempfuZ = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.time    = [1000];%[1, 2, 3];

for tau = 1:length(val.tau)
    for light = 1:length(val.gamma)
        path =  ['Output/COMPARING_T' num2str(val.T) '/Run_CA_N5_Tau' num2str(val.tau(tau)) '_EiE_T' num2str(val.T()) '_tfuA' num2str(val.tempfuA) '_tfuZ' num2str(val.tempfuZ) '_gamma' num2str(val.gamma(light)) '/Network' num2str(val.time) '.mat'];
        load(path)
        prdpry = eco.p.prdpry;
        P = eco.P;
        FN = P;
        Nutri = yout(53,1);
        yout = reshape(yout(53,2:end),[sqrt(length(P)),sqrt(length(P))]);

        ind = FN < 0;
        FN(ind) = 0;
        clear ind

        for u = 1:size(prdpry,1)
            for m = 1:length(P)
                pry_P(u,m) = prdpry(u,m)*P(m);
            end
        end
        GGain  = eco.g_refg .* pry_P;
        EDGE = GGain.*eco.p.lambda;
        EDGE = full(EDGE);
        EDGE = cat(2,EDGE,eco.Nuptk.*P);

        FNS{tau,light} = FN;
        EDGES{tau,light} = EDGE;

        clear u refgiu p
    end
end

clearvars -except EDGES FNS val exct

%% PERCENTAGES
k = 1;
FNcol = zeros(length(FNS{1,1}),1);
for i = 1:size(FNS,1)
    for j = 1:size(FNS,2)
        FNcol(:,k) = FNS{i,j};
        k = k +1;
    end
end

xls = xlsread('nodecompa.xlsx');
FNcol(:,end+1) = zeros(length(FNcol),1);
FNcol(xls(:,3),end) = xls(:,4);            %CHANGE xls(:,4) FOR BLUE AND xls(:,5) FOR GREEN

FNcol(FNcol < exp(-10)) = 0;
FNcol = cat(1,FNcol,sum(FNcol));

for i = 1:size(FNcol,1)-1
    for j = 1:size(FNcol,2)
        FNperc(i,j) = FNcol(i,j)/FNcol(end,j)*100;
    end
end

clear i j k

%% PCA
pcaFN = pca(FNperc);
%% UFF--------------------------------------------------------------
for i = 1:size(FNperc-1,2)
    pytha(i) = sqrt(sum((FNperc(:,end) - FNperc(:,i)) .^ 2));
end

%% SCATTER
scatcolor = zeros(length(pcaFN),3);
scatcolor(end,:) = [1,0,0];
scatter(pcaFN(:,1),pcaFN(:,2),[],scatcolor)

clear scatcolor

%% HEATMAP
for i = 1:size(pcaFN,1)-1
    pytha(i) = sqrt(sum((pcaFN(end,1:2) - pcaFN(i,1:2)) .^ 2));
end

pythareshape = reshape(pytha,[6,9]);
h = heatmap(pythareshape');
h.Title = ['Comparison'];
h.XDisplayLabels = ['0.20'; '0.35'; '0.50'; '0.65'; '0.80'; '0.95'];
h.YDisplayLabels = ['0.20'; '0.40'; '0.60'; '0.80' ;'1.00' ;'1.25' ;'1.66' ;'2.50' ;'5.00'];
h.XLabel = 'gamma';
h.YLabel = 'tau';

clear i

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPARISON DALELIO FN MULTIDIMENSIONAL EUCLIDIAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
clear
clc

%% DATA AGGREGATOR 
val.Nsupply = [5];%[3,10];
val.tau     = [0.2 0.4 0.6 0.8 1 1/0.8 1/0.6 1/0.4 1/0.2];%[0.2, 0.4, 0.6, 0.8, 1, 1/0.8, 1/0.6, 1/0.4, 1/0.2];
val.T       = [25];%[5, 10, 15, 20, 25, 30];
val.gamma   = [0.2:0.15:1];
val.tempfuA = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.tempfuZ = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.time    = [1000];%[1, 2, 3];

for tau = 1:length(val.tau)
    for light = 1:length(val.gamma)
        path =  ['Output/COMPARING_T' num2str(val.T) '/Run_CA_N5_Tau' num2str(val.tau(tau)) '_EiE_T' num2str(val.T()) '_tfuA' num2str(val.tempfuA) '_tfuZ' num2str(val.tempfuZ) '_gamma' num2str(val.gamma(light)) '/Network' num2str(val.time) '.mat'];
        load(path)
        prdpry = eco.p.prdpry;
        P = eco.P;
        FN = P;
        Nutri = yout(53,1);
        yout = reshape(yout(53,2:end),[sqrt(length(P)),sqrt(length(P))]);

        ind = FN < 0;
        FN(ind) = 0;
        clear ind

        for u = 1:size(prdpry,1)
            for m = 1:length(P)
                pry_P(u,m) = prdpry(u,m)*P(m);
            end
        end
        GGain  = eco.g_refg .* pry_P;
        EDGE = GGain.*eco.p.lambda;
        EDGE = full(EDGE);
        EDGE = cat(2,EDGE,eco.Nuptk.*P);

        FNS{tau,light} = FN;
        EDGES{tau,light} = EDGE;

        clear u refgiu p
    end
end

clearvars -except EDGES FNS val exct

%% PERCENTAGES
k = 1;
FNcol = zeros(length(FNS{1,1}),1);
for i = 1:size(FNS,1)
    for j = 1:size(FNS,2)
        FNcol(:,k) = FNS{i,j};
        k = k +1;
    end
end

xls = xlsread('nodecompa.xlsx');
FNcol(:,end+1) = zeros(length(FNcol),1);
FNcol(xls(:,3),end) = xls(:,4);            %CHANGE xls(:,4) FOR BLUE AND xls(:,5) FOR GREEN

FNcol(FNcol < exp(-10)) = 0;
FNcol = cat(1,FNcol,sum(FNcol));

for i = 1:size(FNcol,1)-1
    for j = 1:size(FNcol,2)
        FNperc(i,j) = FNcol(i,j)/FNcol(end,j)*100;
    end
end

clear i j k


%% MULTIDIMENSIONAL EUCLIDEAN DISTANCE

pytha = FNperc(:,end);
pythasum = zeros(size(FNperc,1),1);
pythares = zeros(size(FNperc,2)-1,1);
for i = 1:size(FNperc,2)-1
    for j = 1:size(FNperc,1)
        pythasum(j) = (pytha(j,end) - FNperc(j,i)) .^ 2;
    end
    pythares(i) = sum(pythasum);
end
pytha = sqrt(pythares);
        
%% HEATMAP
pythareshape = reshape(pytha,[6,9]);
h = heatmap(pythareshape');
h.Title = ['Comparison'];
h.XDisplayLabels = ['0.20'; '0.35'; '0.50'; '0.65'; '0.80'; '0.95'];
h.YDisplayLabels = ['0.20'; '0.40'; '0.60'; '0.80' ;'1.00' ;'1.25' ;'1.66' ;'2.50' ;'5.00'];
h.XLabel = 'gamma';
h.YLabel = 'tau';

clear i


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLUE-RED MAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

%% DATA AGGREGATOR 

val.Nsupply = [5];%[3,10];
val.tau     = [1];%[0.2, 0.4, 0.6, 0.8, 1, 1/0.8, 1/0.6, 1/0.4, 1/0.2];%[0.8, 1, 1/0.8, 1/0.6, 1/0.4];%
val.T       = [5, 10, 15, 20, 25, 30];%20;%[5, 10, 15, 20, 25, 30];
val.gamma   = [0.2:0.15:1];
val.tempfuA = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.tempfuZ = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.time    = [3];%[1, 2, 3, 10, 20, 100, 200, 700, 1000];

tau = 1;
totexp = 1;
for temp = 1:length(val.T)
    for light = 1:length(val.gamma)
        
        
path =  ['Output/Run_CA_N5_Tau' num2str(val.tau(tau)) '_EiE_T' num2str(val.T(temp)) '_tfuA' num2str(val.tempfuA) '_tfuZ' num2str(val.tempfuZ) '_gamma' num2str(val.gamma(light)) '/Network' num2str(val.time) '.mat'];
load(path)
prdpry = eco.p.prdpry;
P = eco.P;
FN = P;
Nutri = yout(53,1);
yout = reshape(yout(53,2:end),[sqrt(length(P)),sqrt(length(P))]);

val.sizes   = [length(yout)];

ind = FN < 0;
FN(ind) = 10^-40;
clear ind

P;
for u = 1:size(prdpry,1)
    for m = 1:length(P)
        pry_P(u,m) = prdpry(u,m)*P(m);
    end
end
GGain  = eco.g_refg .* pry_P;
EDGE = P.*GGain.*eco.p.lambda;
EDGE = full(EDGE);

EDGE = cat(2,EDGE,eco.Nuptk.*P);

FNS{temp,light} = FN;
EDGES{temp,light} = EDGE;

clear u refgiu p
    end
end

clearvars -except EDGES FNS val


%% BLUE-RED COLOR MAP
map(1,:) = [0:0.01:1];
map(2,:) = [0:0.01:1];
map(3,:) = [0:0.01:1];
map(1,:) = ones(1,length(map));
map = cat(2,map,flip(flip(map(:,1:end-1),1),2));
map=map';


%% GAMMA-TEMPERATURE COMPARISON
%% HEATMAPS
% 1 = Temperature; 2 = Gamma; 3 = tau;
STATtoCOMPARE =2;

for k = 1:length(val.T)

    if STATtoCOMPARE == 1
        %COMPARE TEMPERATURE
        rowcomp1 = k;
        colcomp1 = 1;
        rowcomp2 = k;
        colcomp2 = size(EDGES,2);
    end

    if STATtoCOMPARE == 2
        %COMPARE GAMMA
        rowcomp1 = 1;
        colcomp1 = k;
        rowcomp2 = size(EDGES,1);
        colcomp2 = k;
    end

    for i = 1:size(EDGES{1,1},1)
        for j = 1:size(EDGES{1,1},2)
            EDGEHEAT(i,j) = (EDGES{rowcomp1,colcomp1}(i,j)-EDGES{rowcomp2,colcomp2}(i,j))/(EDGES{rowcomp1,colcomp1}(i,j)+EDGES{rowcomp2,colcomp2}(i,j));
        end
    end
    EDGEHEAT(isinf(EDGEHEAT)|isnan(EDGEHEAT)) = 0;

    h = heatmap(EDGEHEAT,'Colormap',map);
    h.XLabel = 'Prey No.';
    h.YLabel = 'Predator No.';
    
    if STATtoCOMPARE == 1
        h.Title = {'tau 1','BLUE = T5 bigger','RED = T30 bigger'};
    end
    if STATtoCOMPARE == 2
        h.Title = {'tau 1','BLUE = gamma0.2 bigger','RED = gamma0.95 bigger'};
    end
    
    set(gcf, 'Position', [200 40 750 750]);
    set(gcf, 'Color', 'w');
    
    if STATtoCOMPARE == 1
        fname=['Output/Run_CA_N5_Tau' num2str(val.tau()) '_EiE_T5-30_gamma' num2str(val.gamma(k)) '.png'];
        saveas(gcf,fname);
    end
    if STATtoCOMPARE == 2
        fname=['Output/Run_CA_N5_Tau' num2str(val.tau()) '_EiE_T' num2str(val.T(k)) '_gamma0.2-0.95.png'];
        saveas(gcf,fname);
    end

end

clear i j k rowcomp1 colcomp1 rowcomp2 colcomp2 STATtoCOMPARE fnames




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTINCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear
clc

%% KNOCKOUT COMPARISON

val.Nsupply = [5];%[3,10];
val.tau     = [1];%[0.2, 0.4, 0.6, 0.8, 1, 1/0.8, 1/0.6, 1/0.4, 1/0.2];
val.T       = [20];%20;%[5, 10, 15, 20, 25, 30];
val.gamma   = [0.1];%[0.2:0.15:0.5];
val.tempfuA = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.tempfuZ = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.time    = [1];%[1, 2, 3];

exct = {'AUT' 'AUT-HETER' 'AUT-MIX' 'HETER' 'HETER-MIX' 'MIX'};

totexp = 1;
for e = 1:length(exct)
    for time = 1:length(val.time)
        path =  ['Output/STRESS_Knock-Out/Run_CP_N5_Tau1_EiE_' exct{e} '_T20_tfuA0.05_tfuZ0.05_gamma0.1/Network' num2str(val.time) '.mat'];
        load(path)
        prdpry = eco.p.prdpry;
        P = eco.P;
        FN = P;
        Nutri = yout(53,1);
        yout = reshape(yout(53,2:end),[sqrt(length(P)),sqrt(length(P))]);

        val.sizes   = [length(yout)];


        ind = FN < 0;
        FN(ind) = 10^-40;
        clear ind

        P;
        for u = 1:size(prdpry,1)
            for m = 1:length(P)
                pry_P(u,m) = prdpry(u,m)*P(m);
            end
        end
        GGain  = eco.g_refg .* pry_P;
        EDGE = P.*GGain.*eco.p.lambda;
        EDGE = full(EDGE);

        EDGE = cat(2,EDGE,eco.Nuptk.*P);

        FNS{e,time} = FN;
        EDGES{e,time} = EDGE;

        clear u refgiu p
    end
end

clearvars -except EDGES FNS val exct


path =  ['Output/STRESS_Knock-Out/ssRun_CA_N5_Tau1_EiE_T20_tfuA0.05_tfuZ0.05_gamma0.1/Network3.mat'];
load(path)
prdpry = eco.p.prdpry;
P = eco.P;
FN = P;
Nutri = yout(53,1);
yout = reshape(yout(53,2:end),[sqrt(length(P)),sqrt(length(P))]);

val.sizes   = [length(yout)];

ind = FN < 0;
FN(ind) = 10^-40;
clear ind

P;
for u = 1:size(prdpry,1)
    for m = 1:length(P)
        pry_P(u,m) = prdpry(u,m)*P(m);
    end
end
GGain  = eco.g_refg .* pry_P;
EDGE = P.*GGain.*eco.p.lambda;
EDGE = full(EDGE);

EDGE = cat(2,EDGE,eco.Nuptk.*P);

FNS{7,1} = FN;
EDGES{7,1} = EDGE;

clearvars -except EDGES FNS val exct

%% BLUE-RED COLOR MAP
map(1,:) = [0:0.01:1];
map(2,:) = [0:0.01:1];
map(3,:) = [0:0.01:1];
map(1,:) = ones(1,length(map));
map = cat(2,map,flip(flip(map(:,1:end-1),1),2));
map=map';

%% ELIMINATE AUTO-HETERO
INDX = reshape(1:(val.sizes*val.sizes),val.sizes,val.sizes)';
AUTO = sort(reshape(INDX(:,1),1,[]));
MIXO = sort(reshape(INDX(:,2:(val.sizes-1)),1,[]));
HETE = sort(reshape(INDX(:,val.sizes),1,[]));

%% COMPARISON

for k = 1:length(EDGES)-1
for i = 1:size(EDGES{1,1},1)
    for j = 1:size(EDGES{1,1},2)
        EDGEHEAT(i,j) = (EDGES{k}(i,j)-EDGES{length(EDGES)}(i,j))/(EDGES{k}(i,j)+EDGES{length(EDGES)}(i,j));
    end
end
EDGEHEAT(isinf(EDGEHEAT)|isnan(EDGEHEAT)) = 0;
EDGEHEAT(size(EDGEHEAT,1)-1,1) = 1;
EDGEHEAT(size(EDGEHEAT,1)-1,2) = -1;
EDGEHEATcell{k} = EDGEHEAT;

h = heatmap(EDGEHEAT(MIXO,:),'Colormap',map);
h.XLabel = 'Prey No.';
h.YLabel = 'Predator No.';
h.YData = MIXO;
h.Title = {'10 years after the exctintion event','BLUE = Gain higher after extinction','RED = Gain higher before extinction'};
set(gcf, 'Position', [200 40 750 750]);
set(gcf, 'Color', 'w');
fname=['Output/Run_CP_N5_Tau1_EiE_' exct{k} '_T20_gamma0.1' num2str(val.time) '.png'];
saveas(gcf,fname);

end








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DALELIO ZONE COMPARISON
%% PARTITIONING
INDX = reshape(1:(val.sizes*val.sizes),val.sizes,val.sizes)';

AUTO = sort(reshape(INDX(:,1),1,[]));
MIXO = sort(reshape(INDX(:,2:(val.sizes-1)),1,[]));
HETE = sort(reshape(INDX(:,val.sizes),1,[]));

SMAL = sort(reshape(INDX(1:4,:),1,[]));
MEDI = sort(reshape(INDX(5:6,:),1,[]));
BIG  = sort(reshape(INDX(7:val.sizes,:),1,[]));

memb = ismember(AUTO,SMAL); AUTO_SMAL = AUTO(memb);
memb = ismember(AUTO,MEDI); AUTO_MEDI = AUTO(memb);
memb = ismember(AUTO,BIG); AUTO_BIG = AUTO(memb);

memb = ismember(MIXO,SMAL); MIXO_SMAL = MIXO(memb);
memb = ismember(MIXO,MEDI); MIXO_MEDI = MIXO(memb);
memb = ismember(MIXO,BIG); MIXO_BIG = MIXO(memb);

memb = ismember(HETE,SMAL); HETE_SMAL = HETE(memb);
memb = ismember(HETE,MEDI); HETE_MEDI = HETE(memb);
memb = ismember(HETE,BIG); HETE_BIG = HETE(memb);

columnname = {'SMALL','MEDIUM','BIG'};
rowname = {'STRATEGY','SIZE','AUTO SIZE','MIXO SIZE','HETERO SIZE','/'}';
CATEG = {AUTO,MIXO,HETE;SMAL,MEDI,BIG;AUTO_SMAL,AUTO_MEDI,AUTO_BIG;MIXO_SMAL,MIXO_MEDI,MIXO_BIG;HETE_SMAL,HETE_MEDI,HETE_BIG};
CATEG = cat(1,CATEG,columnname);
CATEG = cat(2,CATEG,rowname);

clear columnname rowname AUTO MIXO HETE MEDI SMAL BIG
clear memb AUTO_SMAL AUTO_MEDI AUTO_BIG MIXO_SMAL MIXO_MEDI MIXO_BIG HETE_SMAL HETE_MEDI HETE_BIG

%% COLLATE

for i = 1:size(EDGES,1)
    for j = 1:size(EDGES,2)
        for k = 1:size(CATEG,1)-1
            for m = 1:size(CATEG,2)-1
                graz{i,j}{k,m} = sum(EDGES{i,j}(CATEG{k,m},1:end-1),2);
            end
        end
    end
end

for i = 1:size(EDGES,1)
    for j = 1:size(EDGES,2)
        for k = 1:size(CATEG,1)-1
            for m = 1:size(CATEG,2)-1
                grazup{i,j}{k,m} = sum(EDGES{i,j}(CATEG{k,m},1:end),2);
            end
        end
    end
end

columnname = {'SMALL','MEDIUM','BIG'};
rowname = {'STRATEGY','SIZE','AUTO SIZE','MIXO SIZE','HETERO SIZE','/'}';

for k = 1:size(CATEG,1)-1
    for m = 1:size(CATEG,2)-1
        graz{k,m} = cat(1,graz{k,m},columnname);
        graz{k,m} = cat(2,graz{k,m},rowname);
    end
end

for k = 1:size(CATEG,1)-1
    for m = 1:size(CATEG,2)-1
        grazup{k,m} = cat(1,grazup{k,m},columnname);
        grazup{k,m} = cat(2,grazup{k,m},rowname);
    end
end

GAIN = cat(2,graz,grazup);

columnname = {'gamma0.2','gamma.350','gamma0.5','gamma0.2','gamma.350','gamma0.5'};
GAIN = cat(1,GAIN,columnname);
columnname = {'Graz','Graz','Graz','Grazup','Grazup','Grazup'};
GAIN = cat(1,GAIN,columnname);
rowname = {'tau0.8','tau1','tau1.25','tau1.6667','tau2.5','/','/'}';
GAIN = cat(2,GAIN,rowname);


%for i = 1:size(EDGES,1)
%    for j = 1:size(EDGES,2)
%        EDGESsumgraz(i,j) = sum(sum(EDGES{i,j}));
%    end
%end

%for i = 1:size(EDGES,1)
%    for j = 1:size(EDGES,2)
%        EDGESsumgrazup(i,j) = sum(sum(EDGES{i,j}(:,1:end-1)));
%    end
%end

%EDGESUM = cat(2,EDGESsumgraz,EDGESsumgrazup);
%columnname = {'gamma0.2','gamma.350','gamma0.5','gamma0.2','gamma.350','gamma0.5'};
%EDGESUM = cat(2,EDGESUM,columnname);
%columnname = {'Graz','Graz','Graz','Grazup','Grazup','Grazup'};
%EDGESUM = cat(2,EDGESUM,columnname);
%rowname = {'tau0.8','tau1','tau1.25','tau1.6667','tau2.5','/','/'}';
%EDGESUM = cat(2,EDGESUM,rowname);

clear columnname rowname i j k m graz grazup EDGESsumgraz EDGESsumgrazup
%% COMPARISON PERCENTAGE

for i = 1:size(GAIN,1)-2
    for j = 1:size(GAIN,2)-1
        for k = 1:size((GAIN{1,1}),1)-1
                sumo = sum(GAIN{i,j}{k,1}) + sum(GAIN{i,j}{k,2}) + sum(GAIN{i,j}{k,3});
            if sumo > 0
                PERC{i,j}(k,1) = sum(GAIN{i,j}{k,1}) / sumo * 100;
                PERC{i,j}(k,2) = sum(GAIN{i,j}{k,2}) / sumo * 100;
                PERC{i,j}(k,3) = sum(GAIN{i,j}{k,3}) / sumo * 100;
            else
                PERC{i,j}(k,:) = [0,0,0];
            end
        end
    end
end

for i = 1:size(GAIN,1)-2
    for j = 1:size(GAIN,2)-1
        for k = 1:size((GAIN{1,1}),1)-1
                sumo = sum(GAIN{i,j}{k,1}) + sum(GAIN{i,j}{k,2}) + sum(GAIN{i,j}{k,3});
            if sumo > 0
                SUMO{i,j}(k,1) = sum(GAIN{i,j}{k,1});
                SUMO{i,j}(k,2) = sum(GAIN{i,j}{k,2});
                SUMO{i,j}(k,3) = sum(GAIN{i,j}{k,3});
            else
                SUMO{i,j}(k,:) = [0,0,0];
            end
        end
    end
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OBSOLETE
% look for the "&&" to see whitch values are modified due to graphic necessities
clear
clc

%% DATA AGGREGATOR 

val.Nsupply = [5];%[3,10];
val.tau     = [0.8, 1, 1/0.8, 1/0.6, 1/0.4];%[0.2, 0.4, 0.6, 0.8, 1, 1/0.8, 1/0.6, 1/0.4, 1/0.2];
val.T       = [20];%20;%[5, 10, 15, 20, 25, 30];
val.gamma   = [0.2:0.15:0.5];
val.tempfuA = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.tempfuZ = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.time    = [3];%[1, 2, 3];

totexp = 1;
for tau = 1:length(val.tau)
    for light = 1:length(val.gamma)
        
        
path =  ['Output/daleliozone/Run_CA_N5_Tau' num2str(val.tau(tau)) '_EiE_T' num2str(val.T()) '_tfuA' num2str(val.tempfuA) '_tfuZ' num2str(val.tempfuZ) '_gamma' num2str(val.gamma(light)) '/Network' num2str(val.time) '.mat'];
load(path)
prdpry = eco.p.prdpry;
P = eco.P;
FN = P;
Nutri = yout(53,1);
yout = reshape(yout(53,2:end),[sqrt(length(P)),sqrt(length(P))]);

val.sizes   = [length(yout)];

ind = FN < 0;
FN(ind) = 10^-40;
clear ind

P;
for u = 1:size(prdpry,1)
    for m = 1:length(P)
        pry_P(u,m) = prdpry(u,m)*P(m);
    end
end
GGain  = eco.g_refg .* pry_P;
EDGE = GGain.*eco.p.lambda;
EDGE = full(EDGE);

EDGE = cat(2,EDGE,eco.Nuptk);

FNS{tau,light} = FN;
EDGES{tau,light} = EDGE;

clear u refgiu p
    end
end

%comp26: temperature on row, light on column
%comp27: tau on row, light on column
%patho = ['comp27.mat'];
%save(patho, 'EDGES', 'FNS')
clearvars -except EDGES FNS val


%% BLUE-RED COLOR MAP
map(1,:) = [0:0.01:1];
map(2,:) = [0:0.01:1];
map(3,:) = [0:0.01:1];
map(1,:) = ones(1,length(map));
map = cat(2,map,flip(flip(map(:,1:end-1),1),2));
map=map';


%% GAMMA-TEMPERATURE COMPARISON
%% HEATMAPS

for k = 1:length(val.T)
%rowcomp1 = k;
%colcomp1 = 1;
%rowcomp2 = k;
%colcomp2 = 6;

rowcomp1 = 1;
colcomp1 = k;
rowcomp2 = 6;
colcomp2 = k;

for i = 1:size(EDGES{1,1},1)
    for j = 1:size(EDGES{1,1},2)
        EDGEHEAT(i,j) = (EDGES{rowcomp1,colcomp1}(i,j)-EDGES{rowcomp2,colcomp2}(i,j))/(EDGES{rowcomp1,colcomp1}(i,j)+EDGES{rowcomp2,colcomp2}(i,j));
    end
end
EDGEHEAT(isinf(EDGEHEAT)|isnan(EDGEHEAT)) = 0;

h = heatmap(EDGEHEAT,'Colormap',map);
h.XLabel = 'Prey No.';
h.YLabel = 'Predator No.';
h.Title = 'BLUE = T5 bigger / RED = T30 bigger';
set(gcf, 'Position', [200 40 750 750]);
set(gcf, 'Color', 'w');
fname=['Output/Run_CA_N5_Tau' num2str(val.tau()) '_EiE_T5-30_gamma' num2str(val.gamma(k)) '.png'];
saveas(gcf,fname);

end



%% MONTAGE (WORK IN PROGRESS)

val.Nsupply = 5;%[3,10];
val.tau     = 1;%0.2;%[0.2, 0.4, 0.6, 0.8, 1, 1/0.8, 1/0.6, 1/0.4, 1/0.2];
val.T       = [5, 10, 15, 20, 25, 30];
val.gamma   = 0.5;%[0.2:0.15:1];
val.tempfuA = 0.05;%[0.01, 0.05, 0.1, 0.5, 1];
val.tempfuZ = 0.05;%[0.01, 0.05, 0.1, 0.5, 1];
val.time    = 3;%[1, 2, 3];

totexp = 1;
for temp = 1:length(val.T)
    for light = 1:length(val.gamma)
        for tau = 1:length(val.tau)
            img = ['Output/Run_CA_N5_Tau' num2str(val.tau(tau)) '_EiE_T' num2str(val.T(temp)) '_tfuA' num2str(val.tempfuA) '_tfuZ' num2str(val.tempfuZ) '_gamma_' num2str(val.gamma(light)) 'Network' num2str(val.time) 'D.jpg'];
            im(totexp,1) = convertCharsToStrings(img);
            %subbo = subplot(length(val.T),length(val.gamma),totexp);
        
            totexp = totexp + 1;
        end
    end
end

montage(im,'Size', [2 3],'ThumbnailSize', [])
text(.5,1.05,['Food Network'],'horiz','center','vert','top','FontSize',20,'units','normalized')
text(.53,0,['Uptake + Graze'],'horiz','center','vert','top','FontSize',15,'units','normalized')
texts = text(-0.05,.51,['Trophic level (Size Cohort)'],'horiz','center','vert','top','FontSize',15,'units','normalized');
set(texts,'Rotation',90);
set(gcf, 'Position', [10 60 1500 750]);
set(gcf, 'Color', 'w')
