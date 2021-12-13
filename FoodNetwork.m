clear
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
%% SELECT GRAPH, NET, DUMMY
JUMPGibbs = 1;
JUMPDalel = 1;
DUMMY = 0;
DALDUMMY = 1;

%% SELECT THRESHOLDS

%If the EDGEquant = 0 the code breks
FNquant = 0; 
EDGEquant = exp(-20);

val.Nsupply = [5];%[3,10];
val.tau     = [0.2];%[0.2, 0.4, 0.6, 0.8, 1, 1/0.8, 1/0.6, 1/0.4, 1/0.2];
val.T       = [20];%[5, 10, 15, 20, 25, 30];
val.gamma   = [0.2];%[0.2:0.15:1];
val.tempfuA = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.tempfuZ = [0.05];%[0.01, 0.05, 0.1, 0.5, 1];
val.time    = [1000];%[1, 2, 3];

totexp = 1;
for temp = 1:length(val.T)
    for light = 1:length(val.gamma)
        for tau = 1:length(val.tau)
%% LOAD DATA
path =  ['Output/Run_CA_N5_Tau' num2str(val.tau(tau)) '_EiE_T' num2str(val.T(temp)) '_tfuA' num2str(val.tempfuA) '_tfuZ' num2str(val.tempfuZ) '_gamma' num2str(val.gamma(light)) '/Network' num2str(val.time) '.mat'];
load(path)
prdpry = eco.p.prdpry;
P = eco.P;
FN = P;
Nutri = yout(53,1);
yout = reshape(yout(53,2:end),[sqrt(length(P)),sqrt(length(P))]);

ind = FN < 0;
FN(ind) = 10^-40;
clear ind

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FNs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
%% FN POSITION
FNindex=0;
for i = 1:sqrt(length(P))
   po = repelem(i,1,sqrt(length(P)));
   FNindex = cat(2,FNindex,po);
end
FNindex(1) = [];

poso=repmat(1:sqrt(length(P)),1,sqrt(length(P)));
FNindex = cat(1,FNindex,poso);

FNindex = FNindex';

clear po poso

%% FN THRESHOLD
FNthresh = find(FN > FNquant);

%% FN COLOUR
auto = [0     1     0];           %GREEN
mixo(1,:) = [0.9    1     0];     %YELLOW
mixo(2,:) = [0.9    1     0];     %YELLOW
mixo(3,:) = [0.9    1     0];     %YELLOW
mixo(4,:) = [0.9    1     0];     %YELLOW
mixo(5,:) = [1    0.8     0];     %ORANGE
mixo(6,:) = [1    0.8     0];     %ORANGE
mixo(7,:) = [1    0.8     0];     %ORANGE
mixo(8,:) = [1    0.8     0];     %ORANGE
mixo(9,:) = [1    0.8     0];     %ORANGE
hete = [1     0     0];           %RED

FNcolour = cat(1,auto,mixo,hete);
FNcolour = repmat(FNcolour,length(yout),1);

clear auto mixo hete

%% FN LABEL
%FNlabel = 1:122;
FNlabel = {'CY' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'HBA'  ...
           'CY' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'HBA'  ...
           'CY' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'MBA' 'HBA'  ...
           'CY-DI' 'MDC-MBA' 'MDC-MBA' 'MDC-MBA' 'MDC-MBA' 'MDC-MBA' 'MDC-MBA' 'MDC-MBA' 'MDC-MBA' 'MDC-MBA' 'HDC-HBA'  ...
           'CY-DI' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'HDC'  ...
           'DI' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'HDC'  ...
           'DI' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'HDC'  ...
           'DI' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'HDC-CO'  ...
           'DI' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'HDC-CO'  ...
           'DI' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'HDC-CO'  ...
           'DI' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'MDC' 'HDC-CO' ' '};       

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDGEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
%% EDGE WIDTH GOOD
P;
for u = 1:size(prdpry,1)
    for m = 1:length(P)
        pry_P(u,m) = prdpry(u,m)*P(m);
    end
end
GGain  = eco.g_refg .* pry_P;
EDGE = P.*GGain.*eco.p.lambda;
EDGE = reshape(EDGE,[length(EDGE)*length(EDGE),1]);


clear u refgiu p

%% EDGE THRESHOLD
EDGEthresh = (find(EDGE > EDGEquant))';

%% EDGE INDEX MAKER
edgesecond = 0;
for i = 1:length(P)
   EDGEindex = repelem(i,1,length(P));
   edgesecond = cat(2,edgesecond,EDGEindex);
end
edgesecond(1) = [];

edgefirst = repmat(1:length(P),1,length(P));
EDGEindex = cat(1,edgefirst, edgesecond);
EDGEindex = EDGEindex';

clear edgefirst edgesecond i m 

%% EDGE COLOUR
uptake = [0.3922    0.8314    0.0745];
graze = [0.6118         0         0];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THRESHOLD RECONCILIATOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
%% FNs and EDGEs compatible indexes
a = ismember(EDGEindex,FNthresh);
b = find(1 == a(:,1));
c = find(1 == a(:,2));
o = 1;
for i = 1:length(b)
   for j = 1:length(b)
      if b(i) == c(j)
         e(o) = b(i);
         o = o+1;
      end 
   end
end
EDGEthresh = EDGEthresh(ismember(EDGEthresh,e));
EDGEthresh = EDGEthresh';

EDGEindReco = cat(2,EDGEindex(EDGEthresh,1),EDGEindex(EDGEthresh,2));
o = 1;
for i = 1:length(EDGEindReco)
   for j = 1:length(EDGEindex)
      if EDGEindReco(i,:) == EDGEindex(j,:)
         EDGEReco(o) = EDGE(j);
         o = o+1;
      end 
   end
end

clear i j o a b c d e r f 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NUTRIENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
%% NUTRIENT SIZE/POSITION/THRESHOLD
FN=cat(1,FN,Nutri); 

FNindex=cat(1,FNindex,[1 .5]);

nutupk = eco.Nuptk.*P;
Nuptkthresh = find(nutupk > EDGEquant);
Nuptkthresh = cat(2,Nuptkthresh,length(FN)+zeros(length(Nuptkthresh),1));
EDGEindReco=cat(1,EDGEindReco,Nuptkthresh);
EDGEReco=cat(2,EDGEReco,eco.Nuptk(Nuptkthresh(:,1))');

clear Nuptk Nuptkthresh nutupk



%% SORTING THE THRESHOLDS TO MAKE GRAPH HAPPY :(
EDGEindReco_sort = cat(2,EDGEindReco,EDGEReco');
EDGEindReco_sort = full(EDGEindReco_sort);

sorto = graph(EDGEindReco_sort(:,1),EDGEindReco_sort(:,2),1:length(EDGEindReco_sort),length(FN));
sorting = table2array(sorto.Edges(:,2));
EDGEindReco_sort = EDGEindReco_sort(sorting,:);

clear EDGEReco sorto sorting 

%% NUTRIENT FN/EDGE COLOUR

FNcolour=cat(1,FNcolour,[0.6196 0.4549 0]);

EDGEcolour = zeros(length(EDGEindReco_sort),3);
for w = 1:length(EDGEindReco_sort)
    if length(FN) == EDGEindReco_sort(w,2)
        EDGEcolour(w,:) = uptake;
    else
        EDGEcolour(w,:) = graze;
    end
end

clear w uptake graze EDGEindReco

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
%% SQRT CORRECTION FOR SCALING
FN(:,2) = nthroot(FN(:,1),3)*50;
ind = FN(:,2) < 0;
FN(ind,2) = 10^-40;
EDGEindReco_sort(:,4) = nthroot(EDGEindReco_sort(:,3),3)*30;

clear ind

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DUMMIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
%% DALELIO DUMMY NETWORK
if DALDUMMY == 1
    
    %NODES
    xls = xlsread('nodecompa.xlsx');
    FN(:,2) = zeros(length(FN),1)+0.000001;
    FN(xls(:,3),2) = nthroot(xls(:,4),3)*3;
    
    %EDGES
    %xls = xlsread('edgecompa.xlsx');
    xls = xlsread('edgecompa-33sort.xlsx');
    xls = xls(find(xls(:,3)),:);
    EDGEindReco_sort = xls;
    EDGEindReco_sort(:,4) = sqrt(xls(:,3));
    
    uptake = [0.3922    0.8314    0.0745];
    mixgraze = [1    0.8     0];
    zoograze = [0.6118         0         0];
    
    sorto = graph(EDGEindReco_sort(:,1),EDGEindReco_sort(:,2),1:length(EDGEindReco_sort),length(FN));
    sorting = table2array(sorto.Edges(:,2));
    EDGEindReco_sort = EDGEindReco_sort(sorting,:);
    
    

    het = [11 22 33 44 55 66 77 88 99 110 121];
    %colors edges
    clear EDGEcolour
    for w = 1:length(EDGEindReco_sort)
        if EDGEindReco_sort(w,2) == 122
            EDGEcolour(w,:) = uptake;
        end
         y = EDGEindReco_sort(w,1);
        if y == 11 || y == 22 || y == 33 || y == 44 || y == 55 || y == 66 || y == 77 || y == 88 || y == 99 || y == 110 || y == 121
            EDGEcolour(w,:) = zoograze;
        else
            EDGEcolour(w,:) = mixgraze;
        end
    end

end
clear i sort sorting uptake graze y

%% DUMMY NETWORK
if DUMMY == 1
    EDGEindReco_sort(2:end,:) = [];
    EDGEcolour(2:end,:) = [];
    FN(:,2) = zeros(1,length(FN))+10;
    EDGEindReco_sort(2:end,:) = [];
    EDGEcolour(2:end,:) = [];
    clear i
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NETWORK GIBBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&
%%
if JUMPGibbs == 1
figure(1)

% GRAPH PLOT MAKER
NWG = graph(EDGEindReco_sort(:,1),EDGEindReco_sort(:,2),0,length(FN));
NWGobject = plot(NWG,'LineWidth',EDGEindReco_sort(:,4),'MarkerSize', FN(:,2), ...
           'NodeColor', FNcolour,'EdgeColor', EDGEcolour,'XData', FNindex(:,1), 'YData', FNindex(:,2));
           %'NodeLabel', FNlabel);
       
% AXIS AND TEXT MAKER
set(gca, 'XDir','reverse')
camroll(-90)
axis off
text(.5,1,['Food Network'],'horiz','center','vert','top','FontSize',20,'units','normalized')

%text(1,1,[' CY = Cyanobacteria \newline CO = Copepods \newline DI = Diatoms \newline HBA = Hete. Bacteria \newline HDC = Hete. Dino/Cilia \newline MBA = Mixo. Bacteria \newline MDC = Mixo. Dino/Cilia'],'horiz','center','vert','top','FontSize',10,'units','normalized')

text(.121,.15,['N'],'horiz','center','vert','top','FontSize',20,'units','normalized')

text(.53,0.05,['Feeding Strategy'],'horiz','center','vert','top','FontSize',15,'units','normalized')
text(.16,0.05,['Autotrophs'],'horiz','center','vert','bottom','FontSize',10,'units','normalized')
text(.53,0.05,['Mixotrophs'],'horiz','center','vert','bottom','FontSize',10,'units','normalized')
text(.88,0.05,['Heterotrophs'],'horiz','center','vert','bottom','FontSize',10,'units','normalized')

texts = text(.045,.51,['Size'],'horiz','center','vert','top','FontSize',15,'units','normalized');
textsdnan = text(.075,.15,['Nanoplankton'],'horiz','center','vert','top','FontSize',10,'units','normalized');
textsdmic = text(.075,.51,['Microplankton'],'horiz','center','vert','top','FontSize',10,'units','normalized');
textsdmes = text(.075,.87,['Mesoplankton'],'horiz','center','vert','top','FontSize',10,'units','normalized');
set(textsdnan,'Rotation',90);
set(textsdmic,'Rotation',90);
set(textsdmes,'Rotation',90);
set(texts,'Rotation',90);

set(gcf, 'Position', [200 60 1300 750]);
set(gcf, 'Color', 'w');
fname=['Output/Run_CA_N5_Tau' num2str(val.tau(tau)) '_EiE_T' num2str(val.T(temp)) '_tfuA' num2str(val.tempfuA) '_tfuZ' num2str(val.tempfuZ) '_gamma_' num2str(val.gamma(light)) 'Network' num2str(val.time) 'G.png'];
saveas(gcf,fname);

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NETWORK D'ALELIO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% &&
%% NUTRIENT POSITION
if JUMPDalel == 1
FNindex(end,:) = [0 .1];

%% LIVING BIOMASS POSITIONS
% transforms the quantity of biomass GIVED between levels in coordinates for the nodes
FNindex(:,3) = zeros(length(FNindex),1);
for i = 1:length(FNindex)
   for j = 1:length(EDGEindReco_sort)
        if i == EDGEindReco_sort(j,1)
             u = EDGEindReco_sort(j,3);
             FNindex(i,3) = FNindex(i,3) + u;
        end
   end
end

FNindex(:,4) = zeros(length(FNindex),1);
for i = 1:length(FNindex)
   for j = 1:length(EDGEindReco_sort)
        if i == EDGEindReco_sort(j,2)
             k = EDGEindReco_sort(j,3);
             FNindex(i,4) = FNindex(i,4) + k;
        end
   end
end

FNindex(:,5) = FNindex(:,3) + FNindex(:,4);

FNindex(isinf(FNindex)|isnan(FNindex)) = 0;

%% GRAPHIC FEATURES
% create small distances between different trophic strategies of same sizes
for v = 1:sqrt(length(P))
    smallgap(v,:) =  0 +0.09 * v;
end
smallgap = cat(1,smallgap,smallgap,smallgap,smallgap,smallgap,smallgap,smallgap,smallgap,smallgap,smallgap,smallgap,0);
FNindex(:,1) = FNindex(:,1) + smallgap;

clear i j u k v

%% GRAPH MAKER

figure(2)
% GRAPH PLOT MAKER
NWD = graph(EDGEindReco_sort(:,1),EDGEindReco_sort(:,2),0,length(FN));
NWDobject = plot(NWD,'LineWidth',EDGEindReco_sort(:,4),'MarkerSize', FN(:,2), ...
           'NodeColor', FNcolour,'EdgeColor', EDGEcolour,'XData', FNindex(:,3), 'YData', FNindex(:,1));
           %'NodeLabel', FNlabel);

% AXIS AND TEXT MAKER
axis on
text(.5,1.05,['Food Network'],'horiz','center','vert','top','FontSize',20,'units','normalized')

%text(.92,.92,[' CY = Cyanobacteria \newline CO = Copepods \newline DI = Diatoms \newline HBA = Hete. Bacteria \newline HDC = Hete. Dino/Cilia \newline MBA = Mixo. Bacteria \newline MDC = Mixo. Dino/Cilia'],'horiz','center','vert','top','FontSize',10,'units','normalized')

text(.53,-0.05,['Nutrient Uptake (Green) + Graze Gain (Red)'],'horiz','center','vert','top','FontSize',15,'units','normalized')

texts = text(-0.05,.51,['Trophic level (Size Cohort)'],'horiz','center','vert','top','FontSize',15,'units','normalized');
set(texts,'Rotation',90);

set(gcf, 'Position', [200 40 1300 750]);
set(gcf, 'Color', 'w');
fname=['Output/Run_CA_N5_Tau' num2str(val.tau(tau)) '_EiE_T' num2str(val.T(temp)) '_tfuA' num2str(val.tempfuA) '_tfuZ' num2str(val.tempfuZ) '_gamma_' num2str(val.gamma(light)) 'Network' num2str(val.time) 'D.png'];
saveas(gcf,fname);

end
%%
clearvars -except totexp temp light tau val FNquant EDGEquant JUMPGibbs JUMPDalel DUMMY DALDUMMY

totexp = totexp + 1;
        end
    end
end
%%

