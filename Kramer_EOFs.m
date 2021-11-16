%Sasha Kramer
%sasha.kramer@lifesci.ucsb.edu
%UCSB IGPMS

%%%Code to perform Empirical Orthogonal Function (EOF) analysis on HPLC pigments and plot loadings + modes
%Lines 7-72 also appear in Kramer_cluster.mat
%Map the directory where you will load your data:
cd /Users/skramer/Documents/UCSB/Research/Data/HPLC_Aph_Rrs/

%Load your samples (formatted here as a .mat file):
load Global_HPLC_all.mat  

%The order of columns in this matrix (called "Global_HPLC") is:
%Tchla,Tchlb,Tchlc,ABcaro,ButFuco,HexFuco,Allo,Diadino,Diato,Fuco,Perid,Zea,MVChla,DVchla,Chllide,MVChlb,DVchlb,Chlc12,Chlc3,Lut,Neo,Viola,Phytin,Phide,Pras

%Initial quality control - set all values below detection (based on NASA GSFC limits) to 0
%0.001 = Chlc3, Chlc12, Chllide, Viola, Diadino, Diato, Allo, Zea, Lut, ABcaro
indPH = find(Global_RHPLC(:,24) ~= 0 & Global_RHPLC(:,24) <= 0.003); %Phide
Global_RHPLC(indPH,24) = 0; clear indPH

indPE = find(Global_RHPLC(:,11) ~= 0 & Global_RHPLC(:,11) <= 0.003); %Perid
Global_RHPLC(indPE,11) = 0; clear indPE

indBF = find(Global_RHPLC(:,5) ~= 0 & Global_RHPLC(:,5) <= 0.002); %ButFuco
Global_RHPLC(indBF,5) = 0; clear indBF

indFU = find(Global_RHPLC(:,10) ~= 0 & Global_RHPLC(:,10) <= 0.002); %Fuco
Global_RHPLC(indFU,10) = 0; clear indFU

indNE = find(Global_RHPLC(:,21) ~= 0 & Global_RHPLC(:,21) <= 0.002); %Neo
Global_RHPLC(indNE,21) = 0; clear indNE

indPR = find(Global_RHPLC(:,25) ~= 0 & Global_RHPLC(:,25) <= 0.002); %Pras
Global_RHPLC(indPR,25) = 0; clear indPR

indHF = find(Global_RHPLC(:,6) ~= 0 & Global_RHPLC(:,6) <= 0.002); %HexFuco
Global_RHPLC(indHF,6) = 0; clear indHF

indCB = find(Global_RHPLC(:,16) ~= 0 & Global_RHPLC(:,16) <= 0.003); %MVchlb
Global_RHPLC(indCB,16) = 0; clear indCB

indDC = find(Global_RHPLC(:,14) ~= 0 & Global_RHPLC(:,14) <= 0.002); %DVchla
Global_RHPLC(indDC,14) = 0; clear indDC

indPY = find(Global_RHPLC(:,23) ~= 0 & Global_RHPLC(:,23) <= 0.003); %Phytin
Global_RHPLC(indPY,23) = 0; clear indPY

%Check percent of pigments below detection (if >80% of values are below detection, I
%remove the pigment from further analysis):
for i = 1:25
    belowD = find(Global_RHPLC(:,i) <= 0.001);
    j(i) = length(belowD);
    percent(i) = 100*(j(i)./145); %need to replace the denominator with your # of samples (here, 145)
end
clear belowD i j

%First remove degradation products: Chllide,Phytin,Phide
deg = [15,23,24];
Rpigcluster1 = Global_RHPLC;
Rpigcluster1(:,deg) = [];

%Remove redundant pigs (Tchlb, Tchlc, ABcaro, MVchla):
%In my dataset, Lut, Pras, and DVchlb were also below detection >80% of the
%time so I remove them here
deg2 = [2,3,4,8,9,13,16,19,22]; %remember that the order of pigments changed when you removed degradation pigments above
Rpigcluster2 = Rpigcluster1;
Rpigcluster2(:,deg2) = [];
label2 = {'TChla','ButFuco','HexFuco','Allo','Fuco','Perid','Zea','DVchla','MVchlb','Chlc12','Chlc3','Neo','Viola'};

%Normalize to chlorophyll-a:
normchl = Rpigcluster2(:,2:end)./Global_RHPLC(:,1);
normlabel = label2(2:end);

%Create a color vector for future plotting (I choose the # of colors based
%on the results of the hierarchical cluster analysis):
colors = jet(5);

%Define a vector for calculating EOFs based on pigments normalized to Tchla:
%I sort the pigment ratios in order of the cluster results:
normchlS = [normchl(:,8),normchl(:,11),normchl(:,12),normchl(:,3),normchl(:,1),normchl(:,2),normchl(:,4),normchl(:,9),normchl(:,10),normchl(:,5),normchl(:,6),normchl(:,7)];
pigEOF = normchlS;

%Mean center pigments before calculating EOFs
pigmeans = nanmean(pigEOF);
pigstd = nanstd(pigEOF);
for i = 1:size(pigEOF,1);
    pigs_center(i,:) = (pigEOF(i,:) - pigmeans)./pigstd;
end
clear i

%Calculate EOFs: Amplitude Functions (AFs_pig) and loadings (EOFs_pig)
[EOFs_pig,AFs_pig,eigvalues_pig] = pca(pigs_center,'Centered',false,'Rows','complete'); 

%Look at the percent variance explained by all modes
var_explained = (eigvalues_pig(1:end)/sum(eigvalues_pig))';
var_cutoff = 0.031; %change this cutoff based on how many modes you want to keep

%For plotting later on:
expvar = var_explained(var_explained>var_cutoff)';
expvar_percent = expvar.*100;
perc_exp = [expvar_percent cumsum(expvar_percent)]
expvar_str = num2str(expvar_percent,'%.1f%%s'); 

%Remove modes that don't explain much variance
%Keeping top 6 modes so everything less than 3% gets cut...
EOFs_pigC = EOFs_pig(:,var_explained>var_cutoff);
AFs_pigC = AFs_pig(:,var_explained>var_cutoff);
latlon = [Global_Rlat,Global_Rlon];
AFs_pig_latlon = [latlon,AFs_pigC];

%Amplitude functions:
AF1p = AFs_pigC(:,1);
AF2p = AFs_pigC(:,2);
AF3p = AFs_pigC(:,3);
AF4p = AFs_pigC(:,4);
AF5p = AFs_pigC(:,5);
AF6p = AFs_pigC(:,6);

%Correlating AFs of each mode with pigments and pig:chl ratios
%Add latitude and longitude to check correlation with modes
corrmat1 = [AF1p,AF2p,AF3p,AF4p,AF5p,AF6p,pigEOF,latlon]; %latlon = your latitude and longitude vector
[R_p1,P_p1] = corrcoef(corrmat1,'rows','pairwise');

corrmat = [AF1p,AF2p,AF3p,AF4p,AF5p,AF6p,pigEOF];
[R_pig,P_pig] = corrcoef(corrmat,'rows','pairwise');
%check here if EOF's are correlated with each other - confirm that they are not!

%Trim R_pig to only include correlations of pigments with AFs
R_pig_trim = R_pig(1:6,7:end); %6 is number of modes you want to keep/plot, 7 is that +1 --> change if you want to look at more modes
P_pig_trim = P_pig(1:6,7:end); %trims to only show correlations of AF's with pigments

%Prepare to make some plots: put pigs & labels in same order as clusters
plotlabel = {'MVchlb','Neo','Viola','Allo','ButFuco','HexFuco','Fuco','Chlc12','Chlc3','Perid','Zea','DVchla'};

%Rearrange EOF loadings and corrcoef's to match cluster results
EOFs_rearr = EOFs_pigC;

R_pigR = R_pig_trim';
R_rearr = R_pigR;

%Here's a way to plot the loadings (repeat for modes 1-X):
ys = [-1 1];
xs = [0 (size(pigEOF,2)+1)]; 
figure(5),clf
hold on
b1 = bar([1:4],EOFs_rearr(1:4,1)); %change the column of the matrix here to plot other modes
b1b = bar([5:6],EOFs_rearr(5:6,1));
b1c = bar([7:9],EOFs_rearr(7:9,1));
b1d = bar([10],EOFs_rearr(10,1));
b1e = bar([11:12],EOFs_rearr(11:12,1));
set(b1,'FaceColor',colors(3,:))
set(b1b,'FaceColor',colors(2,:)) 
set(b1c,'FaceColor',colors(4,:)) 
set(b1d,'FaceColor',colors(5,:)) 
set(b1e,'FaceColor',colors(1,:))
set(gca,'XTick',1:length(EOFs_rearr),'XTickLabel',[]);
ylim(ys);
xlim(xs);
text(0.3,0.8,['Mode 1: ',expvar_str(1,1:5)],'FontSize',26,'FontWeight','bold') %add percent var to top right within plot
x = [1:1:(size(pigEOF,2)+1)]; %number of bars
y = ones(1,size(pigEOF,2))*-0.95; 
xticklabels(plotlabel)
ax = gca;
ax.YGrid = 'on'; %adding horizontal grid lines
ax.XTickLabelRotation=90; 
ax.YTick = [-1:0.5:1];
set(gca,'fontsize',24)
box on
for i1=1:length(R_rearr)
    text(x(i1)-0.25,-0.85,num2str(R_rearr(i1,1)*100,'%0.0f'),'FontSize',15,'fontweight','bold')
end
clear b1*

%Here's an example of plotting the amplitude functions in space:
figure('Color','white'),clf
ax = worldmap([-80 80],[-180 180]);
load coastlines
plotm(coastlat,coastlon,'k','LineWidth',1.5)
gridm('off')
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8])
hold on
scatterm(Global_Rlat,Global_Rlon,70,AF1p,'filled'),colorbar
caxis([-2 2])
h = colorbar;
cptcmap('GMT_polar')
ylabel(h,'AF1 magnitude')
setm(ax,'mlabelparallel',-90,'fontsize',20,'fontweight','bold','mlabellocation',[-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180],'plabellocation',[-80 -60 -40 -20 0 20 40 60 80])
set(h,'fontsize',18)
clear ax h coastlat coastlon c i b1 b2 b3 l