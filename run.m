clear
addpath('par-data','measure','function');
dataNAME='CCV_Per0.1';%For SEIC, it can be run by simply replacing the dataset name.
SEIC_SC=0;%1 means running SEIC-SC, 0 means run SEIC only.
%%%%%%%%%%%%%%%%【Loading and pre-processing】%%%%%%%%%%%%%%%%
disp(['Loading:    ' dataNAME])
load(dataNAME);
disp('pre-processing...')
Y=truelabel{1};
N=size(data{1},2);
V=length(data);
for v=1:V
    data{v}(isnan(data{v})) = 0;
    data{v} = NormalizeData(data{v});
    data{v}=data{v}';
end
%%%%%%%%%%%%%%%%【Parameters】%%%%%%%%%%%%%%%%
clu=length(unique(Y));
par=5*(V^2)/(N);
if SEIC_SC==1%Parameters for SEIC-SC
    Yclu=floor(length(Y)/clu);
    anchor=clu*Yclu;%anchor=[clu*5 clu*10 clu*Yclu];
    sparse=0.05;%sparse=[0.01 0.05 0.1 0.35];
end
%%%%%%%%%%%%%%%%%【Running】%%%%%%%%%%%%%%%%%%
for repeat=1:10%Repeat run ten times
    disp(['The ' num2str(repeat) ' run....................'])
    tic
    Results = SEIC(data,clu,par);
    time1=toc;
    tic
    labels = litekmeans(Results, clu,'Replicates',10);
    time2=toc;
    res(repeat,:)=  Clustering8Measure(Y, labels);
    if SEIC_SC==1
        tic
        labels = good(Results,clu,anchor,sparse);
        time3=toc;
        RES(repeat,:)=  Clustering8Measure(Y, labels);
    end
end
res=mean(res);
disp(['Running time:' num2str(time1)])
disp(['Kmeans time:' num2str(time2)])
disp(['【SEIC】 ' 'ACC:' num2str(res(1)*100) '%     NMI:' num2str(res(2)*100) '%     PUR:' num2str(res(6)*100) '%'])
if SEIC_SC==1
    RES=mean(RES);
    disp(['Sparse anchor graph spectral clustering time:' num2str(time3)])
    disp(['【SEIC-SC】 ' 'ACC:' num2str(RES(1)*100) '%     NMI:' num2str(RES(2)*100) '%     PUR:' num2str(RES(6)*100) '%'])
end
