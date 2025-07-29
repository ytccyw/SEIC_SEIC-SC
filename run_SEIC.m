clear
addpath('par-data','measure','function');
dataNAME='BDGP_fea_Per0.1';%For SEIC, it can be run by simply replacing the dataset name.
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
clu=length(unique(Y));
par=5*(V^2)/(N);
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
end
res=mean(res);
disp(['Running time:' num2str(time1)])
disp(['Kmeans time:' num2str(time2)])
disp(['【SEIC】 ' 'ACC:' num2str(res(1)*100) '%     NMI:' num2str(res(2)*100) '%     PUR:' num2str(res(6)*100) '%'])