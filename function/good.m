function groups = good(em,n,maodian,baifenbi)
warning off;
N = size(em,1);
[~, ~, ~, ~, dis] = litekmeans(em,maodian);
Mean = mean(dis,2)+eps;
zhi = exp(-(dis.^2)./Mean);
% zhi=1./(1+dis);
rowmin = min(zhi,[],2);rowmax = max(zhi,[],2);zhi=rescale(zhi,"InputMin",rowmin,"InputMax",rowmax);
sumS=sum(zhi,2);
S=zeros(size(zhi));
[aaa,bbb]=sort(zhi,2,'descend');
if baifenbi<0.8
    sumS=sumS*baifenbi;
    for a=1:N
        t=1;temp=0;
        while temp<sumS(a)
            jiajia=aaa(a,t);
            S(a,bbb(a,t))=jiajia;
            t=t+1;temp=temp+jiajia;
        end
    end
else
    sumS=sumS-sumS*baifenbi;
    S=zhi;
    for a=1:N
        t=maodian;temp=0;
        while temp<sumS(a)
            jiajia=aaa(a,t);
            S(a,bbb(a,t))=0;
            t=t-1;temp=temp+jiajia;
        end
    end
end

if (numel(S)-nnz(S))/numel(S)>0.9
    S=sparse(S);
end
groups = Tcut_for_bipartite_graph(S,n);
