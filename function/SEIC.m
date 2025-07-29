function [EM] = SEIC(M,clu,alpha)
%% default paremeters setting
V=length(M);
for v=1:V
    dim{v} = size(M{v});
end
n=dim{1}(1);
tol        = 1e-2;
max_iter   = 30;
rho        = 1.7;
mu         = 1e-4;
max_mu     = 1e10;
%----------------------
Y3=zeros([n,clu,V]);
T=randn([n,clu,V]);
ebusenno=zeros([n,clu,V]);
for v=1:V
    Z{v}=zeros(n,clu);
    H{v}=zeros(clu,dim{v}(2));
end
%% FFT setting
Ftemp = zeros([n,clu,V]);
for i = 1:2
    Gamma{i} = zeros([n,clu,V]);
    G{i} = porder_diff(T,i);
    Eny = diff_element([n,clu,V],i);
    Ftemp   = Ftemp + Eny;
end
clear Eny
%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;
    for v=1:V       
        %------update Zv
        ZZZ=2*alpha*M{v}*H{v}';
        [uu,~,vv] = svd(ZZZ+mu*(-Y3(:,:,v)/mu+T(:,:,v)+ebusenno(:,:,v)),'econ');
        clear ZZZ
        Z{v}=uu*vv';
        %------update Fv
        [uu,~,vv] = svd(Z{v}'*M{v},'econ');
        H{v}=uu*vv';
    end
    %------update T
    Z_tensor = cat(3, Z{:,:});
    H_tensor = zeros([n,clu,V]);
    for i = 1:2
        H_tensor = H_tensor + porder_diff_T(G{i}+Gamma{i}/mu,i); 
    end
    T = real( ifftn(  fftn(Z_tensor-ebusenno+Y3/mu+H_tensor)    ./(1+Ftemp)     )    );
    clear H_tensor
    %------update Gi
    for i = 1:2
        TT=porder_diff(T,i);
        G{i} = prox_htnn_F(TT-Gamma{i}/mu,2*mu);
        Gamma{i} = Gamma{i}+mu*(G{i}-TT);
    end
    clear TT
    %------update ebusenno
    ebusenno = prox_l1(Z_tensor-T+Y3/mu,alpha/mu);
    %% Stop criterion
    dY   = Z_tensor-T-ebusenno;
    chg=max(abs(dY(:)));
    disp(['iter:' num2str(iter)])
    if ((chg < tol)&&(iter>=15))||(iter==max_iter)
        EM=[];
        for v = 1:V
            EM=cat(2,EM,T(:,:,v));
        end
        break
    end
    %% Update mulipliers: Y3, and mu
    Y3 = Y3+mu*dY;
    clear dY
    mu = min(rho*mu,max_mu);
end
end




