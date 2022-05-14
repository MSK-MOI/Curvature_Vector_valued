function [all_curv] = vec_curv_fun(start_pt,end_pt,alpha,gamma,dataset)
%vec_curv 
% alpha: the mass that remains at the same place for Markov random walk
% alpha>=0: fixed value alpha 
% alpha=-1: alpha=pi
% alpha=-2: alpha=1-pi
% alpha=-3: alpha=rho_i/sum_j~i(rho_j)
%
% gamma: distance between layers
% gamma>=0: fixed value gamma
% gamma=-1: gamma as correlation between two layers' each node(rescaled...need further test)
% maybe other options ...
%
% (Recomended? alpha=-1, gamma=1)
% Jiening, Phong 7/30/2021
% -------------------------------------------------------

if dataset == 1
    load('../new_datasets/metabric_format.mat','adj','CNA','RNA')
    output_name = 'Metabric';
elseif dataset == 2
    load('../new_datasets/brca_tcga.mat','adj','CNA','Methyl','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'brca';
elseif dataset == 3
    load('../new_datasets/brca_tcga.mat','adj','CNA','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'brca';
elseif dataset == 4
    load('../new_datasets/hnsc_tcga.mat','adj','CNA','RNA','Methyl')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'hnsc';
elseif dataset == 5
    load('../new_datasets/hnsc_tcga.mat','adj','CNA','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'hnsc';
elseif dataset == 6
    load('../new_datasets/laml_tcga.mat','adj','CNA','RNA','Methyl')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'laml';
elseif dataset == 7
    load('../new_datasets/laml_tcga.mat','adj','CNA','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'laml';
elseif dataset == 8
    load('../new_datasets/luad_tcga.mat','adj','CNA','RNA','Methyl')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    Methyl = Methyl + 0.001;
    Methyl = Methyl./max(Methyl(:));
    Methyl = Methyl - 0.0005;
    output_name = 'luad';
elseif dataset == 9
    load('../new_datasets/luad_tcga.mat','adj','CNA','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'luad';
elseif dataset == 10
    load('../new_datasets/lusc_tcga.mat','adj','CNA','RNA','Methyl')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'lusc';
elseif dataset == 11
    load('../new_datasets/lusc_tcga.mat','adj','CNA','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'lusc';
elseif dataset == 12
    load('../new_datasets/ov_tcga.mat','adj','CNA','RNA','Methyl')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'ov';
elseif dataset == 13
    load('../new_datasets/ov_tcga.mat','adj','CNA','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'ov';
elseif dataset == 14
    load('../new_datasets/paad_tcga.mat','adj','CNA','RNA','Methyl')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'paad';
elseif dataset == 15
    load('../new_datasets/paad_tcga.mat','adj','CNA','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'paad';
elseif dataset == 16
    load('../new_datasets/sarc_tcga.mat','adj','CNA','RNA','Methyl')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'sarc';
elseif dataset == 17
    load('../new_datasets/sarc_tcga.mat','adj','CNA','RNA')
    CNA = CNA + 2.1;
    RNA = RNA + 0.1;
    output_name = 'sarc';
end

p = gcp('nocreate');
if isempty(p)
    numcores = feature('numcores');
    myCluster = parcluster('local');
    myCluster.NumWorkers = numcores;
    saveProfile(myCluster); 
    parpool('local',numcores) 
end 

[Ng,Np]=size(CNA);% Ng: # of genes; Np: # of patients
if ~exist('Methyl','var')
    Nc=2;% Nc: # of channels w/o methylation
else
    Nc=3;% Nc: # of channels
end
[u,v]=find(triu(adj));
edge_list=sortrows([u,v],1);
Ne=size(edge_list,1);%Ne: # of edges

if gamma==-1
    if Nc==3
        gamma_i=zeros([Ng,2]);
        [~,pval] = corr(CNA',RNA');
        gamma_i(:,1)=diag(pval);
        [~,pval] = corr(CNA',Methyl');
        gamma_i(:,2)=diag(pval);
    elseif Nc==2
        gamma_i=zeros([Ng,1]);
        [~,pval] = corr(CNA',RNA');
        gamma_i(:,1)=diag(pval);
    end
    gamma_i=gamma_i*50+1;
elseif gamma>=0
    gamma_i=gamma*ones(Ng,Nc-1);
end


all_curv=zeros(Np,Ne);
all_curv_c=zeros(Np,Ne);
all_W = zeros(Np,Ne);

%parpool('local',40);
for pi=start_pt:end_pt
    filename = [output_name num2str(Nc) '_' num2str(alpha) '_' num2str(gamma) '_' num2str(pi) '.mat'];
    if ~isfile(filename)
        rho=zeros(Ng,Nc);
        inv_mea=zeros(Ng,Nc);
        if Nc==3
            rho(:,1)=CNA(:,pi);
            rho(:,2)=RNA(:,pi);
            rho(:,3)=1-Methyl(:,pi);
        elseif Nc==2
            rho(:,1)=CNA(:,pi);
            rho(:,2)=RNA(:,pi); 
        end

        if alpha==-1||alpha==-2
            inv_mea=invariant(adj,rho);
        end
        Adj_w=zeros([size(adj),Nc]);
        Adj_sw=zeros([size(adj),Nc]);
        bigAdj=zeros(Ng*Nc);
        for i=1:Nc
            Adj_temp=adj.*rho(:,i)';
            sum_temp=adj*rho(:,i);
            Adj_temp=Adj_temp./sum_temp;
            Adj_w(:,:,i)=Adj_temp;
            Adj_temp=(Adj_temp+Adj_temp')/2;
            Adj_temp(Adj_temp~=0)=1./sqrt(Adj_temp(Adj_temp~=0));
            bigAdj(1+Ng*(i-1):Ng*i,1+Ng*(i-1):Ng*i)=Adj_temp;
            Adj_sw(:,:,i)=Adj_temp;
        end

        if Nc==3
            bigAdj(1:Ng,1+Ng:Ng*2)=diag(gamma_i(:,1));
            bigAdj(1+2*Ng:3*Ng,1+Ng:Ng*2)=diag(gamma_i(:,2));
            bigAdj(1+Ng:2*Ng,1:Ng)=diag(gamma_i(:,1));
            bigAdj(1+Ng:2*Ng,1+2*Ng:3*Ng)=diag(gamma_i(:,2));
        elseif Nc==2
            bigAdj(1:Ng,(1+Ng):(Ng*2))=diag(gamma_i(:,1));
            bigAdj((1+Ng):(2*Ng),1:Ng)=diag(gamma_i(:,1));
        end
        G=graph(bigAdj);
        D=distances(G);
        parfor k=1:Ne
           u=edge_list(k,1);
           v=edge_list(k,2);
           setu0=find(adj(u,:)==1);
           setu=[setu0,u];
           setv0=find(adj(v,:)==1);
           setv=[setv0,v]; 
           rho_u=[];
           rho_v=[];
           for i=1:Nc
               fprintf('%i',alpha)
               if alpha==-1
                   alpha_i=inv_mea(u,i);
               elseif alpha==-2
                   alpha_i=1-inv_mea(u,i);
               elseif alpha==-3
                   alpha_i=rho(u,i)/sum(rho(setu,i));
               elseif alpha>=0
                   alpha_i=alpha;
               end
               rho_u=[rho_u,(1-alpha_i)*Adj_w(u,setu0,i),alpha_i];
               if alpha==-1
                   alpha_i=inv_mea(v,i);
               elseif alpha==-2
                   alpha_i=1-inv_mea(v,i);
               elseif alpha==-3
                   alpha_i=rho(v,i)/sum(rho(setv,i));
               elseif alpha>=0
                   alpha_i=alpha;
               end
               rho_v=[rho_v,(1-alpha_i)*Adj_w(v,setv0,i),alpha_i];
           end
           cc=sum(Adj_sw(u,v,:));
           if Nc==3
               DD=D([setu,setu+Ng,setu+2*Ng],[setv,setv+Ng,setv+2*Ng]);
           elseif Nc==2
               DD=D([setu,setu+Ng],[setv,setv+Ng]);
           end
           d=Kantorovich(DD,rho_u,rho_v);
           curv=1-d/cc;
           all_curv(pi,k)=curv;
           all_curv_c(pi,k)=cc-d;
           all_W(pi,k)=d;
        end
        curv_i=all_curv(pi,:);
        curv_c_i=all_curv_c(pi,:);
        W_i=all_W(pi,:);
        pi
        save(filename,'curv_i','curv_c_i','pi','W_i');
    end
end
end

