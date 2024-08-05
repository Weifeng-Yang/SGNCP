%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions

function [core,var,loss,timerun,bts,rate]=GSNTD(var,ngmar,maxiteropt,Rdims,stopindex,r,p,beta)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(size(ngmar));
for i=1:num-1
    Rdims(i)=10;
end
for i=1:num
    var{i}=rand(size(ngmar,i),Rdims(i));
end
core=zeros(Rdims);
core=tenrand(size(core));
core=tensor(core);

LK=zeros(1,num);
LCK=0;
L=ones(1,num);
LC=1;
tk=1;
beta=repmat(beta,num,1);
beta(num)=0;

Lapk=LLaplace(ngmar);
varK=var;
coreK=core;
wk=(tk-1)/(tk);

returnloss=norm(tensor(ngmar));
for i=1:num-1
    alpha(i)=0;
end
alpha(num)=1000;
loss(1)=computeloss(ngmar,core,var,p,alpha,Lapk,beta);


rate=[0];
t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);

LtempK=LK;
LCtempK=LCK;
vartempK=varK;
coretempK=coreK;    
varK=var;
coreK=core;
LK=L;
LCK=LC;    
[core,var,LC,L,btc,bt]=GSNTDupdate(core,var,coretempK,vartempK,num,ngmar,r,wk,L,LC,LtempK,LCtempK,p,Lapk,alpha,beta);
loss(i+1)=computeloss(ngmar,core,var,p,alpha,Lapk,beta);


% Judging whether to extrapolate
if(loss(i+1)>=loss(i))
    var=varK;
    core=coreK;
    L=LK;
    LC=LCK;
    [core,var,LC,L]=GSNTDupdate(core,var,coretempK,vartempK,num,ngmar,r,0,L,LC,LtempK,LCtempK,p,Lapk,alpha,beta);
    loss(i+1)=computeloss(ngmar,core,var,p,alpha,Lapk,beta);
end


%% Check if termination condition is met

fprintf("GSNTD\n");

fprintf("nonzero:%d\n",nnz(core));  

check1=norm(tensor(core));
check2=norm(tensor(coreK));
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));  
    check1=check1+norm(var{j}-varK{j},'fro');
    check2=check2+norm(varK{j},'fro');
end

bts{i}=[btc,bt];
t2=clock;
timerun(i+1)=etime(t2,t1);
% Res=abs(loss(i+1)-loss(i));
Res=check1/check2;
fprintf("cri：%d\n",Res);
stop=stopcheck(Res,timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end

tk=(1+sqrt(1+4*tk^2))/2;

wk=(tk-1)/(tk);

if(p<2)
    wk=0;
end

end
end


function loss=computeloss(ngmar,core,var,p,alpha,Lapk,beta)
    loss=compute(core,var,ngmar);
    for i=1:length(size(ngmar))-1
        temp=abs(var{i}).^p;
        loss=loss+beta(i)*sum(temp(:));
    end
    for i=1:length(size(ngmar))
        loss=loss+alpha(i)/2*trace(var{i}'*Lapk{i}*var{i});
    end
end
