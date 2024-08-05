%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
function [var,loss,timerun,bts]=L1SGCP(var,ngmar,maxiteropt,stopindex,r,lamda)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(size(ngmar));

LK=zeros(1,num);
L=ones(1,num);
tk=1;
Lapk=LLaplace(ngmar);
bts=[];
wk=zeros(1,num);
varK=var;
for j=1:num
wk(j)=(tk-1)/(tk);
end

for i=1:num
    alpha(i)=0;
end
alpha(num)=0.1;



for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end




loss(1)=computeloss(ngmar,var,lamda,Lapk,alpha);


t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
vv=var;

for j=1:num
    wk(j)=min(wk(j),0.9999*sqrt(LK(j)/L(j)));
    %% Update parameters
    vv{j}=var{j}+wk(j)*(var{j}-varK{j});
    varK{j}=var{j};
    LK(j)=L(j);
    [V,L(j)]=gradGSNCP(vv,ngmar,j,num,r,Lapk,alpha(j));
        if(j<num)
        var{j}=PROXL1(V,lamda,1/(L(j)*r));
        else
        var{j}=PROXn1(V);
        end
end
loss(i+1)=computeloss(ngmar,var,lamda,Lapk,alpha);


%% Judging whether to extrapolate
if(loss(i+1)>loss(i))
    var=varK;
    for j=1:num
    [V,L(j)]=gradGSNCP(var,ngmar,j,num,r,Lapk,alpha(j));
        if(j<num)
        var{j}=PROXL1(V,lamda,1/(L(j)*r));
        else
        var{j}=PROXn1(V);
        end
    end
    loss(i+1)=computeloss(ngmar,var,lamda,Lapk,alpha);
end



%% Check if termination condition is met
fprintf("SGCP\n");
check1=0;
check2=0;
for j=1:num
    fprintf("nonzero Rows:%d\n",nnz(var{j}));
    check1=check1+norm(var{j}-varK{j},'fro');
    check2=check2+norm(varK{j},'fro');
end

% %% Count the number of non-zero Rows in each matrix.
% for j=1:num
%     fprintf("nonzero Rows:%d\n",sum(any(var{j},2)));
% end

bts{i}=wk;
t2=clock;
timerun(i+1)=etime(t2,t1);
% Res=abs(loss(i+1)-loss(i));
% Res=norm(tensor(S2-S1))/returnloss;
% S1=S2;
Res=check1/check2;
% Res=norm(tensor(S2-S1))/returnloss;
fprintf("Rel：%d\n",Res);
stop=stopcheck(Res,timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(4);
    break;
end
tk=(1+sqrt(1+4*tk^2))/2;
for j=1:num
wk(j)=(tk-1)/(tk);
end
end
end


function loss=computeloss(ngmar,var,lamda,Lapk,alpha)
    loss=computeCP(var,ngmar);
    for i=1:length(size(ngmar))
    loss=loss+lamda/2*norm(var{i},1)+alpha(i)/2*trace(var{i}'*Lapk{i}*var{i});
    end
end



function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end


