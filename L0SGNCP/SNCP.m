%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
function [var,loss,timerun]=SNCP(var,ngmar,maxiteropt,stopindex,r,alphan,lamda)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(size(ngmar));

LK=zeros(1,num);
L=ones(1,num);
tk=1;
bts=[];
wk=zeros(1,num);
varK=var;
for j=1:num
wk(j)=(tk-1)/(tk);
end

for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}~=0));
end




loss(1)=computeloss(ngmar,var,lamda);


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
    if(j<num)
    [V,L(j)]=gradSNCP(vv,ngmar,j,num,r,alphan,lamda);
    else
    [V,L(j)]=gradSNCP(vv,ngmar,j,num,r,alphan,0);
    end
    var{j}=PROXn1(V);
end
loss(i+1)=computeloss(ngmar,var,lamda);


%% Judging whether to extrapolate
if(loss(i+1)>loss(i))
    var=varK;
    for j=1:num
    if(j<num)
    [V,L(j)]=gradSNCP(var,ngmar,j,num,r,alphan,lamda);
    else
    [V,L(j)]=gradSNCP(vv,ngmar,j,num,r,alphan,0);
    end
    var{j}=PROXn1(V);
    end
    loss(i+1)=computeloss(ngmar,var,lamda);
end



%% Check if termination condition is met
fprintf("SNCP\n");
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


function loss=computeloss(ngmar,var,lamda)
    loss=computeCP(var,ngmar);
    for i=1:length(size(ngmar))
    loss=loss+lamda/2*norm(var{i},1);
    end
end



function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end


function [U,L]=gradSNCP(var,ngmar,n,num,r,alpha,lamda)
  [Xtemp,temp]=krob2(var,n,num,ngmar);
   a=0;
   ck=norm(temp,'fro');
   L=ck;
   tao=1/(r*ck);
   if(temp==0) 
       while(a==0)
           a=rand(1);
       end
       tao=1/a;
   end
   A=zeros(size(var{n}));
   mar=var{n}*temp;
   for j=1:size(var{n},2)
       ar=var{n}(:,j);
       temp1=Xtemp(:,j)-mar(:,j)-lamda;
       A(:,j)=ar+tao*temp1/(temp(j,j)+alpha);
   end   
   U=A;
end