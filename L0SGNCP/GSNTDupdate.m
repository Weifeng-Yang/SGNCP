function [core,var,LC,L,btc,bt]=GSNTDupdate(core,var,coreK,varK,num,ngmar,r,wk,L,LC,LK,LCK,p,Lapk,alpha,beta)

    btc=min(wk,0.99*sqrt(LCK/LC));
    core=core+btc*(core-coreK);
    [V,LC]=gradGSNTDcore(core,var,ngmar,r,num,Lapk,alpha(num));
    core=PROXn1(V);
    
    for j=1:num
    bt(j)=min(wk,0.99*sqrt(LK(j)/L(j)));
    var{j}=var{j}+bt(j)*(var{j}-varK{j});
    if(j<num)
    [V,L(j)]=gradGSNTD(core,var,ngmar,r,j,num,Lapk,alpha(j),p,beta(j));
    else
    [V,L(j)]=gradGSNTD(core,var,ngmar,r,j,num,Lapk,alpha(j),p,0);    
    end
    var{j}=PROXn1(V);
    end

end


function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end

function [V,L]=gradGSNTD(core,var,ngmar,r,n,num,Lapk,alpha,p,beta)
core=tensor(core);
index=1:num;
index(n)=[];
coreg=ttm(core,var,index);
tempB=double(tenmat(coreg,n));
temp=tempB*tempB';
Xn=double(tenmat(ngmar,n));
if(n==num)
corenum=double(tenmat(core,num));
coretemp=corenum*corenum';
U=var{n}*temp-Xn*tempB'+alpha*Lapk{n}*var{n}*coretemp;
L=norm(temp,'fro')+alpha*norm(Lapk{n},'fro')+alpha*norm(coretemp,'fro'); %%
else
U=var{n}*temp-Xn*tempB'+alpha*Lapk{n}*var{n}+beta*p*var{n}.^(p-1);
L=norm(temp,'fro')+alpha*norm(Lapk{n},'fro')+(beta*p)^(1/(p-1)); %%
end
V=var{n}-1/(r*L)*U;
end


function [V,LC,U]=gradGSNTDcore(core,var,ngmar,r,num,Lapk,alpha)
ngmar=tensor(ngmar);
core=tensor(core);
Vc=core;
Vx=ngmar;
temp=1;

for i=1:num
    vart=var{i}'*var{i};
    Vc=ttm(Vc,vart,i);
    Vx=ttm(Vx,var{i}',i);
    temp=temp*norm(vart,'fro');
end
U=Vc-Vx+alpha*ttm(core,var{num}'*Lapk{num}*var{num},num);
LC=temp+alpha*norm(var{num}'*Lapk{num}*var{num},'fro');
V=core-1/(r*temp)*U;

end
