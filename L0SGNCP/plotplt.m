function check=plotplt(datas,trigger,iter)
color=["-*","-p","->","-+","-x","-o"];
mess={};
check=0;
for j=1:length(trigger)
i=j;
data=datas{i};
num=1;
lossdata=data{1};
trdata=data{2};

if (abs(trdata(end)-round(trdata(end)))>0.1)
    trdata(end)=ceil(trdata(end));
end
trdata(end)=floor(trdata(end))+0.0001;

if(~isempty(lossdata))
    mess=pltplt(lossdata,trdata,num,trigger(i),mess,color(i),iter);
end
end
end

function mess=pltplt(lossdata,trdata,num,trigger,mess,color,iter)
% figure(1)
% color=["-o","-*","-+","-x","-d","-s","-p","-h"];
% subplot(1,3,iter);
ss=40;
for i=1:num
    maker_idx = 1:ss:length(lossdata);
            loss=lossdata;
%     semilogy(loss,color,'linewidth',2,'MarkerIndices',maker_idx,'MarkerSize',10);
     if(trigger==1)
    mes{i}='PAL for solving $\ell_1$-GSNTD';
    elseif(trigger==2)
    mes{i}='PAL for solving $\ell_0$-SNTD';
    elseif(trigger==3)
    mes{i}='PAL for solving $\ell_0$-GSNTD';
     elseif(trigger==4)
    mes{i}='iPAL for solving $\ell_1$-GSNTD';
    elseif(trigger==5)
    mes{i}='iPAL for solving $\ell_0$-SNTD';
    elseif(trigger==8)
    mes{i}='APGL for solving $\ell_0$-GSNCP';
    end
    hold on;
end
ls=length(mess);
for j=1:length(mes)
    mess{j+ls}=mes{j};
end
% h=legend(mess,'Interpreter','latex');
% set(gca,'FontSize',35);
% set(h,'FontSize',30)
% xlabel('Number of iterations','FontSize',35);
% ylabel('Objective funciton value','FontSize',35);
% prettyAxes().gbase2()
% 
subplot(1,3,iter);
for i=1:num
        maker_idx = 1:ss:length(lossdata);
        loss=lossdata;
    plot(trdata,loss,color,'linewidth',3,'MarkerIndices',maker_idx,'MarkerSize',12);
    hold on;
end
h=legend(mess,'Interpreter','latex');
set(gca,'FontSize',30);
set(h,'FontSize',30)
xlabel('Time (seconds)','FontSize',35);
ylabel('Objective funciton value','FontSize',35);
prettyAxes().gbase2()
    plt0=0;
end

