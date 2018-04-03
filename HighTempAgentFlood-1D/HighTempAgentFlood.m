% 未知数为：p(1)...p(14),qout出口流量
% 所有单位为国际单位
% clc
% clear
% close all
function [OR,xx] = main2(alpha, Miu)

NumX = 15;
MiuW = 1E-3;    %水粘度
MiuO = Miu*1E-3;
% MiuO = 68.37E-3;    %油粘度
K = 2.62E-12;   %渗透率
dx = 2.16E-2;   %网格宽度
dt = 0.0002;  %时间步长，定步长计算
nt = 277312;     %迭代时间次数
qont = zeros(nt,1);       %每次迭代计算后的产油量
Por = 0.4143;   %孔隙度
DenW = 1000;    %水密度
DenO = 900;     %油密度

A = zeros(NumX,NumX);       %未知数系数矩阵
B = zeros(NumX,1);      %右端项矩阵

fiw = zeros(NumX,1);    %流动系数 Flow index for water
fio = zeros(NumX,1);    %流动系数 Flow index for oil

%OR = zeros(nt,1);           %累计产量
PVtotal = NumX*dx*pi*(2.5E-2)^2/4*Por;      %总孔隙体积

%初始化
p = zeros(NumX,1);
p(end) = 0.2E6;         %出口压力保持0.2 MPa
Sw = 0.6515 + zeros(NumX,1);    %含水饱和度
SwN = zeros(NumX,1);        %存储下个时间步的饱和度值
c = ones(NumX,1);      %水在水相中的百分数
c(1) = 0.9;
cN = zeros(NumX,1);     %存储下个时间步的水百分数
qin = 8.3E-6/(pi*(2.5E-2)^2/4);      %dt时间的入口注入量,单位m/s
q = zeros(NumX,1);
q(1) = qin;
q(end) = qin;
Rs = zeros(NumX,1);     %吸附

%%%%%饱和度显示%%%%%%%%%
VA = zeros(2*(NumX+1),2);       %网格中的所有点
VA(1:(NumX+1),1) = 1:(NumX+1);
VA((NumX+2):(2*(NumX+1)),1) = 1:(NumX+1);
VA((NumX+2):(2*(NumX+1)),2) = 1;

FA = zeros(NumX,4);             %网格中的所有块
FA(:,1) = 1:NumX;
FA(:,2) = 2:(NumX+1);
FA(:,3) = (NumX+3):(2*NumX+2);
FA(:,4) = (NumX+2):(2*NumX+1);
% figure(1)
% patch('Faces',FA,'Vertices',VA,'FaceVertexCData',(1-Sw),'FaceColor','flat','Edgecolor','none');
% colormap(jet),colorbar
% xlim([1 (NumX+1)])
% set(gca,'ytick',[0 1])
% set(gca,'yticklabel',{'0','1'})
% xtk = 1:(NumX+1);
% set(gca,'xtick',xtk)

for ii = 1:nt
	for i = 1:NumX
	    fiw(i) = K*Kr(Sw(i),1)/MiuW;
	    fio(i) = K*Kr(Sw(i),2)/MiuO;
	    Rs(i) = 1E-3*(2*c(i))/(1 + 2*c(i));      % 吸附量
	end
	
	
	% the first grid
	i = 1;
	A(i,i) = -(fio(i) + fiw(i));
	A(1,i+1) = fio(i) + fiw(i);
	B(i) = fio(i)*(Pc(Sw(i+1),alpha*(1-c(i+1))) - Pc(Sw(i),alpha*(1-c(i))))- dx*(q(i)+Rs(i))/DenW;
	
	% the 2nd to NumX-2 grid
	for i = 2:(NumX-2)
	    A(i,i-1) = fio(i-1) + fiw(i-1);
	    A(i,i) = -(fio(i-1) + fiw(i-1) + fio(i) + fiw(i));
	    A(i,i+1) = fio(i) + fiw(i);
	    B(i) = fio(i)*(Pc(Sw(i+1),alpha*(1-c(i+1))) - Pc(Sw(i),alpha*(1-c(i)))) - fio(i-1)*(Pc(Sw(i),alpha*(1-c(i))) - ...
	        Pc(Sw(i-1),alpha*(1-c(i-1)))) - dx*Rs(i)/DenW;
	end
	
	% the NumX-1 grid
	i = NumX-1;
	A(i,i-1) = fio(i-1) + fiw(i-1);
	A(i,i) = -(fio(i-1) + fiw(i-1) + fio(i) + fiw(i));
	A(i,i+1) = 0;
	B(i) = fio(i)*(Pc(Sw(i+1),alpha*(1-c(i+1))) - Pc(Sw(i),alpha*(1-c(i)))) - fio(i-1)*(Pc(Sw(i),alpha*(1-c(i))) - ...
	        Pc(Sw(i-1),alpha*(1-c(i-1)))) - dx*Rs(i)/DenW - (fio(i) + fiw(i))*p(i+1);
	    
	% the last grid
	i = NumX;
	A(i,i-1) = fio(i-1) + fiw(i-1);
	A(i,i) = fio(i-1) + fiw(i-1);
	B(i) = - fio(i-1)*(Pc(Sw(i),alpha*(1-c(i))) - Pc(Sw(i-1),alpha*(1-c(i-1))))...
	    - dx*Rs(i)/DenW -dx*q(i)/(fio(i)+fiw(i))*(fio(i)/DenO + fiw(i)/DenW)...
        + dx*fio(i)*fiw(i)*(Pc(Sw(i-1),alpha*(1-c(i-1))) - Pc(Sw(i),alpha*(1-c(i))))...
	    /(fio(i)+fiw(i))*(1/DenO+1/DenW);
	
	% solve matrix equation
	Sol = A\B;
	p(1:(NumX-1)) = Sol(1:(NumX-1)) ;
	
	%显示求饱和度
	% the first grid
	i = 1;
	SwN(i) = Sw(i) + dt/Por/dx*(fiw(i)*(p(i+1)-p(i))/dx + (q(i)+Rs(i))/DenW);
	
	% the 2nd to NumX-1 grid
	for i = 2:(NumX-1)
	    SwN(i) = Sw(i) + dt/Por/dx*(fiw(i)*(p(i+1)-p(i))/dx - fiw(i-1)*(p(i)-p(i-1))/dx);
	end
	
	% the last grid
	i = NumX;
	qwo = (q(end) - fio(i)*(- Pc(i-1,alpha*(1-c(i-1))) + Pc(i,alpha*(1-c(i)))))*fiw(i)/(fiw(i)+fio(i));
	SwN(i) = Sw(i) + dt/Por/dx*(-fiw(i-1)*(p(i-1)-p(i))/dx + (qwo+Rs(i))/DenW);
	
	%显示求水质量分数
	% the first grid
	i = 1;
	cN(i) = (Sw(i)*c(i) + dt/Por/dx*(fiw(i)*c(i)*(p(i+1)-p(i))/dx + (q(i)*c(i))/DenW)...
	    )/SwN(i);
	% the 2nd to NumX-1 grid
	for i =2:(NumX-1)
	    cN(i) = (Sw(i)*c(i) + dt/Por/dx*(fiw(i)*c(i)*(p(i+1)-p(i))/dx - fiw(i-1)*c(i-1)*(p(i)-p(i-1))/dx))...
	        /SwN(i);
	end
	% the last grid
	i = NumX;
	cN(i) = (Sw(i)*c(i) + dt/Por/dx*(-fiw(i-1)*c(i-1)*(p(i-1)-p(i))/dx + (qwo*c(i)/DenW)...
	    ))/SwN(i);
	
    qont(ii) = q(end) - qwo;
    OR(ii) = sum(qont(1:ii));
	Sw = SwN;
	c = cN;
end
%patch('Faces',FA,'Vertices',VA,'FaceVertexCData',(1-Sw),'FaceColor','flat','Edgecolor','none');

% figure(2)
xx = 1:nt;
xx = 2*xx/nt;
OR = OR/10;
ORend = OR(end);
%plot(xx,OR)