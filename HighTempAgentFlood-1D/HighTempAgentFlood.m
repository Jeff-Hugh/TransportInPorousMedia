% δ֪��Ϊ��p(1)...p(14),qout��������
% ���е�λΪ���ʵ�λ
% clc
% clear
% close all
function [OR,xx] = main2(alpha, Miu)

NumX = 15;
MiuW = 1E-3;    %ˮճ��
MiuO = Miu*1E-3;
% MiuO = 68.37E-3;    %��ճ��
K = 2.62E-12;   %��͸��
dx = 2.16E-2;   %������
dt = 0.0002;  %ʱ�䲽��������������
nt = 277312;     %����ʱ�����
qont = zeros(nt,1);       %ÿ�ε��������Ĳ�����
Por = 0.4143;   %��϶��
DenW = 1000;    %ˮ�ܶ�
DenO = 900;     %���ܶ�

A = zeros(NumX,NumX);       %δ֪��ϵ������
B = zeros(NumX,1);      %�Ҷ������

fiw = zeros(NumX,1);    %����ϵ�� Flow index for water
fio = zeros(NumX,1);    %����ϵ�� Flow index for oil

%OR = zeros(nt,1);           %�ۼƲ���
PVtotal = NumX*dx*pi*(2.5E-2)^2/4*Por;      %�ܿ�϶���

%��ʼ��
p = zeros(NumX,1);
p(end) = 0.2E6;         %����ѹ������0.2 MPa
Sw = 0.6515 + zeros(NumX,1);    %��ˮ���Ͷ�
SwN = zeros(NumX,1);        %�洢�¸�ʱ�䲽�ı��Ͷ�ֵ
c = ones(NumX,1);      %ˮ��ˮ���еİٷ���
c(1) = 0.9;
cN = zeros(NumX,1);     %�洢�¸�ʱ�䲽��ˮ�ٷ���
qin = 8.3E-6/(pi*(2.5E-2)^2/4);      %dtʱ������ע����,��λm/s
q = zeros(NumX,1);
q(1) = qin;
q(end) = qin;
Rs = zeros(NumX,1);     %����

%%%%%���Ͷ���ʾ%%%%%%%%%
VA = zeros(2*(NumX+1),2);       %�����е����е�
VA(1:(NumX+1),1) = 1:(NumX+1);
VA((NumX+2):(2*(NumX+1)),1) = 1:(NumX+1);
VA((NumX+2):(2*(NumX+1)),2) = 1;

FA = zeros(NumX,4);             %�����е����п�
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
	    Rs(i) = 1E-3*(2*c(i))/(1 + 2*c(i));      % ������
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
	
	%��ʾ�󱥺Ͷ�
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
	
	%��ʾ��ˮ��������
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