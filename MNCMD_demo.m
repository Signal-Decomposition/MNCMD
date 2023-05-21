% This algorithm uses  many basic functions provided by Shiqian Chen.
% Special Thanks to Shiqian Chen and his team.
clc
clear all
close all
FS=14;

L = 1000; 
fs=500;
t=(0:L-1)*(1/fs);

rng(8); 
STD=0.1;
noise=addnoise(length(t),0,STD);
f_channel1 = cos(2*pi*2*t) +cos(2*pi*24*t.^2)+ 2*(cos(2*pi*36*t))+noise; 
f_channel2 = cos(2*pi*24*t.^2) + 2*(cos(2*pi*36*t))+noise;
IFreal=[2*ones(1,L);36*ones(1,L);48*t];

x = [f_channel1;f_channel2]; 

Sig=f_channel1+f_channel2;
SampFreq=fs;
window = 256;
Nfrebin = 1024;

[Spec,f] = STFT(Sig',SampFreq,Nfrebin,window);

bw = SampFreq/80;
beta1 = 1e-4; 
num = 3; 
delta = 20;
[fidexmult, tfdv] = extridge_mult(Sig, SampFreq, num, delta, beta1,bw,Nfrebin,window);

beta = 1e-5;
iniIF = curvesmooth(f(fidexmult),beta);
s=x;
fs=fs;
alpha = 5e-6;
beta = 1e-5; 
eIF=iniIF;
var=STD^2;
tol = 1e-7;
tic
[IFmset,IA,smset]=MNCMD(s,fs,eIF,alpha,beta,var,tol);
toc
IF=IFmset(:,:,end);
IMF=smset(:,:,:,end);
[~,index]=sort(mean(IF,2));
IF=IF(index,:);
IMF=IMF(index,:,:);


figure
for i=1:3
    plot(t,IFreal(i,:),'b','linewidth',1.5)
    hold on
    plot(t,IF(i,:),'r--','linewidth',1.5)
end
legend('true','estimation')
set(gca,'xtick',[0:0.4:2])
set(gca,'ytick',[0:20:100])
xlabel('time / s')
ylabel('frequency / Hz')
text(1.5,8,'g_1','FontSize',FS)
text(1.5,42,'g_2','FontSize',FS)
text(1.5,67,'g_3','FontSize',FS)
set(gca,'FontSize',FS)

figure
subplot(421)
ylim([-5,5])
plot(x(1,:),'k')
title('x_1')
set(gca,'FontSize',FS)
subplot(423)
plot(IMF(1,:,1),'b')
ylim([-2,2])
ylabel('g_1')
set(gca,'FontSize',FS)
subplot(425)
plot(IMF(2,:,1),'b')
ylim([-5,5])
ylabel('g_2')
set(gca,'FontSize',FS)
subplot(427)
plot(IMF(3,:,1),'b')
ylim([-2,2])
ylabel('g_3')
xlabel('samples')
set(gca,'FontSize',FS)

subplot(422)
plot(x(2,:),'k')
ylim([-5,5])
title('x_2')
set(gca,'FontSize',FS)
subplot(424)
plot(IMF(1,:,2),'b')
ylim([-2,2])
set(gca,'FontSize',FS)
subplot(426)
plot(IMF(2,:,2),'b')
ylim([-5,5])
set(gca,'FontSize',FS)
subplot(428)
plot(IMF(3,:,2),'b')
ylim([-2,2])
xlabel('samples')
set(gca,'FontSize',FS)