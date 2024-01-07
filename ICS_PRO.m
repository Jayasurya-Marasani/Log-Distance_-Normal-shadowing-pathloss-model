clc;clear;close all;
f=1.5e9;
d0=100;
sigma=3;
d=(1:2:31).^2;
n=[2 3 6]; 
for k=1:3
    logdist(k,:)=logdist_or_lognorm_PL_M(f,d,d0,n(k));
    lognorm(k,:)=logdist_or_lognorm_PL_M(f,d,d0,n(k),sigma);
end
subplot(2,1,1)
semilogx(d,logdist(1,:),'r-o',d,logdist(2,:),'b-o',d,logdist(3,:),'g-o')
grid on
axis([1 1000 40 110])
title(['Log-distance pathloss model for f=',num2str(f/1e6),'MHz'])
xlabel('Distance(m)'),ylabel('Path loss(dB)')
legend('n=2','n=3','n=6',"location","northwest")
subplot(2,1,2)
semilogx(d,lognorm(1,:),'r-o',d,lognorm(2,:),'b-o',d,lognorm(3,:),'g-o')
grid on
axis([1 1000 40 110])
title(['Log-normal pathloss model for f=',num2str(f/1e6),'MHz & \sigma=', num2str(sigma), 'dB'])
xlabel('Distance(m)'),ylabel('Path loss(dB)')
legend('n=2','n=3','n=6',"location","northwest")
function PL=logdist_or_lognorm_PL_M(fc,d,d0,n,sigma)
    lamda=299792458/fc;
    PL=-20*log10(lamda/(4*pi*d0))+10*n*log10(d/d0);
    if nargin>4
        PL=PL+sigma*randn(size(d));
    end
end