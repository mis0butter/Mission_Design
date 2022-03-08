clear
clc

mission='Lageos';


if isequal(mission,'Topex')
a=7705;
e=0.0010;
I=65.99;
end

if isequal(mission,'GRACE')
a=6820;
e=0.0016;
I=89.02;
end

if isequal(mission,'ERS-1')
a=7156;
e=0.0010;
I=98.6;
end

if isequal(mission,'Lageos')
a=12271;
e=0.0040;
I=109.83;
end

a=a*1000;
%% Constants
mu=3.986004415e14;
ae=6378136.3;
we=7.292115e-5;
g=9.81;
J2=1.082e-3;
ws=1.99096871e-7;
%% Orbital Rates
nb=sqrt(mu/a^3);

Obd=-3/2*nb*(ae/a)^2*J2*cosd(I)/(1-e^2)^(1/2);
wbd=-3/4*nb*(ae/a)^2*J2*(1-5*cosd(I)^2)/(1-e^2)^2;
Mbd=nb*(1-3/4*(ae/a)^2*J2*(1-3*cosd(I)^2)/(1-e^2)^(3/2));
ubd=wbd+Mbd;

%% Periods
Tp=2*pi/nb;
Ta=2*pi/Mbd;
Tn=2*pi/ubd;
TD=2*pi/(we+Obd);
TS=2*pi/(ws+Obd);




