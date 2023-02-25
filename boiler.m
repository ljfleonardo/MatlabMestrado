% CONTINUOUS SYSTEM DRID!
%Nonlinear third order model tor druIII-doVDcolll8r-riser
%Author K J 1BtrolD 870806
% function [dl,am]=boiler(pow, qfw, tfw, qs, states)

% INPUT pov qtw tfw qs
% OUTPUT dl am
% STATE p Vw xr
% p = states(1);
% Vw = states(2);
% xr = states(3);
%Parameters;
adrum= 20;
vdrum = 40;
Vr= 37;
Vdc= 19;
k= 0.01;
%Initials
p = 7.576;
Vw=13.621;
xr=0.091263;


% p=1; Vw=1;xr=1;
pow=1; qfw=1;tfw=1;qs = qfw;
%Properties if steam abd water ub saturated state
a01= 2.728e6;
a11 = -1.792e4;
a21 = -924.0;
hs = a01+(a11+a21*(p-10))*(p-10);
dhsdp = a11+2*a21*(p-10);

a02 = 55.43;
a12 = 7.136;
a22 = 0.224;
rs = a02+(a12+a22*(p-10))*(p-10);
drsdp = a12+2*a22*(p-10);

a03 = 1.408e6;
a13 = 4.565e4;
a23= -1010.0;
hw = a03+(a13+a23*(p- 10))*(p-10);
dhwdp = a13+2*a23*(p-10);

a04 = 691.35;
a14 = -1.867;
a24 = 0.081;
rw = a04+(a14+a24*(p-10))*(p-10);
drwdp = a14+2*a24*(p-10);

a05 = 311.0;
a15 = 7.822;
a25= -0.32;
ts = a05+(a15+a25*(p-10))*(p-10);
dtsdp = a15+2*a25*(p-10);

%Properties of water in subcritical state
%hd = hw+(a06+a16*(p- 10))*(td-ts)
% dhddp = dhwdp+a16*(td-ts)-(a06+a16*(p-10»*dtsdp
%cp - a06+a16*(p-10)
a06 = 5900;
a16 = 250;
%rd· rw+(a07+a17*(p-10»*(td-ts)
%drddp • drwdp+a17*(td-ts)-(a07+a17*(p-10»*dtsdp
%drddt - a07+a17*(p-10)
%a07 = 2.4
%a17= 0 . 2
hfw = hw+(a06+a16*(p-10)*(tfw-ts));
hc = hs-hw;
hr = xr*hs+(1-xr)*hw;
%Average ateam quality volume ratio
s2 = rs/(xr*(rw-rs));
s3 = 1+xr*(rw/rs-1);
am = rw/(rw-rs)*(1-s2*log(s3));
damdx = rw*s2*(log(s3)/(xr*(rw-rs))-1/s3/rs);
%Drum level
lw = Vw/adrum;
lr = am*Vr/adrum;
dl = lr+lw;

%Equations for derivatives of state variables
Vst = Vdrum - Vw + am*Vr;
Vwt = Vw + Vdc + (1-am)*Vr;
e11 = Vst*(hs*drsdp+rs*dhsdp)+Vwt*(hw*drwdp+rw*dhwdp);
e12 = hw*rw-hs*rs;
e13 = (hs*rs-hw*rw)*Vr*damdx;
b1 = pow*1e6+qfw*hfw-qs*hs;
e21 = Vst*drsdp+Vwt*drwdp;
e22 = rw-rs;
e23 = (rs-rw)*Vr*damdx;
b2 = qfw-qs;
e31 = ((1-xr)*hc*drsdp+rs*dhsdp)*am*Vr+(rw*dhwdp-xr*hc*drwdp)*(1-am)*Vr;
e32 = 0;
e33 = ((1-xr)*rs+xr*rw)*hc*Vr*damdx;
b3 = pow*1e6-qdc*xr*hc;

%Solve linear equation for derivatives of state
p1 = e21/e11;
e221 = e22-e12*p1;
e231 = e23-e13*p1;
b21 = b2-b1*p1;

p2 = e31/e11;
e321 = -e12*p2;
e331 = e33-e13*p2;
b31 = b3-b1*p2;
p3 = e321/e221;
e332 = e331-e231*p3;
b32 = b31-b21*p3;

dxr = b32/e332;
dVw = (b21-e231*dxr)/e221;
dp = (b1-e12*dVw-e13*dxr)/e11;

%Circulation flov
s1 = 2*(rw-rs)*Vr*am/k;
qdc = sqrt(s1);
qr= qdc-(am*drsdp+(1-am)*drwdp)*Vr*dp+(rw-rs)*Vr*damdx*dxr;
%Total condensation flov
qc = (rs*Vst*dhsdp+rw*Vwt*dhwdp)*dp/hc;
%Condensation flov in risera
qcr = (rs*am*Vr*dhsdp+rw*(1-am)*Vr*dhsdp)*dp/hc;






%pow Pover from fuel [MW]
%qfw Feedvater flov [ltg/a]
%tfw Feedvater temperature [deg C]
%qs Steam flov [ltg/s]
%dl Drum level [m]
%am Steam quality volume ratio
%qc Condensate flov (total) [ltg/s]
%qcr Condensate flov (risers) [kg/a]
%p Drum pressure [MPa]
%Vv Drum vater volume [m*m*m]
%xr Steam quality at riser outlet
%Properties of steam and vater in saturated state









%Properties of vater in subcritical state













% end