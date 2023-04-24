clc
clear all


[Y_s,nu_s,h,d31,d32,h_t,a,G_a,Y_a,h_a,beta_1,beta_2,V,...
    E_z,Y3_p,Y1_p,Y2_p,nu12_p,alpha1,alpha2,mu19,mu39,mu29,NoD,y1,y2,Nfact,prob]=Material_Properties;

%prob=2; % 1 for plane stress  %2 for plane strain

if prob==2
[s11bar_s,s11bar_t,d31bar,lambda19_bar]=Plane_strain(nu_s,Y_s,Y1_p,nu12_p,Y2_p,Y3_p,d31,mu19,mu39,d32);
else

[s11bar_s,s11bar_t,d31bar,lambda19_bar]=Plane_stress(Y_s,Y1_p,Y3_p,Y2_p,d31,nu12_p,mu19,mu29,mu39);
end

[tauT,sigma_zT,XT,Shearforce,Moment,MxtT,kappa_tT,w_tT,w_sT,MxsT,QxtT]=Results(alpha2,s11bar_t,s11bar_s,beta_1,beta_2,h_t,alpha1,h,lambda19_bar,d31bar,E_z,a,NoD,y1,y2,Nfact);
%d31bar,tauT,sigmaT,XTE_z,l