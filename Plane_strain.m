function[s11bar_s,s11bar_t,d31bar,lambda19_bar]=Plane_strain(nu_s,Y_s,Y1_p,nu12_p,Y2_p,Y3_p,d31,mu19,mu39,d32)
%Y1_p,nu12_p,Y2_p,Y3_p,
s11bar_s=(1-nu_s^2)/Y_s;
s11p=1/Y1_p; %8.05e-12;  %1.22e-12;%1/Y1_p;12.3e-12;% 
s12p=-nu12_p/Y2_p; %-2.35e-12; %-10.01e-12;%-4.05e-12;%
s22p=1/Y2_p; %8.05e-12; %1.22e-12;%12.3e-12;% 
s13p=-nu12_p/Y3_p; %-5.24e-12; %9.22e-12;%-5.31e-12;%
s23p=-nu12_p/Y3_p; %-5.24e-12; %9.22e-12;%-5.31e-12;%-5.31e-12;%
s11bar_t=s11p-s12p^2/s22p;            
d31bar=d31-d32*s12p/s22p;
lambda19_bar=-(mu19*(s11p-s12p^2/s22p)+mu39*(s13p-s23p*s12p/s22p));
end