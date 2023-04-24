function[s11bar_s,s11bar_t,d31bar,lambda19_bar]=Plane_stress(Y_s,Y1_p,Y3_p,Y2_p,d31,nu12_p,mu19,mu29,mu39)

s11bar_s=1/Y_s;
s11p=1/Y1_p;
s11bar_t=s11p;
s12p=-nu12_p/Y2_p;
d31bar=d31;
s13p=-nu12_p/Y3_p;
lambda19_bar=-(s11p*mu19+s12p*mu29+s13p*mu39);
end