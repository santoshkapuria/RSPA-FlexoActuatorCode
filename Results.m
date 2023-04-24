function[tauT,sigma_zT,XT,Shearforce,Moment,MxtT,kappa_tT,w_tT,w_sT,MxsT,QxtT]=Results(alpha2,s11bar_t,s11bar_s,beta_1,beta_2,h_t,alpha1,h,lambda19_bar,d31bar,E_z,a,NoD,y1,y2,Nfact)


Gamma_11=(alpha2*s11bar_t/(beta_1*h_t))+(alpha1*s11bar_s/(beta_1*h));
Gamma_12=(6*s11bar_t/(beta_1*h_t^2))-(6*s11bar_s/(beta_1*h^2));
Gamma_22=(12/beta_2)*((s11bar_s/(h^3))+(s11bar_t/(h_t^3)));
Gamma_21=(6/beta_2)*((s11bar_t/h_t^2)-(s11bar_s/h^2));
Gamma=(Gamma_22*Gamma_11-(Gamma_12*Gamma_21));
P=[1 0 -Gamma_11 0 Gamma_22 0 -Gamma 0];
r=roots(P);

%%%% Finding the Roots


e=-1;  
b=-Gamma_11;
c=Gamma_22;
d=-Gamma;

F=-(2*b^3+9*e*b*c+27*d*e^2)/(27*e^3);
E=-(b^2+3*c*e)/(3*e^2);
G=(F^2/4)+(E^3/27);

w=(G^(1/2)-F/2)^(1/3);
del=(G^(1/2)+F/2)^(1/3);
%del=nthroot(G^(1/2)+F/2,3);
LAMBDA=w-del;
p=(w-del+(b/(3*e)))^(1/2);
alpha_j=(-LAMBDA/2)+b/(3*e);
beta_j=3^(1/2)*(w+del)/2;
theta=atan(beta_j/alpha_j)/2;
alpha=(alpha_j^2+beta_j^2)^(1/4)*cos(theta);
gamma=(alpha_j^2+beta_j^2)^(1/4)*sin(theta);


%%%% Coefficients of linear equations

delta_1=(-p^3+Gamma_11*p)/Gamma_12;
delta_2=(-alpha^3+alpha*Gamma_11+3*alpha*gamma^2)/Gamma_12;
delta_3=(-gamma^3-gamma*Gamma_11+3*alpha^2*gamma)/Gamma_12;

delta_11=p;
delta_12=alpha;
delta_13=-gamma;
delta_21=p^2*delta_1;
delta_22=delta_2*(alpha^2-gamma^2)+2*delta_3*alpha*gamma;%(alpha^3-3*alpha*beta^2-alpha*gamma_11)*(alpha^2-beta^2)/gamma_12+(beta^3-3*alpha^2*beta+beta*gamma_11)*(2*alpha*beta)/gamma_12;
delta_23=-2*delta_2*alpha*gamma+delta_3*(alpha^2-gamma^2);%(alpha^3-3*alpha*beta^2-alpha*gamma_11)*(-2*alpha*beta)/gamma_12+(beta^3-3*alpha^2*beta+beta*gamma_11)*(alpha^2-beta^2)/gamma_12;
delta_31=p^3*delta_1-Gamma_21;%((p^6-(gamma_11*p^4))/gamma_12)-gamma_21;  %-2.2929e12;
delta_32=delta_2*alpha*(alpha^2-3*gamma^2)+delta_3*gamma*(3*alpha^2-gamma^2)-Gamma_21;%((alpha^6-beta^6+15*alpha^2*beta^4-15*alpha^4*beta^2)/gamma_12-gamma_11*(alpha^4+beta^4-6*alpha^2*beta^2)/gamma_12)-gamma_21; %4.7482e15;%
delta_33=delta_2*gamma*(gamma^2-3*alpha^2)+delta_3*alpha*(alpha^2-3*gamma^2);%(20*alpha^3*beta^3-6*alpha*beta^5-6*alpha^5*beta)/gamma_12-gamma_11*(4*alpha*beta^3-4*alpha^3*beta)/gamma_12; %3.6446e15;%
% delta_11=p;
% delta_12=alpha;
% delta_13=-gamma;
% delta_21=(p^5-p^3*Gamma_11)/Gamma_12;
% delta_22=(alpha^3-3*alpha*gamma^2-alpha*Gamma_11)*(alpha^2-gamma^2)/Gamma_12+(gamma^3-3*alpha^2*gamma+gamma*Gamma_11)*(2*alpha*gamma)/Gamma_12;
% delta_23=(alpha^3-3*alpha*gamma^2-alpha*Gamma_11)*(-2*alpha*gamma)/Gamma_12+(gamma^3-3*alpha^2*gamma+gamma*Gamma_11)*(alpha^2-gamma^2)/Gamma_12;
% delta_31=((p^6-(Gamma_11*p^4))/Gamma_12)-Gamma_21;  %-2.2929e12;
% delta_32=((alpha^6-gamma^6+15*alpha^2*gamma^4-15*alpha^4*gamma^2)/Gamma_12-Gamma_11*(alpha^4+gamma^4-6*alpha^2*gamma^2)/Gamma_12)-Gamma_21; %4.7482e15;%
% delta_33=(20*alpha^3*gamma^3-6*alpha*gamma^5-6*alpha^5*gamma)/Gamma_12-Gamma_11*(4*alpha*gamma^3-4*alpha^3*gamma)/Gamma_12; %3.6446e15;%

% delta_1=(p^3-Gamma_11*p)/Gamma_12
% delta_2=(alpha^3-alpha*Gamma_11-3*alpha*gamma^2)/Gamma_12
% delta_3=(gamma^3+gamma*Gamma_11-3*alpha^2*gamma)/Gamma_12


Q1t=delta_1/p;
Q2t=(alpha*delta_2-gamma*delta_3)/(alpha^2+gamma^2);
Q3t=(alpha*delta_3+gamma*delta_2)/(alpha^2+gamma^2);

M1t=(2*delta_1-h_t*p)/(2*p^2);
M2t=(2*delta_2*(alpha^2-gamma^2)-alpha*(h_t*(alpha^2+gamma^2)+4*delta_3*gamma))/(2*(alpha^2+gamma^2)^2);
M3t=(2*delta_3*(alpha^2-gamma^2)-gamma*(h_t*(alpha^2+gamma^2)-4*delta_2*alpha))/(2*(alpha^2+gamma^2)^2);

% w1tt=-6*s11bar_t*(2*delta_1-h_t*p)/(p^4*h_t^3)
% w2tt=-(6*s11bar_t/(h_t^3*(alpha^2+gamma^2)^4))*(2*delta_2*(alpha^4+gamma^4-6*alpha^2*gamma^2)-8*delta_3*alpha*gamma*(alpha^2-gamma^2)-alpha*h_t*(alpha^4-gamma^2*(2*alpha^2+3*gamma^2)))
% w3tt=-(6*s11bar_t/(h_t^3*(alpha^2+gamma^2)^4))*(2*delta_3*(alpha^4+gamma^4-6*alpha^2*gamma^2)+8*delta_2*alpha*gamma*(alpha^2-gamma^2)+gamma*h_t*(gamma^4-alpha^2*(2*gamma^2+3*alpha^2)))
% % 

w1t=-12*s11bar_t*M1t/(p^2*h_t^3);
w2t=(-12*s11bar_t/h_t^3)*(M2t*(alpha^2-gamma^2)-M3t*2*alpha*gamma)/((alpha^2+gamma^2)^2);
w3t=(-12*s11bar_t/h_t^3)*(M3t*(alpha^2-gamma^2)+M2t*2*alpha*gamma)/((alpha^2+gamma^2)^2);


M1s=(2*delta_1+h*p)/(2*p^2);
M2s=(2*delta_2*(alpha^2-gamma^2)+alpha*(h*(alpha^2+gamma^2)-4*delta_3*gamma))/(2*(alpha^2+gamma^2)^2);
M3s=(2*delta_3*(alpha^2-gamma^2)+gamma*(h*(alpha^2+gamma^2)+4*delta_2*gamma))/(2*(alpha^2+gamma^2)^2);
% 
% S1=(3*(2*delta_1-h_t*p))/(p^2*h_t^2);
% S2=(3*(2*delta_2*alpha^2-2*delta_2*gamma^2-alpha^3*h_t-alpha*gamma^2*h_t-4*delta_3*alpha*gamma))/(h_t^2*(alpha^2+gamma^2)^2);
% S3=(3*(2*delta_3*alpha^2-2*delta_3*gamma^2-gamma^3*h_t-alpha^2*gamma*h_t+4*delta_2*alpha*gamma))/(h_t^2*(alpha^2+gamma^2)^2);
% 
% S11=1/(h_t*p);
% S12=alpha/(h_t*(alpha^2+gamma^2));
% S13=gamma/(h_t*(alpha^2+gamma^2));

% M1=(-2*gamma_tau1-gamma_sig1*h_t+2*lambda_1^2)/(gamma_sig1*lambda_1*2);
% M2=(2*alpha^3-2*alpha*gamma_tau1+2*alpha*beta^2-alpha*gamma_sig1*h_t)/(gamma_sig1*2*(alpha^2+beta^2));
% M3=(2*beta^3+2*beta*gamma_tau1+2*alpha^2*beta+beta*gamma_sig1*h_t)/(gamma_sig1*2*(alpha^2+beta^2));


V11=1/p;
V12=alpha/(alpha^2+gamma^2);
V13=gamma/(alpha^2+gamma^2);

M11=-1/(p^2);
M12=(delta_2*gamma^2-delta_2*alpha^2+2*delta_3*gamma*alpha)/((alpha^2+gamma^2)^2);
M13=(delta_3*alpha^2-delta_3*gamma^2+2*delta_2*gamma*alpha)/((alpha^2+gamma^2)^2);
M21=(delta_2*gamma+delta_3*alpha)/(alpha^2+gamma^2);
M23=(delta_2*alpha-delta_3*gamma)/(alpha^2+gamma^2);
% 
% M1h=(2*delta_1+h*p)/(2*p^2);
% M2h=(2*delta_2*alpha^2-2*delta_2*gamma^2+alpha^3*h+alpha*gamma^2*h-4*delta_3*alpha*gamma)/(2*(alpha^2+gamma^2)^2);
% M3h=(2*delta_3*alpha^2-2*delta_3*gamma^2+gamma^3*h+alpha^2*gamma*h+4*delta_2*alpha*gamma)/(2*(alpha^2+gamma^2)^2);
% 
% N1h=-1/p;
% N2h=-alpha/(alpha^2+gamma^2);
% N3h=-gamma/(alpha^2+gamma^2);

% u11=s11bar_t/((p^2)*h_t);
% u12=s11bar_t*(alpha^2-gamma^2)/(h_t*(alpha^2+gamma^2)^2);
% u13=(s11bar_t*2*alpha*gamma)/(h_t*(alpha^2+gamma^2)^2);
% 
% w11=2*S1*s11bar_t/(h_t*p);
% w12=2*s11bar_t*(S2*alpha-S3*gamma)/(h_t*(alpha^2+gamma^2));
% w13=2*s11bar_t*(S3*alpha+S2*gamma)/(h_t*(alpha^2+gamma^2));

% w21=2*S1*s11bar_t/(p^2*h_t);
% w22=2*s11bar_t*S2*(alpha^2-gamma^2)/(h_t*(alpha^2+gamma^2)^2);
% w23=2*s11bar_t*S3*(alpha^2-gamma^2)/(h_t*(alpha^2+gamma^2)^2);
% 
% k1=2*s11bar_t*S1/h_t;
% k2=2*s11bar_t*S2/h_t;
% k3=2*s11bar_t*S3/h_t;

ch=cosh(p*a/2);
sh=sinh(p*a/2);
cc=cosh(alpha*a/2)*cos(gamma*a/2);
ss=sinh(alpha*a/2)*sin(gamma*a/2);
cs=cosh(alpha*a/2)*sin(gamma*a/2);
sc=sinh(alpha*a/2)*cos(gamma*a/2);

A11=delta_11*ch;
A12=delta_12*cc+delta_13*ss;
A13=delta_12*ss-delta_13*cc;
A21=delta_21*ch;
A22=delta_22*cc+delta_23*ss;
A23=delta_22*ss-delta_23*cc;
A31=delta_31*sh;
A32=delta_32*sc+delta_33*cs;
A33=delta_32*cs-delta_33*sc;
 A=[A11 A12 A13;A21 A22 A23;A31 A32 A33];

D1=(6*lambda19_bar/(beta_1*h_t))+d31bar*Nfact(2)/beta_1 ; %%1.5897e+08;  %2.4087e+08;  %%%3.9743e+07;%;%+
D2=(12*lambda19_bar/(h_t^2*beta_2));
D3=0;
 D=[D1;D2;D3];

C=A\(E_z*D);
%  tau=C(1)*sinh(lambda_1*x)+C(2)*sinh(alpha*x)*cos(beta*x)+C(3)*cosh(alpha*x)*sin(beta*x);
% sigma=C(1)*D11*cosh(lambda_1*x)+C(2)*(D12*cosh(alpha*x)*cos(beta*x)+D13*sinh(alpha*x)*sin(beta*x))+C(3)*(D12*sinh(alpha*x)*sin(beta*x)-D13*cosh(alpha*x)*cos(beta*x));

 
 
 

 %%%%%%% expression for curvature

Q0t=C(1)*Q1t*sh+C(2)*(sc*Q2t+cs*Q3t)+C(3)*(Q2t*cs-Q3t*sc);

M0t=C(1)*M1t*ch+C(2)*(M2t*cc+M3t*ss)+C(3)*(M2t*ss-M3t*cc)-a*Q0t/2;

Q0s=-Q0t;

M0s=-C(1)*M1s*ch-C(2)*(M2s*cc+M3s*ss)-C(3)*(M2s*ss-M3s*cc)-a*Q0s/2;

% 
% C5h=-C(1)*Q1t*sh-C(2)*(sc*Q2t+cs*Q3t)+C(3)*(-Q2t*cs+Q3t*sc);
% 
% C6=C(1)*S1*ch+C(2)*(cc*S2+ss*S3)+C(3)*(S2*ss-S3*cc);
% 
% %C7=C(1)*M1*ch+C(2)*(cc*M2-ss*M3)+C(3)*(cc*M3+ss*M2);
% C8=C(1)*S11*ch+C(2)*(cc*S12+ss*S13)+C(3)*(ss*S12-cc*S13);
% 
% C9h=-C(1)*M1h*ch-C(2)*(M2h*cc+M3h*ss)-C(3)*(M2h*ss-M3h*cc)-C5h*a/2;
% 
%  C10h=C(1)*N1h*ch+C(2)*(N2h*cc+N3h*ss)+C(3)*(N2h*ss-N3h*cc);

x=y1:((y2-y1)/NoD):y2;
 
 for i=1:length(x)

tau(i)=C(1)*sinh(p*x(i))+C(2)*sinh(alpha*x(i))*cos(gamma*x(i))+C(3)*cosh(alpha*x(i))*sin(gamma*x(i));
    
sigma_z(i)=C(1)*delta_1*cosh(p*x(i))+C(2)*(delta_2*cosh(alpha*x(i))*cos(gamma*x(i))+delta_3*sinh(alpha*x(i))*sin(gamma*x(i)))+C(3)*(delta_2*sinh(alpha*x(i))*sin(gamma*x(i))-delta_3*cosh(alpha*x(i))*cos(gamma*x(i)));
 
Qxt(i)=C(1)*Q1t*sinh(p*x(i))+C(2)*(sinh(alpha*x(i))*cos(gamma*x(i))*Q2t+cosh(alpha*x(i))*sin(gamma*x(i))*Q3t)+C(3)*(Q2t*cosh(alpha*x(i))*sin(gamma*x(i))-Q3t*sinh(alpha*x(i))*cos(gamma*x(i)))-Q0t;

Mxt(i)=C(1)*M1t*cosh(p*x(i))+C(2)*(M2t*cosh(alpha*x(i))*cos(gamma*x(i))+M3t*sinh(alpha*x(i))*sin(gamma*x(i)))+C(3)*(M2t*sinh(alpha*x(i))*sin(gamma*x(i))-M3t*cosh(alpha*x(i))*cos(gamma*x(i)))-Q0t*x(i)-M0t;

kappa_t(i)=-12*((Mxt(i)*s11bar_t/(h_t^3))-(lambda19_bar*E_z/(h_t^2)));

w_t(i)=C(1)*w1t*cosh(p*x(i))+C(2)*(w2t*cosh(alpha*x(i))*cos(gamma*x(i))+w3t*sinh(alpha*x(i))*sin(gamma*x(i)))+C(3)*(w2t*sinh(alpha*x(i))*sin(gamma*x(i))-w3t*cosh(alpha*x(i))*cos(gamma*x(i)))+((Q0t*2*s11bar_t*x(i)^3)/(h_t^3))+6*x(i)^2*(M0t*s11bar_t+E_z*lambda19_bar*h_t)/(h_t^3);

w_s(i)=w_t(i)-(beta_2*sigma_z(i));
%G(i)=((Q0t*2*s11bar_t*x(i)^3)/(h_t^3))-6*x(i)^2*(M0t*s11bar_t+E_z*lambda19_bar*h_t)/(h_t^3);
Mxs(i)=-C(1)*M1s*cosh(p*x(i))-C(2)*(M2s*cosh(alpha*x(i))*cos(gamma*x(i))+M3s*sinh(alpha*x(i))*sin(gamma*x(i)))-C(3)*(M2s*sinh(alpha*x(i))*sin(gamma*x(i))-M3s*cosh(alpha*x(i))*cos(gamma*x(i)))-Q0s*x(i)-M0s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qh(i)=-C(1)*Q1t*sinh(p*x(i))-C(2)*(sinh(alpha*x(i))*cos(gamma*x(i))*Q2t+cosh(alpha*x(i))*sin(gamma*x(i))*Q3t)+C(3)*(-Q2t*cosh(alpha*x(i))*sin(gamma*x(i))+Q3t*sinh(alpha*x(i))*cos(gamma*x(i)))-C5h;
% 
% sigmaxa(i)=C(1)*S11*cosh(p*x(i))+C(2)*(S12*cosh(alpha*x(i))*cos(gamma*x(i))+S13*sinh(alpha*x(i))*sin(gamma*x(i)))+C(3)*(S12*sinh(alpha*x(i))*sin(gamma*x(i))-S13*cosh(alpha*x(i))*cos(gamma*x(i)))-C8;
% 
% sigmaxb(i)=C(1)*S1*cosh(p*x(i))+C(2)*(cosh(alpha*x(i))*cos(gamma*x(i))*S2+sinh(alpha*x(i))*sin(gamma*x(i))*S3)+C(3)*(sinh(alpha*x(i))*sin(gamma*x(i))*S2-cosh(alpha*x(i))*cos(gamma*x(i))*S3)-C6;
% 
% K(i)=(sigmaxb(i)-(6*lambda19_bar*E_z/(h_t*s11bar_t)))*2*s11bar_t/h_t;
% 
% k(i)=C(1)*k1*cosh(p*x(i))+C(2)*(k3*sinh(alpha*x(i))*sin(gamma*x(i))+k2*cosh(alpha*x(i))*cos(gamma*x(i)))+C(3)*(k2*sinh(alpha*x(i))*sin(gamma*x(i))-k3*cosh(alpha*x(i))*cos(gamma*x(i)))-2*C6*s11bar_t/h_t-(12*E_z*lambda19_bar/(h_t^2));
% 
% %  Mx(i)=C(1)*M1*cosh(lambda_1*x(i))+C(2)*(cosh(alpha*x(i))*cos(beta*x(i))*M2-sinh(alpha*x(i))*sin(beta*x(i))*M3)+C(3)*(cosh(alpha*x(i))*cos(beta*x(i))*M3+sinh(alpha*x(i))*sin(beta*x(i))*M2)-1.7764e-15*x(i)-C7;
% % 
% %   k(i)=(Mx(i)-(mubar*E_z*h_t/s11bar_t))*12*s11bar_t/(h_t^3);
% 
%  % ut(i)=C(1)*u11*sinh(p*x(i))+C(2)*(u12*cos(gamma*x(i))*sinh(alpha*x(i))+u13*cosh(alpha*x(i))*sin(gamma*x(i)))+C(3)*(u12*cosh(alpha*x(i))*sin(gamma*x(i))-u13*cos(gamma*x(i))*sinh(alpha*x(i)))+(E_z*d31bar-(C8*s11bar_t))*x(i);
% 
%  
% 
%   %wtx(i)=C(1)*w11*sinh(p*x(i))+C(2)*(w12*cos(gamma*x(i))*sinh(alpha*x(i))+w13*cosh(alpha*x(i))*sin(gamma*x(i)))+C(3)*(w12*cosh(alpha*x(i))*sin(gamma*x(i))-w13*cos(gamma*x(i))*sinh(alpha*x(i)))-(x(i)*(12*E_z*mubar+(2*C6*h_t*s11bar_t))/(h_t^2));
%   
%  wt(i)=C(1)*w21*cosh(p*x(i))+C(2)*(w22*cosh(alpha*x(i))*cos(gamma*x(i))+w23*sinh(alpha*x(i))*sin(gamma*x(i)))+C(3)*(w22*sinh(alpha*x(i))*sin(gamma*x(i))-w23*cosh(alpha*x(i))*cos(gamma*x(i)))-(x(i)^2*(6*E_z*lambda19_bar+C6*h_t*s11bar_t)/(h_t^2));
%  
%  ws(i)=wt(i)-(sigma(i)*h_a/Y_a);
% 
%  Nx(i)=h_t*sigmaxa(i);
% 
%  Mx(i)=((h_t^2)/6)*sigmaxb(i);
% 
%  Mxh(i)=-C(1)*M1h*cosh(p*x(i))-C(2)*(M2h*cosh(alpha*x(i))*cos(gamma*x(i))+M3h*sinh(alpha*x(i))*sin(gamma*x(i)))-C(3)*(M2h*sinh(alpha*x(i))*sin(gamma*x(i))-M3h*cosh(alpha*x(i))*cos(gamma*x(i)))-C5h*x(i)-C9h;
%  
% 
%  Nxh(i)=C(1)*N1h*cosh(p*x(i))+C(2)*(N2h*cosh(alpha*x(i))*cos(gamma*x(i))+N3h*sinh(alpha*x(i))*sin(gamma*x(i)))+C(3)*(N2h*sinh(alpha*x(i))*sin(gamma*x(i))-N3h*cosh(alpha*x(i))*cos(gamma*x(i)))-C10h;
 end
XT=x'*Nfact(1);
 tauT=tau'*Nfact(3);
sigma_zT=sigma_z'*Nfact(3);
  
QxtT=Qxt';
 MxtT=Mxt';
 kappa_tT=kappa_t';
 w_tT=w_t';
 w_sT=w_s';
 MxsT=Mxs'*Nfact(5);

%GT=G';

%  figure(5)
%  plot(x/a,-Nxh*Nfact(4))
% 
% NxhT=-Nxh'*Nfact(4);
% W=ws';
 %Ws=h*ws'/(h_t*a/2);

%%%%%% Evaluation of Shear force and Moment
y=y1;
v1=C(1)*V11*cosh(p*y)+C(2)*(V12*cosh(alpha*y)*cos(gamma*y)+V13*sinh(alpha*y)*sin(gamma*y))+C(3)*(V12*sinh(alpha*y)*sin(gamma*y)-V13*cosh(alpha*y)*cos(gamma*y));
m1=C(1)*M11*delta_1*(cosh(p*y)-p*y*sinh(p*y))+C(2)*(M12*cosh(alpha*y)*cos(gamma*y)-M13*sinh(alpha*y)*sin(gamma*y))...
    +C(2)*(M21*y*cosh(alpha*y)*sin(gamma*y)+M23*y*cos(gamma*y)*sinh(alpha*y))+...
    C(3)*(M12*sinh(alpha*y)*sin(gamma*y)+M13*cosh(alpha*y)*cos(gamma*y))...
    +C(3)*(-M21*y*cos(gamma*y)*sinh(alpha*y)+M23*y*cosh(alpha*y)*sin(gamma*y));
y=y2;
v2=C(1)*V11*cosh(p*y)+C(2)*(V12*cosh(alpha*y)*cos(gamma*y)+V13*sinh(alpha*y)*sin(gamma*y))+C(3)*(V12*sinh(alpha*y)*sin(gamma*y)-V13*cosh(alpha*y)*cos(gamma*y));
m2=C(1)*M11*delta_1*(cosh(p*y)-p*y*sinh(p*y))+C(2)*(M12*cosh(alpha*y)*cos(gamma*y)-M13*sinh(alpha*y)*sin(gamma*y))+C(2)*(M21*y*cosh(alpha*y)*sin(gamma*y)+M23*y*cos(gamma*y)*sinh(alpha*y))+...
    C(3)*(M12*sinh(alpha*y)*sin(gamma*y)+M13*cosh(alpha*y)*cos(gamma*y))+C(3)*(-M21*y*cos(gamma*y)*sinh(alpha*y)+M23*y*cosh(alpha*y)*sin(gamma*y));


Shearforce=v2-v1;
Moment=m2-m1;
end

%%%%%% if G vlaue is negative 

% if G<=0
% 
% for i=1:3
% L=-(27)^(1/2)*F/(2*(-E)^(3/2));
% 
% M(i)=cos((acos(L)+2*pi*(i-1))/3);
% 
% Gam(i)=2*(-E/3)^(1/2)*M(i) ;
% 
% alpha(i)=abs((Gam(i)+(b/(3*e)))^(1/2)) 
% 
% end
% A1=alpha(1)*cosh(alpha(1)*a/2);
% A2=alpha(2)*cos(alpha(2)*a/2);
% A3=alpha(3)*cos(alpha(3)*a/2);
% B1=((alpha(1)^5-(gamma_tau1*alpha(1)^3))/gamma_sig1)*cosh(alpha(1)*a/2);
% B2=((alpha(2)^5+(gamma_tau1*alpha(2)^3))/gamma_sig1)*cos(alpha(2)*a/2);
% B3=((alpha(3)^5+(gamma_tau1*alpha(3)^3))/gamma_sig1)*cos(alpha(3)*a/2);
% C1=((alpha(1)^6-(gamma_tau1*alpha(1)^4))/gamma_sig1)*sinh(alpha(1)*a/2);
% C2=-((alpha(2)^6+(gamma_tau1*alpha(2)^4))/gamma_sig1)*sin(alpha(2)*a/2);
% C3=-((alpha(3)^6+(gamma_tau1*alpha(3)^4))/gamma_sig1)*sin(alpha(3)*a/2);
% 
% D11=(alpha(1)^3-(gamma_tau1*alpha(1)))/gamma_sig1;
% D12=(alpha(2)^3+(gamma_tau1*alpha(2)))/gamma_sig1;
% D13=(alpha(3)^3+(gamma_tau1*alpha(3)))/gamma_sig1;
% Coeff=[A1 A2 A3;B1 B2 B3;C1 C2 C3];
% D1=(6*mubar/(phi*h_t))+d31bar/phi;%%1.5897e+08;  %2.4087e+08;  %%%3.9743e+07;%;%+
% D2=(12*mubar/h_t^2)*(Y_a/h_a);
% D3=0;
%  D=[D1;D2;D3];
% 
%  C=Coeff\(E_z*D);
% 
%  x=0:0.00001:a/2;
%  for i=1:length(x)
% tau(i)=C(1)*sinh(alpha(1)*x(i))+C(2)*sin(alpha(2)*x(i))+C(3)*sin(alpha(3)*x(i));
% 
% sigma(i)=C(1)*D11*cosh(alpha(1)*x(i))-C(2)*D12*cos(alpha(2)*x(i))-C(3)*D13*cos(alpha(3)*x(i));
% 
%  end
%  %  XT=x'/a;
% %  tauT=tau';
% %  sigmaT=sigma';
%  figure(1)
%  plot(x/a,tau)
%  figure(2)
%  plot(x/a,sigma)
% 
% else


