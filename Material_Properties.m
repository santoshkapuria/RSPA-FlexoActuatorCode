 

function [Y_s,nu_s,h,d31,d32,h_t,a,G_a,Y_a,h_a,beta_1,beta_2,V,...
    E_z,Y3_p,Y1_p,Y2_p,nu12_p,alpha1,alpha2,mu19,mu39,mu29,NoD,y1,y2,Nfact,prob]=Material_Properties


%%%%%%%%%%%%%%%% SUBSTRAIT PLATE INPUT DATA%%%%%%%%%%%%%%%%%%%%
Y_s=73e09;        %Youngs Modulus
nu_s=0.17;          % Poisons Ratio
h=150e-6;      % Thickness
%%%%%%%%%%%%%%%%%%%%% ACTUATOR INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
 Y1_p=150e9;%81.3e09;%70e9;%70e9;%% %%%%Youngs Modulus
Y2_p=150e9;%81.3e09;%70e9;%81.3e09;%70e9;%70e9;%150e9;%81.3e09;%70e9;%150e9;%70e9;%81.3e09;%70e9;%81.3e09;%70e9;%81.3e09;%70e9;%81.3e09;%150e9;%70e9;%81.3e09;%70e9;%%150e9;%%  %
Y3_p=150e9;%64.5e09;%70e9;%70e9;%150e9;%70e9;%150e9;%70e9;%150e9;%70e9;%70e9;%70e9;%70e9;%150e9;%70e9;%70e9;%%150e9;%% 
 nu12_p=0.3;      % Poisons 

d31=-47.616e-12;%0;%-123e-12;%-5.2;%0;%-5.2;%;%0;%;-95.231e-12;%0;%-210.2e-12;% % Mechanical Strain Coupling constant
d32=-47.616e-12;%0;%-123e-12;%-5.2;%-47.616e-12;%0;%-123e-12;%-5.2;%%-47.616e-12;%0;%-123e-12;%-123e-12;%-47.616e-12;%-123e-12;%-47.616e-12;%-95.231e-12;%0;%-123e-12;%-210.2e-12;%-47.616e-12;%-5.2;%
mu19=8.5e-6;   %0;%0;%0;%0;%0;%10e-6;%
mu29=100e-6;%0;%0;%0;%0;%0;
mu39=100e-6;%0;%100e-6;%0;%100e-6;%0;%100e-6;%0;%100e-6; %0;%0;%5e-6;%  %0;

h_t=5e-6; %5e-6	8e-6	10e-6	12e-6	15e-6	18e-6	20e-6	22e-6	25e-6	28e-6	30e-6	32e-6	35e-6	38e-6	40e-6	42e-6	45e-6	
           %48	50e-6	52e-6	55e-6	58e-6	60e-6	62e-6	65e-6	68e-6	70e-6	72e-6	75e-6	78e-6	80e-6	82e-6	85e-6	88e-6	
           %90	92e-6	95e-6	98e-6	100e-6	120e-6	130e-6	140e-6	150e-6	160e-6	170e-6	180e-6	190e-6	200e-6
 %thickness of PZT
a=1.5e-3;       %length of PZT
%%%%%%%%%%%%%%%%%%%%% AdDHESIVE INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
G_a=0.5e09;%1.07e9;%1.67e9;%;% % shear modelus of adhesive 
Y_a=1.40e09;%3e9;%4.7e9;% % youngs modulus of adhesive
h_a=40e-6;   %thickness of adhesive
beta_1=h_a/G_a;
beta_2=h_a/Y_a;
V=20;              %Input Voltage
E_z=-V/h_t;      %Electric field
alpha1=4;   %PZT pasted on one side of plate 
alpha2=4;
NoD=700;
y1=0;
y2=a/2;
%xT=(x/a)';
Nfact=[1/a;1;1;1;1];%[1/a;1;1/G_a;1/(G_a*a/2);1/(G_a*(a/2)^2)]; % [1/a;1;1/G_a;1/(G_a*a/2);1/(G_a*(a/2)^2)];% 
prob=2; %2  % 1 for plane stress  %2 for plane strain


end