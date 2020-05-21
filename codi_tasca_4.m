% Aquest codi ha estat realitzat per:
% -Alexis Leon Delgado
% -David Morante Torra
% -Ferran Rubio Vallhonrat
% -Juan Garrido Moreno

% El codi determina l'impuls i empenta específics a partir del set de
% paràmetres de la tasca 3, imposant tovera adaptada. Finalment, determina
% el guany relatiu i la relació d'àrees que ho fa possible.

clear all;
close all;
clc;

%Introducció de les dades del motor a analitzar. Totes elles expressades
%segons el sistema internacional d'unitats.
M_0=0.85; z=11000; T_SL=288.15; P_SL=101325;
pi_d=0.98; eta_f=0.89; eta_lpc=0.88; eta_hpc=0.86;
pi_b=0.99; eta_b=0.99; eta_hpt=0.91; eta_lpt=0.92;
eta_mH=0.993; eta_mL=0.997; pi_np=0.99; pi_ns=0.99;
T_t4=1450; gamma_c=1.4; gamma_t=1.3; R_g=287; h=43000000; g=9.81;

%% CÀLCULS PREVIS

%Es calculen els coeficients de pressió per l'aire del compressor (Cp_c) i
%de la turbina (Cp_t), segons les gammes depenent de si ha passat o no per
%la cambra de combustió.
Cp_t=(gamma_t/(gamma_t-1))*R_g;
Cp_c=(gamma_c/(gamma_c-1))*R_g;

%S'obtenen també els valors de pressió i temperatura del flux lliure a
%l'altura de funcionament. Amb aquests es calula la velocitat d'entrada de
%l'aire al motor.
[T_0,P_0]=atmosferaISA(z,T_SL,P_SL);
u_0=M_0*sqrt(gamma_c*R_g*T_0);

%% Set de paràmetres

%Es proposa un grup de paràmetres que dóna relativament bons resultats per
%l'impuls i l'empenta específics.
alpha=10;
pi_f=1.5737;
pi_lpc=5.3;
pi_hpc=8.5;
OPR=pi_lpc*pi_hpc;

%Per tal de fer l'anàlisi del motor, es divideix el procés en les etapes
%que es troba l'aire al seu pas aquest. La notació emprada correspon a
%la que apareix a l'esquema del motor.

%% (0)-(2): INLET

theta_0=1+0.5*(gamma_c-1)*M_0^2; % theta_0=T_t0/T_0
delta_0=theta_0^(gamma_c/(gamma_c-1)); % delta_0=P_t0/P_0

%Es calcula la temperatura d'estancament del flux lliure, que roman
%constant durant la presa d'aire.
T_t0=T_0*theta_0;
T_t2=T_t0;

%Es calculen també les pressions d'estancament. Aquesta és menor a la presa
%d'aire degut a l'eficiència d'aquesta.
P_t0=P_0*delta_0;
P_t2=P_t0*pi_d;

%% (2)-(2.5): LPC

%Es calcula la pressió al compressor de baixa pressió gràcies al rendiment
%isentròpic d'aquest. A partir d'aquest, i de la gamma de l'aire s'obté
%l'eficiència isentròpica que permet conèixer la temperatura.
P_t25=P_t2*pi_lpc;
tau_lpc=1+(pi_lpc^((gamma_c-1)/gamma_c)-1)/eta_lpc;
T_t25=T_t2*tau_lpc;

%% (2.5)-(3): HPC

%El procediment de càlcul al compressor d'alta pressió és el mateix que al
%de baixa, partint de la temperatura i la pressió d'aquest últim.
P_t3=P_t25*pi_hpc;
tau_hpc=1+(pi_hpc^((gamma_c-1)/gamma_c)-1)/eta_hpc;
T_t3=T_t25*tau_hpc;

%% (3)-(4): Cambra de Combustió

%Amb les propietats físiques de l'aire a l'entrada i la sortida d'aquesta i
%l'eficiència isentròpica es pot conèixer la relació de cabal de
%combustible respecte el cabal d'aire. I amb la relació de pressions
%d'estancament, la pressió a la sortida.
f=(Cp_t*T_t4-Cp_c*T_t3)/(eta_b*h-Cp_t*T_t4);
P_t4=P_t3*pi_b;

%% (4)-(4.5): HPT

%Per obtenir la relació de temperatures a la turbina d'alta pressió s'usen
%l'eficiència mecànica en aquesta, així com les temperatures al compressor
%i a l'entrada de la HPT, tenint en compte la presència de combustible. A
%partir d'aquesta relació, ja es poden calcular la relació de pressions i
%la temperatura i pressió a la sortida d'aquest element.
tau_hpt=1-(1/eta_mH)*(Cp_c/Cp_t)*(1/(1+f))*(T_t25/T_t4)*(tau_hpc-1);
pi_hpt=(1-(1-tau_hpt)/eta_hpt)^(gamma_t/(gamma_t-1));
T_t45=T_t4*tau_hpt;
P_t45=P_t4*pi_hpt;

%% (4.5)-(5): LPT

%La turbina de baixa pressió segueix el mateix concepte i procediment que
%la turbina d'alta pressió.
tau_f=1+(pi_f^((gamma_c-1)/gamma_c)-1)/eta_f;
tau_lpt=1-(1/eta_mL)*(Cp_c/Cp_t)*(1/(1+f))*(T_t0/T_t4)*(1/tau_hpt)*((tau_lpc-1)+alpha*(tau_f-1));
pi_lpt=(1-(1-tau_lpt)/eta_lpt)^(gamma_t/(gamma_t-1));
T_t5=T_t45*tau_lpt;
P_t5=P_t45*pi_lpt;

%% (5)-(9): PRIMARY NOZZLE

%Les temperatures d'estancament de la tovera primària es mantenen
%constants, mentre que la pressió disminueix lleugerament a causa de
%l'eficiència de la tovera.
P_t9=P_t5*pi_np;
T_t9=T_t5;

%Relació de pressió d'estancament a la tovera i a l'entrada del motor.
pi_primari_total=P_t9/P_0;  

% Es realitza la hipòtesi de tovera adaptada, però com prèviament s'ha
% comprovat que certa, per la tovera primària no s'ha modificat el codi
% respecte la tasca 3.
M_9=sqrt(2/(gamma_t-1)*((pi_primari_total)^((gamma_t-1)/gamma_t)-1));

if M_9<1 % Hipòtesi certa
    P_9=P_0;
    pi_primari=1;
else    %Hipòtesi falsa: tovera ofegada, el Mach no pot ser superior a la unitat.
    M_9=1;
    pi_primari=pi_primari_total*(2/(gamma_t+1))^(gamma_t/(gamma_t-1));
    P_9=pi_primari*P_0;
    % On: pi_sec=P_19/P_0. S'utilitza gamma_c ja que és fluxe primari
end
T_9=T_t9/(1+0.5*(gamma_t-1)*M_9^2);
u_9=M_9*sqrt(gamma_t*R_g*T_9);

%% (0)-(1.3): SECONDARY INLET

%A la presa secondària d'aire, la temperatura d'estancament és major que la
%del flux lliure.
T_t13=T_0*theta_0*tau_f;

%% (1.3)-(1.9): SECONDARY NOZZLE

%Com l'aire no experimenta cap altre procés després de la presa secondària,
%la temperatura d'estancament a la tovera secondària és la mateixa que a la
%presa. Tot i això cal diferenciar si aquesta es troba xocada o no.
T_t19=T_t13;

% pi_secundari_total=P_t19/P_0
pi_secundari_total=delta_0*pi_d*pi_f*pi_ns;

%A continuació, s'imposa tovera adaptada. Prèviament, analitzant els
%resultats de la tasca 3, s'ha vist que a diferència de la tovera primària,
%la tovera secundària estava ofegada.
M_19_mat=sqrt(2/(gamma_c-1)*(pi_secundari_total^((gamma_c-1)/gamma_c)-1));
P_19=P_0;
pi_secundari=1;

T_19=T_t19/(1+(gamma_c-1)/2*M_19_mat^2);
u_19=M_19_mat*sqrt(gamma_c*R_g*T_19);

%% Càlcul del l'empenyiment específic

%Es plantegen les equacions d'empenyiment dels dos tipus de flux i es fan
%relatives al flux d'aire que veu el motor. L'empenta total que veu el
%motor és la suma de les dues.
Specific_Thrust_primari=(1+f)*u_9-u_0+(1+f)*T_9/u_9*R_g*(1-1/pi_primari);
Specific_Thrust_secundari=alpha*(u_19-u_0)+alpha*(T_19*R_g)/u_19*(1-1/pi_secundari);
Specific_Thrust=Specific_Thrust_primari+Specific_Thrust_secundari;

%% Càlcul de l'impuls específic

%L'impuls específic es calcula com la relació entre la força produïda i el
%combustible emprat.
I_sp=Specific_Thrust/(f*g)

% Guany relatiu d'impuls específic
I_sp_ap3=5.655999018289032e+03;
delta_I_sp_ap4=I_sp-I_sp_ap3
G_I_sp_ap4=I_sp/I_sp_ap3

%% Càlcul de la relació d'àrees
M_18_mat=1;
% Respectius mass flow parameters
MFP_18=sqrt(gamma_c)*(M_18_mat/(1+0.5*(gamma_c-1)*M_18_mat^2)^((gamma_c+1)/(2*(gamma_c-1))));
MFP_19=sqrt(gamma_c)*(M_19_mat/(1+0.5*(gamma_c-1)*M_19_mat^2)^((gamma_c+1)/(2*(gamma_c-1))));
%Còmput relació d'àrees tovera secundària
relacio_area=MFP_18/MFP_19;


