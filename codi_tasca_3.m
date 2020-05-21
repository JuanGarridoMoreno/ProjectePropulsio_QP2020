% Aquest codi ha estat realitzat per:
% -Alexis Leon Delgado
% -David Morante Torra
% -Ferran Rubio Vallhonrat
% -Juan Garrido Moreno

% En base al set preliminar de valors de la tasca 1, aquest programa troba
% el set que fa màxim l'impuls específic d'acord a l'interval especificat
% per a cada variable a la tasca 2

clear all
close all
clc

format long

% Introducció de les dades del motor a analitzar. Totes elles expressades
% segons el sistema internacional d'unitats.
M_0=0.85; z=11000; T_SL=288.15; P_SL=101325;
pi_d=0.98; eta_f=0.89; eta_lpc=0.88; eta_hpc=0.86;
pi_b=0.99; eta_b=0.99; eta_hpt=0.91; eta_lpt=0.92;
eta_mH=0.993; eta_mL=0.997; pi_np=0.99; pi_ns=0.99;
T_t4=1450; gamma_c=1.4; gamma_t=1.3; R_g=287; h=43000000; g=9.81;

% Inicialització dels valors de comparació
I_sp_com=0; 
Specific_Thrust_com=0;

% Incialització dels valors obtinguts a la tasca 1
Specific_Thrust_ap1=1.317407200719343e+03; 
I_sp_ap1=5.213737456185124e+03;

%% CÀLCULS PREVIS

% Es calculen els coeficients de pressió per l'aire del compressor (Cp_c) i
% de la turbina (Cp_t), segons les gammes depenent de si ha passat o no per
% la cambra de combustió.
Cp_t=(gamma_t/(gamma_t-1))*R_g;
Cp_c=(gamma_c/(gamma_c-1))*R_g;

% S'obtenen també els valors de pressió i temperatura del flux lliure a
% l'altura de funcionament. Amb aquests es calula la velocitat d'entrada de
% l'aire al motor.
[T_0,P_0]=atmosferaISA(z,T_SL,P_SL);
u_0=M_0*sqrt(gamma_c*R_g*T_0);
%% DEFINICIÓ DE RANGS D'ESTUDI

% El grup de paràmetres escollits a la tasca 1 s'estableix com a nominal
alpha_nominal=8;
pi_f_nominal=1.6;
pi_lpc_nominal=4.3;
pi_hpc_nominal=7.5;
valors_nominals=[alpha_nominal,pi_f_nominal,pi_lpc_nominal,pi_hpc_nominal];

% Es defineix la sensibilitat (variació +- del paràmetre)
sensibilitat_alpha=2;
sensibilitat_pi_f=0.5;
sensibilitat_pi_lpc=1;
sensibilitat_pi_hpc=1;

% Nombre d'elements per paràmetre variat
N=20;

% Definició de límits per a la relació de derivació
alpha_inf=alpha_nominal-sensibilitat_alpha;
alpha_sup=alpha_nominal+sensibilitat_alpha;
% Definició del vector de valors de la relació de derivació
vec_alpha=linspace(alpha_inf,alpha_sup,N);

% Definició de límits per a la relació de pressions del fan
pi_f_inf=pi_f_nominal-sensibilitat_pi_f;
pi_f_sup=pi_f_nominal+sensibilitat_pi_f;
% Definició del vector de valors de la relació de pressions del fan
vec_pi_f=linspace(pi_f_inf,pi_f_sup,N);

% Definició de límits per a la relació de pressions del LPC
pi_lpc_inf=pi_lpc_nominal-sensibilitat_pi_lpc;
pi_lpc_sup=pi_lpc_nominal+sensibilitat_pi_lpc;
% Definició del vector de valors de la relació de pressions del LPC
vec_pi_lpc=linspace(pi_lpc_inf,pi_lpc_sup,N);

% Definició de límits per a la relació de pressions del HPC
pi_hpc_inf=pi_hpc_nominal-sensibilitat_pi_hpc;
pi_hpc_sup=pi_hpc_nominal+sensibilitat_pi_hpc;
% Definició del vector de valors de la relació de pressions del HPC
vec_pi_hpc=linspace(pi_hpc_inf,pi_hpc_sup,N);

% Generació d'una matriu buida de combinacions possibles
valors_possibles=[];

%% ESTUDI DE TOTES LES COMBINACIONS

% Es crea un bucle que analitza totes les combinacions possibles
% Com es tenen 4 variables, es tracta de 4 bucles for anidats
for i=1:length(vec_alpha)
    for j=1:length(vec_pi_f)
        for k=1:length(vec_pi_lpc)
            for l=1:length(vec_pi_hpc)
              alpha=vec_alpha(i);
              pi_f=vec_pi_f(j);
              pi_lpc=vec_pi_lpc(k);
              pi_hpc=vec_pi_hpc(l);

        % Resolució del motor d'acord a la metodologia analítica de la
        % tasca 1 per a cada valor del paràmetre variable
                      
% (0)-(2): INLET

theta_0=1+0.5*(gamma_c-1)*M_0^2; % theta_0=T_t0/T_0
delta_0=theta_0^(gamma_c/(gamma_c-1)); % delta_0=P_t0/P_0

T_t0=T_0*theta_0;
T_t2=T_t0;

P_t0=P_0*delta_0;
P_t2=P_t0*pi_d;

% (2)-(2.5): LPC
P_t25=P_t2*pi_lpc;
tau_lpc=1+(pi_lpc^((gamma_c-1)/gamma_c)-1)/eta_lpc;
T_t25=T_t2*tau_lpc;

% (2.5)-(3): HPC
P_t3=P_t25*pi_hpc;
tau_hpc=1+(pi_hpc^((gamma_c-1)/gamma_c)-1)/eta_hpc;
T_t3=T_t25*tau_hpc;

% (3)-(4): COMB
f=(Cp_t*T_t4-Cp_c*T_t3)/(eta_b*h-Cp_t*T_t4);
P_t4=P_t3*pi_b;

% (4)-(4.5): HPT
tau_hpt=1-(1/eta_mH)*(Cp_c/Cp_t)*(1/(1+f))*(T_t25/T_t4)*(tau_hpc-1);
pi_hpt=(1-(1-tau_hpt)/eta_hpt)^(gamma_t/(gamma_t-1));
T_t45=T_t4*tau_hpt;
P_t45=P_t4*pi_hpt;

% (4.5)-(5): LPT
tau_f=1+(pi_f^((gamma_c-1)/gamma_c)-1)/eta_f;
tau_lpt=1-(1/eta_mL)*(Cp_c/Cp_t)*(1/(1+f))*(T_t0/T_t4)*(1/tau_hpt)*((tau_lpc-1)+alpha*(tau_f-1));
pi_lpt=(1-(1-tau_lpt)/eta_lpt)^(gamma_t/(gamma_t-1));
T_t5=T_t45*tau_lpt;
P_t5=P_t45*pi_lpt;

% (5)-(9): PRIMARY NOZZLE
P_t9=P_t5*pi_np;
T_t9=T_t5;
pi_primari_total=P_t9/P_0; %Poner esto o cadena?
M_9=sqrt(2/(gamma_t-1)*((pi_primari_total)^((gamma_t-1)/gamma_t)-1));
if M_9<1
    P_9=P_0;
    pi_primari=1;
else
    M_9=1;
    pi_primari=pi_primari_total*(2/(gamma_t+1))^(gamma_t/(gamma_t-1)); % On: pi_sec=P_19/P_0. S'utilitza gamma_c ja que és fluxe secundari
    P_9=pi_primari*P_0;
end
T_9=T_t9/(1+0.5*(gamma_t-1)*M_9^2);
u_9=M_9*sqrt(gamma_t*R_g*T_9);

% (0)-(1.3): SECONDARY INLET
T_t13=T_0*theta_0*tau_f;

% (1.3)-(1.9): SECONDARY NOZZLE
T_t19=T_t13;

pi_secundari_total=delta_0*pi_d*pi_f*pi_ns; % pi_secundari_total=P_t19/P_0
M_19=sqrt(2/(gamma_c-1)*((pi_secundari_total)^((gamma_c-1)/gamma_c)-1));
if M_19<1
    P_19=P_0;
    pi_secundari=1;
else
    M_19=1;
    pi_secundari=pi_secundari_total*(2/(gamma_c+1))^(gamma_c/(gamma_c-1));% On: pi_sec=P_19/P_0. S'utilitza gamma_c ja que és fluxe secundari
    P_19=pi_secundari*P_0;
end
T_19=T_t19/(1+(gamma_c-1)/2*M_19^2);
u_19=M_19*sqrt(gamma_c*R_g*T_19);

% Càlcul del Thrust específic
Specific_Thrust_primari=(1+f)*u_9-u_0+(1+f)*T_9/u_9*R_g*(1-1/pi_primari);
Specific_Thrust_secundari=alpha*(u_19-u_0)+alpha*(T_19*R_g)/u_19*(1-1/pi_secundari);
Specific_Thrust=Specific_Thrust_primari+Specific_Thrust_secundari;

% Specific impuls
I_sp=Specific_Thrust/(f*g);

% Degut a la possibilitat d'obtenir valors complexos en el resultat
% a causa de la impossibilitat de complir el balanç de potències,
% els resultats es consideran nuls en aquesta situació (són físicament impossibles)
if imag(I_sp)~=0||imag(Specific_Thrust)~=0
Specific_Thrust=0;
I_sp=0;
end

% Aquest condicional emmagatzema la combinació que presenta màxim impuls
% específic
if I_sp>I_sp_com
    I_sp_com=I_sp;
    alpha_I_sp=alpha;
    pi_f_I_sp=pi_f;
    pi_lpc_I_sp=pi_lpc;
    pi_hpc_I_sp=pi_hpc;
    P_9_I_sp=P_9;
    P_19_I_sp=P_19;
    M_9_I_sp=M_9;
    M_19_I_sp=M_19;
    Specific_Thrust_I_sp=Specific_Thrust;
end

% Aquest condicional emmagatzema totes les combinacions que presenten un
% impuls específic superior al set de paràmetres prelimnar (tasca 1)
if I_sp>I_sp_ap1
    mataux=[alpha pi_lpc pi_hpc pi_f I_sp Specific_Thrust].';
    valors_possibles=[valors_possibles mataux];
end


            end
        end
    end
end

% Impressió de resultats

% Guany relatiu d'impuls específic
delta_I_sp_ap3=I_sp_com-I_sp_ap1;
G_I_sp_ap3=((I_sp_com/I_sp_ap1)-1)*100

% Guany relatiu d'empenta específica
delta_Specific_Thrust_ap3=Specific_Thrust_I_sp-Specific_Thrust_ap1;
G_Specific_Thrust_ap3=((Specific_Thrust_I_sp/Specific_Thrust_ap1)-1)*100

% Overall pressure ratio
O_PR=pi_hpc_I_sp*pi_lpc_I_sp