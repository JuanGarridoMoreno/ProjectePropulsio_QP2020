function [T_0,P_0]=atmosferaISA(z,T_SL,P_SL)
T_0=T_SL-6.5*10^(-3)*z;
P_0=P_SL*(1-(6.5*10^(-3)*z)/(T_SL))^5.252;
end
