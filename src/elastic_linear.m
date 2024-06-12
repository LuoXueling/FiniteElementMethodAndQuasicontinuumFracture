function [sigma,strainEnergy,stiffnessDensity] = elastic_linear(epsilon,s,materialCfg)
K=materialCfg.K;
bulk=materialCfg.bulk;
shear=materialCfg.shear;
lambda=materialCfg.lambda;
E=materialCfg.E;
v=materialCfg.v;
len=length(epsilon);

% CALCULATE Cep FIRST!!!!
Cep=s*K;
strainEnergy(1)=0.5*(K*epsilon)'*epsilon; % positive part
strainEnergy(2)=0; %negative part
sigma=Cep*epsilon; 
stiffnessDensity=Cep;
end
