function modelMatrix = TwoCompartmentForwardModel(Ki, lamda, Vb, HCT, t)
% compute the model matrix for two compartment model ignoring extraction ratio E

k2 = Ki * (1-HCT)/lamda;
modelMatrix = Vb + Ki*exp(-k2*t);