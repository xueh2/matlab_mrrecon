num = 400;

p = randn(num,1);
mean(p)
std(p)
p = p + 0.8;
p2 = randn(num,1);
p2 = p2 - 0.8;
data = [p; p2];
label = ones(2*num, 1);
label(num+1:end) = 2;

foldnum = 5
C_low = -20;
C_high = -6;
C_step = 2;
gamma_low = -35;
gamma_high = -20;
gamma_step = 3;
fine_range = 2;
fine_step = 0.25;
finesearch_flag = 1;

[goodC, goodGamma] = GridSearch_C_gamma(label, data, foldnum,...
    C_low, C_high, C_step, ...
    gamma_low, gamma_high, gamma_step, ...
    finesearch_flag, fine_range, fine_step);

libsvm_options = ['-c ' num2str(2^goodC) ' -g ' num2str(2^goodGamma)]

model = svmtrain_Version2p82(label, data, libsvm_options);
[predicted_label, accuracy, prob_estimates] = svmpredict_Version2p82(label, data, model);