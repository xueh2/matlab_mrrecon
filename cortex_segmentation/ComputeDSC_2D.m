
function current_corrR = ComputeDSC_2D(slice, reflected_sliceRotated_Transi)
% compute the DSC value

num1 = length(find(slice>0));
num2 = length(find(reflected_sliceRotated_Transi>0));

pp = slice & reflected_sliceRotated_Transi;

corrNum = length(find(pp>0));

current_corrR = 2*corrNum / (num1+num2);