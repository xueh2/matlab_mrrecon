function T2 = T1_to_T2(T1_0, T1, T2_0, r1, r2)
% compute T2 from T1 Gd
% all times are in ms
% T2 = T1_to_T2(T1_0, T1, T2_0, r1, r2)

gd = T1_to_Gd(T1_0, T1, r1);
[T1_2, T2] = Gd2Time(T1_0, T2_0, r1, r2, gd);