function plotTimming_GT_UT(timing1, timing2)
% plot GT UT timing

n = size(timing1, 1);

t1 = zeros(n, 1);
t2 = zeros(n, 1);

for ii=1:n
    t1(ii) = timing1{ii, 1};
    t2(ii) = timing2{ii, 1};
end

disp( ['total time: ' num2str(sum(t1)) ' ' num2str(sum(t2))] );

figure;
hold on
plot([1:n], t1, '*b');
plot([1:n], t2, 'or');
legend('baseline', 'test');
hold off
