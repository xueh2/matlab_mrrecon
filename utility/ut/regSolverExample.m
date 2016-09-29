
[A,b_bar,x] = shaw(32);
randn ('seed',41997);
e = 1e-3*randn(size(b_bar)); 
b = b_bar + e;
[U,s,V] = csvd (A);
figure
hold on
subplot (2,1,1); picard (U,s,b_bar);
subplot (2,2,2); picard (U,s,b);
hold off

lambda = [1,3e-1,1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5];
X_tikh = tikhonov(U,s,V,b,lambda);
F_tikh = fil_fac(s,lambda);
iter = 30; reorth = 0;
[X_lsqr,rho,eta,F_lsqr] = lsqr_b(A,b,iter,reorth,s);
subplot (2,2,1); surf (X_tikh), axis ('ij'), title ('Tikhonov solutions')
subplot (2,2,2); surf (log10 (F_tikh)), axis ('ij'), title ('Tikh ¯lter factors, log scale')
subplot (2,2,3); surf (X_lsqr (:,1:17)), axis ('ij'), title ('LSQR solutions')
subplot (2,2,4); surf (log10(F_lsqr)), axis ('ij'), title ('LSQR ¯lter factors, log scale')

subplot (1,2,1); l_curve (U,s,b); axis ([1e-3,1,1,1e3])
subplot (1,2,2); plot_lc(rho,eta,'o'), axis([1e-3,1,1,1e3])

n = 16;
[A,b,x] = i_laplace (n,2);
b = b + 1e-4*randn(size(b));
L = get_l (n,1);
[U,s,V] = csvd (A); [UU,sm,XX] = cgsvd (A,L);
I = 1; 
for i=[3,6,9,12]
subplot (2,2,I); plot (1:n,V(:,i)); axis ([1,n,-1,1])
xlabel (['i = ',num2str(i)]), I = I + 1;
end

subplot (2,2,1), text (12,1.2,'Right singular vectors V(:,i)')

I = 1; 
for i=[n-2,n-5,n-8,n-11]
subplot (2,2,I); plot (1:n,XX(:,i)); axis ([1,n,-1,1])
xlabel (['i = ',num2str(i)]), I = I + 1;
end
subplot (2,2,1), text (10,1.2,'Right generalized singular vectors XX(:,i)')

k_tsvd = 7; k_tgsvd = 6;
X_I = tsvd (U,s,V,b,1:k_tsvd);
X_L = tgsvd (UU,sm,XX,b,1:k_tgsvd);
X_D = dsvd (UU,sm,XX,b,lambda);
X_Tik = tikhonov (UU,sm,XX,b,lambda);
figure; plot (1:n,X_I,1:n,x,'x'), axis ([1,n,0,1.2]), xlabel ('L = I')
figure; plot (1:n,X_L,1:n,x,'x'), axis ([1,n,0,1.2]), xlabel ('L \neq I')

figure; plot (1:n,X_Tik,1:n,x,'x'), axis ([1,n,0,1.2]), xlabel ('Tikhonov, L \neq I')
figure; l_curve (UU,sm,b, 'Tikh', L);

figure; plot (1:n,X_D,1:n,x,'x'), axis ([1,n,0,1.2]), xlabel ('Damped svd, L \neq I')
figure; l_curve (UU,sm,b, 'Tikh', L);

picard(U,s,b);





