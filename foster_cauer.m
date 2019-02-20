clear all;
RC_Foster=importdata('RC_foster.txt')
m=length(RC_Foster);
for i=[1:m]
R_Foster(i)=RC_Foster(i,1);
C_Foster(i)=RC_Foster(i,2);
end
Ztf=0;

%sを伝達関数で定義
    for i=[1:m]
        s=tf('s');
        Ztf=tf(Ztf+(R_Foster(i)/(1+R_Foster(i)*C_Foster(i)*s)));
        [N_tf,D_tf] = tfdata(Ztf);
        chain_length = numel(D_tf{1})-1;
        D=N_tf{1}; D(1)=[]; N=D_tf{1};
    end
[N_v,D_v]=tfdata(Ztf,'v');%simulinkパラメータ用
    
%sをシンボルとして定義
%     for i=[1:m]
%         syms s
%         Ztf=Ztf+(R_Foster(i)/(1+R_Foster(i)*C_Foster(i)*s))
%         [N_tf,D_tf] = numden(Ztf);
%         D = sym2poly(N_tf);
%         N = sym2poly(D_tf);
%         chain_length = numel(N)-1;
%     end

%連分数展開
for i = 1:chain_length
    n = max([N D]);
    N = N/n; D = D/n;
    C_Cauer(i)=N(1)/D(1);
    Dnew = D*D(1);
    Nnew = N*D(1)-[D 0]*N(1); Nnew(1) = [];
    R_Cauer(i) = Dnew(1)/Nnew(1); 
    D = Dnew*Nnew(1)-Nnew*Dnew(1); D(1) = [];
    N = Nnew*Nnew(1);
end