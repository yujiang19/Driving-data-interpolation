function deltdis = optv(sv, ds, x1, x2, a0, a2, vc, vcni)

n1 = [0; vc];
n2 = [10; vcni];
nc1 = [x1; vc+x1*a0];
nc2 = [x1+x2; vcni-a2*(10-x1-x2)];
A = [n1(1)^3 n1(1)^2 n1(1)^1 1; n2(1)^3 n2(1)^2 n2(1)^1 1;...;
    nc1(1)^3 nc1(1)^2 nc1(1)^1 1; nc2(1)^3 nc2(1)^2 nc2(1)^1 1];
Y =  [n1(2); n2(2); nc1(2); nc2(2)];
cofs = inv(A)*Y;
ac = @(xx)cofs(1)*xx.^3 + cofs(2)*xx.^2 + cofs(3)*xx.^1 + cofs(4);
x1 = 0:1:10;
y1 = ac(x1);
s1 = trapz(x1,y1);
deltdis = abs(s1 - sv);

end