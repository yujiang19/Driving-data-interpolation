function deltdis = optv3(sv, ds, x1, x2, a0, vc, vcni, vrr, vcnir, aa_max, k, a_cur, a2)
W1 = x1 - (1.01.^(1:1:k)-1);
W2 = x2 - (1.01.^(1:1:k)-1);
a1 = -a0.*W1;
for j = 1:floor(k/ds)
    vv = vc + a1(j)*ds;
    vrr = [vrr vv];
    vc = vrr(end);
end
if -a0<0
    a2 = min([-a0,-a_cur,-a2]).*W2;  
else
    a2 = aa_max.*W2;
end
for j = 1:floor(k/ds)
    vv = vcni + a2(j)*ds;
    vcnir = [vcnir vv];
    vcni = vcnir(end);
end
vcnir = flip(vcnir);
[m indexv] = min(abs(vrr-vcnir));
vr = [vrr(1:indexv) vcnir(indexv+1:end)];
vr(vr<0)=0;
t1 = 0:1:10;
s1 = trapz(t1,vr);
deltdis = abs(s1 - sv);
end