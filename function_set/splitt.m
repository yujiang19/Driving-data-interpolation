function se_index = splitt(m)
t_length = length(m.M.t);
Vms = m.M.v;
i=1;
n=1;
while(i<=t_length-1)
    if Vms(i)~=0&&Vms(i+1)==0
        stop1_index(n)=i+1;
        n=n+1;
    end
    i=i+1;
end
i=1;
n=2;
sp_l=length(stop1_index);
while(i<sp_l)
    j=stop1_index(i);
    start1_index(n)=j+1;
    n=n+1;
    i=i+1;
end
start1_index(1,1) = 1;
diff_d = stop1_index - start1_index;
v_mean = [];
for i = 1:length(stop1_index)
    v_mean = [v_mean mean(Vms(start1_index(i):stop1_index(i)))];
end
se_index = [start1_index; stop1_index];
f_index = find(v_mean<2);
se_index(:,f_index) = [];
seg_length = se_index(2,:) - se_index(1,:);
le_index = find(seg_length<30);
se_index(:,le_index) = [];
end