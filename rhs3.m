function Z1=rhs3(vertices,qin,nodes)
mk=1/2*det([ones(1,3);vertices']);
qin=qin(nodes);
W=[2 1 1;1 2 1;1 1 2]*2*mk/24;
a=zeros(3,1);
for i=1:3
    a(i,1)=qin(1,i);
end
Z1=W*a;