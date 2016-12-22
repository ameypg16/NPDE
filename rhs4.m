function Z2=rhs4(vertices,qin,nodes)
mk=1/2*det([ones(1,3);vertices']);
qin=qin(nodes);
L1=[ones(1,3);vertices']'\[1;0;0];
L2=[ones(1,3);vertices']'\[0;1;0];
L3=[ones(1,3);vertices']'\[0;0;1];

% element stiffness matrix
M=mk*[L1(2)*L1(2)+L1(3)*L1(3) L1(2)*L2(2)+L1(3)*L2(3) L1(2)*L3(2)+L1(3)*L3(3)
      L2(2)*L1(2)+L2(3)*L1(3) L2(2)*L2(2)+L2(3)*L2(3) L2(2)*L3(2)+L2(3)*L3(3)
      L3(2)*L1(2)+L3(3)*L1(3) L3(2)*L2(2)+L3(3)*L2(3) L3(2)*L3(2)+L3(3)*L3(3)];
a=zeros(3,1);
for i=1:3
    qin(1,i)=a(i,1);
end
Z2=M*a;