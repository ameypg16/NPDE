function U=stima2(vertices,qin,nodes)
mk=1/2*det([ones(1,3);vertices']);
qin=qin(nodes);
S=zeros(6,1);
S(1,1)=qin(1,1)^2; S(4,1)=qin(1,1)*qin(1,2)*2;
S(2,1)=qin(1,2)^2; S(5,1)=qin(1,2)*qin(1,3)*2;
S(3,1)=qin(1,3)^2; S(6,1)=qin(1,1)*qin(1,3)*2;
T11=[24 4 4 2*6 2*2 2*6]*2*mk/720;
T12=[6 6 2 2*2*2 2*2 2*2]*2*mk/720;
T13=[6 2 6 2*2 2*2 2*4]*2*mk/720;
T22=[4 24 4 2*6 2*6 2*2]*2*mk/720;
T23=[2 6 6 2*2 2*4 2*2]*2*mk/720;
T33=[4 4 24 2*2 2*6 2*6]*2*mk/720;
T21=T12;T32=T23;T31=T13;
U=zeros(3,3);
U(1,1)=T11*S;U(1,2)=T12*S;
U(1,3)=T13*S;U(2,2)=T22*S;
U(2,3)=T23*S;U(3,3)=T33*S;
U(2,1)=U(1,2);U(3,2)=U(2,3);U(3,1)=U(1,3);

%%%%%%%% assembled matrix2 %%%%%%%%%%%%%
%%%  for j=1:size(Elem,1)
%       U(Elem(j,:),Elem(j,:))=U(Elem(j,:),Elem(j,:))+3*stima2(Coord(Elem(j,:),:),q1in,Elem(j,:));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5          
% q1in is the frist component of the input vector




