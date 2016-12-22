function X=stima5(vertices,qin,pin,nodes)
mk=1/2*det([ones(1,3);vertices']);
qin=qin(nodes);
pin=pin(nodes);
T11=[24 4 4 2*6 2*2 2*6]*2*mk/720;
T12=[6 6 2 2*2*2 2*2 2*2]*2*mk/720;
T13=[6 2 6 2*2 2*2 2*4]*2*mk/720;
T22=[4 24 4 2*6 2*6 2*2]*2*mk/720;
T23=[2 6 6 2*2 2*4 2*2]*2*mk/720;
T33=[4 4 24 2*2 2*6 2*6]*2*mk/720;
T21=T12;T32=T23;T31=T13;
S=zeros(6,1);
S(1,1)=qin(1,1)*pin(1,1); S(4,1)=qin(1,1)*pin(1,2)+qin(1,2)*pin(1,1);
S(2,1)=qin(1,2)*pin(1,2); S(5,1)=qin(1,2)*pin(1,3)+qin(1,3)*pin(1,2);
S(3,1)=qin(1,3)*pin(1,3); S(6,1)=qin(1,1)*pin(1,3)+qin(1,3)*pin(1,1);
X=zeros(3,3);
X(1,1)=T11*S;X(1,2)=T12*S;
X(1,3)=T13*S;X(2,2)=T22*S;
X(2,3)=T23*S;X(3,3)=T33*S;
X(2,1)=X(1,2);X(3,2)=X(2,3);X(3,1)=X(1,3);
%%%%%%%% assembled matrix2 %%%%%%%%%%%%%
%%%  for j=1:size(Elem,1)
%       U(Elem(j,:),Elem(j,:))=U(Elem(j,:),Elem(j,:))+2*stima2(Coord(Elem(j,:),:),q1in,q2in,Elem(j,:));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5          
% qin=[q1in,q2in] is the input vector


