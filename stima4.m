function W=stima4(vertices)
mk=1/2*det([ones(1,3);vertices']);
%qin=qin(nodes);
W=[2 1 1;1 2 1;1 1 2]*mk*2/24;
%%%%%%%% assembled matrix4 %%%%%%%%%%%%%
%%%  for j=1:size(Elem,1)
%       V(Elem(j,:),Elem(j,:))=V(Elem(j,:),Elem(j,:))+stima4(Coord(Elem(j,:),:));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5          
% qin is the input vector
