function dy=equa(t,y)

%global sigma delta lambda beta good%m omega gamma GA good k ktps
global k mu dphim delta gamma do omega dphio vecd indo

dy = zeros(9,1);    % a column vector

dphi=dphio+dphim*cos(omega*t);
indo=indo+1;
%vecd(indo)=dphi;

%order : Real(Ex) Im(Ex) Real(Ey) Im(Ey) Real(Px) Im(Px) Real(Py) Im(Py) D 

dy(1) = -k*y(1)+mu*y(5);  %Ex real
dy(2) = -k*y(2)+mu*y(6);  %Ex ima
dy(3) = -k*y(3)+mu*y(7)-dphi*y(4);  %Ey real
vecd(indo,1)=-k*y(3)+mu*y(7);
vecd(indo,2)=dphi*y(4);

dy(4) = -k*y(4)+mu*y(8)+dphi*y(3);  %Ey ima
vecd(indo,4)=-k*y(4)+mu*y(8);
vecd(indo,5)=dphi*y(3);
vecd(indo,3)=t;

dy(5) = -(y(5)-delta*y(6)-y(1)*y(9));  %Px real
dy(6) = -(y(6)+delta*y(5)-y(2)*y(9));  %Px ima
dy(7) = -(y(7)-delta*y(8)-y(3)*y(9));  %Py real
dy(8) = -(y(8)+delta*y(7)-y(4)*y(9));  %Py ima
dy(9) = -gamma*(y(9)-do+(y(1)*y(5)+y(2)*y(6)+y(3)*y(7)+y(4)*y(8))); %D inversion de population
