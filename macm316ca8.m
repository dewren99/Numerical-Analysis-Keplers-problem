h = 0.0005;
T=200;
N = T/h; %=400000
q=ones(N+1,2); %3D, 400001 by 2 by 1
size(q);
p=ones(N+1,2); %400001 by 2 by 1
A=ones(N+1,1); %Angular momentum A(t)
H=ones(N+1,1); %Hamiltonian H(t)

A_symplecticEuler = ones(N,1);
H_symplecticEuler = ones(N,1);

e = 0.6;

q(1,1) = 1 - e; %=0.4
q(1,2) = 0;
p(1,1) = 0;
p(1,2) = sqrt((1+e)/(1-e)); %=2

for i = 1:1:N+1
    p_der = (-1)/(((q(i,1)^2 +q(i,2)^2))^(3/2));
    
    p_one = p_der*q(i,1);
    p_two = p_der*q(i,2);
    
    p(i+1,1)=p(i,1)+h*p_one; %Euler's method
    p(i+1,2)=p(i,2)+h*p_two; %Euler's method
    
    q_one = p(i,1);
    q_two = p(i,2);

    q(i+1,1)=q(i,1)+h*q_one; %Euler's method
    q(i+1,2)=q(i,2)+h*q_two; %Euler's method
    
    A(i)=q(i,1)*p(i,2)-q(i,2)*p(i,1); %Angular Momentum A(t)
    H(i)=(1/2)*(p(i,1)^2+p(i,2)^2)-1/(sqrt(q(i,1)^2+q(i,2)^2)); %Hamiltonian H(t)
end

figure(1)
plot(q(:,1),q(:,2),'r.');
title('Fig 1. Part1 q1-q2 plane');

figure(2)
plot(A,'r');
title('Angular Momentum');
ylabel('A(t)')

figure(3)
plot(H,'r');
title('Hamiltonian')
ylabel('H(t)')

for i = 1:1:N
    q(i+1,1) = q(i,1)+h*p(i,1);
    q(i+1,2) = q(i,2)+h*p(i,2);
    
    p_euler=(h)/(((q(i+1,1)^2 +q(i+1,2)^2))^(3/2));
    
    p(i+1,1)=p(i,1)-p_euler*q(i+1,1);
    p(i+1,2)=p(i,2)-p_euler*q(i+1,2);
    
   A_symplecticEuler(i)=q(i,1)*p(i,2)-q(i,2)*p(i,1);
   H_symplecticEuler(i)=(1/2)*(p(i,1)^2+p(i,2)^2)-1/(sqrt(q(i,1)^2+q(i,2)^2));
end

figure(4)
plot(q(:,1),q(:,2),'r.');

figure(5)
plot(A_symplecticEuler,'r');

figure(6)
plot(H_symplecticEuler,'r');
