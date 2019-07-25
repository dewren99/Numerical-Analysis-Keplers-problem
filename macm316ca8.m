h = 0.0005;
T=200;
N = T/h; %=400000
q=ones(N+1,2);
size(q);
p=ones(N+1,2);
A=ones(N+1,1); %Angular momentum A(t)
H=ones(N+1,1); %Hamiltonian H(t)

A_symplecticEuler = ones(N+1,1);
H_symplecticEuler = ones(N+1,1);

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
plot(q(:,1),q(:,2));
title("Fig 1. Position of the moving planet at time t_n using Euler's method")
xlabel('q1')
ylabel('q2')

figure(2)
plot(0:N,A);
title('Fig 2a. Angular Momentum');
xlabel('t')
ylabel('A(t)')

figure(3)
plot(0:N,H);
title('Fig 2b. Hamiltonian')
xlabel('t')
ylabel('H(t)')

for i = 1:1:N+1
    q(i+1,1) = q(i,1)+h*p(i,1);
    q(i+1,2) = q(i,2)+h*p(i,2);
    
    p_euler=(h)/(((q(i+1,1)^2 +q(i+1,2)^2))^(3/2));
    
    p(i+1,1)=p(i,1)-p_euler*q(i+1,1);
    p(i+1,2)=p(i,2)-p_euler*q(i+1,2);
    
   A_symplecticEuler(i)=q(i,1)*p(i,2)-q(i,2)*p(i,1);
   H_symplecticEuler(i)=(1/2)*(p(i,1)^2+p(i,2)^2)-1/(sqrt(q(i,1)^2+q(i,2)^2));
end

figure(4)
plot(q(:,1),q(:,2),'r');
title("Fig 3a. Position of the moving planet at time t_n using symplectic Euler method");
xlabel('q1');
ylabel('q2');

figure(5)
plot(0:N,A_symplecticEuler,'r');
title('Fig 3b. Angular Momentum');
xlabel('t');
ylabel('A(t)');


figure(6)
plot(0:N,H_symplecticEuler,'r');
title('Fig 3c. Hamiltonian');
xlabel('t');
ylabel('H(t)');
