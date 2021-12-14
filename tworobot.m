clear all
clc

T=0.001;N=10000;K1=eye(2,2);K2=eye(2,2);t=0:T:T*(N-1);
m1=1;m2=1.5;l1=0.2;l2=0.3;lc1=0.1;lc2=0.15;J1=0.013;J2=0.045;g=9.8;
I=ones(1,N);
qd=[0.3*I; 0.1*I];dqd=zeros(2,N);ddqd=zeros(2,N);%参考(1)
% qd=[(0.3+0.02*sin(t)); (0.1+0.01*cos(t))];dqd=[0.02*cos(t);-0.01*sin(t)];ddqd=[-0.02*sin(t);-0.01*cos(t)];%参考(2)
q=zeros(2,N);dq=zeros(2,N);ddq=zeros(2,N);%变量

theta=[m1*lc1*lc1+m2*(l1*l1+lc2*lc2)+J1+J2 m2*l1*lc2 m2*lc2*lc2+J2 m1*lc1+m2*l1 m2*lc2];
for i=1:N-1
	q(:,i+1)=dq(:,i)*T+q(:,i);
	dq(:,i+1)=ddq(:,i)*T+dq(:,i);
	M=[theta(1)+2*theta(2)*cos(q(2,i+1)) theta(3)+theta(2)*cos(q(2,i+1));theta(3)+theta(2)*cos(q(2,i+1)) theta(3)];
	C=[-theta(2)*sin(q(2,i+1))*dq(2,i+1) -theta(2)*sin(q(2,i+1))*(dq(2,i+1)+dq(1,i+1));theta(2)*sin(q(2,i+1))*dq(1,i+1) 0];
	ddq(:,i+1)=M\(-K1*(q(:,i+1)-qd(:,i+1))-(C+K2)*(dq(:,i+1)-dqd(:,i+1))+M*ddqd(:,i+1));
end
e=q-qd;de=dq-dqd;dde=ddq-ddqd;
subplot(2,2,1);plot(t,q(1,:));hold on;plot(t,qd(1,:));%grid on;
subplot(2,2,3);plot(t,e(1,:),'-');%grid on;
subplot(2,2,2);plot(t,q(2,:));hold on;plot(t,qd(2,:));%grid on;
subplot(2,2,4);plot(t,e(2,:),'-');%grid on;
figure;
subplot(2,2,1);plot(t,dq(1,:));hold on;plot(t,dqd(1,:));%grid on;
subplot(2,2,3);plot(t,de(1,:),'-');%grid on;
subplot(2,2,2);plot(t,dq(2,:));hold on;plot(t,dqd(2,:));%grid on;
subplot(2,2,4);plot(t,de(2,:),'-');%grid on;

