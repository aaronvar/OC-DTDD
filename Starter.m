% close all;
clc;

nnodes=100; %number of nodes
carr_freq=1.9*10^9; %carrier frequency
beta=3.6; %path loss exponent
sidelength=2000; %rho - dimension of area the points are placed (m)
SIsuppression=110;%SIsuppression_var; %SI suppression level in dB
node_power=ones(nnodes,1)*ind_var;%20;%transmit power of each node
%% initiate model
dist_matrix=zeros(nnodes,nnodes);%distance between each node
node_weight=1/nnodes*ones(nnodes,1);%weight allocated to each node

% place nodes randomly in area - will adjust to be method in report later
randcoords=sidelength*rand(nnodes,2);
% link with pair
Q=false(nnodes,nnodes);%1 if element j,k provides useful signal between node j and k
for i=1:nnodes/2
    Q(i,i+nnodes/2)=1;
    Q(i+nnodes/2,i)=1;
end
Qbar=(~Q);
%calulate distance between each node
for i=1:nnodes
    for j=1:nnodes
        dist_matrix(i,j)=norm(randcoords(i,:)-randcoords(j,:));
    end
end

%define fading coefficient using rayleigh fading, sampling from standard normal dist
const=(physconst('LightSpeed')/(4*pi*carr_freq))^2*node_power;
E=const.*((dist_matrix/1000).^-beta);% - mean of received signal with fading (depends on distance and frequency)
fade_coeff=raylrnd(E);%random coefficient based on Rayleigh distribution with mean E
gamma=node_power.*(abs(fade_coeff).^2); %received signal power from nodes j to node k
%Add self interference with SI suppression
for i=1:nnodes
    gamma(i,i)=sum(gamma(isfinite(gamma(:,i))))/SIsuppression;
end

D=gamma.*Q;%matrix of desired signal powers
I=gamma.*Qbar;%matrix of undesired signal powers
%% scheduling schemes
% %Scheme 1: Centralised scheme similar to OC-DTDD
exit=true;
lambda=4;prev_lambda=1;
%Initiate t, w and l randomly
n=0;A=zeros(nnodes,1);B=zeros(nnodes,1);
t_state=round(rand(nnodes,1));w_state=round(rand(nnodes,1));l_state=round(rand(nnodes,1));
t_state_sum=zeros(nnodes,1);
counter=1;
% %iterative loop to calculate t_state
while exit
    if (abs(lambda-prev_lambda)<1)
        exit=false;
    else
        n=n+1;
        %calculate tx_sum for each node combination and sum to get tx for
        % each node
        t_state_prev=t_state;
        for i=1:nnodes
            for k=1:nnodes
                t_state_sum(k)=(node_weight(k)*l_state(k)/log(2))*(D(i,k)*(1-w_state(k)/sqrt(t_state_prev.'*D(:,k)))+I(i,k));
                if isnan(t_state_sum(k))||isinf(t_state_sum(k))
                    t_state_sum(k)=0;
                end
            end
            if sum(t_state_sum)>0.0000000001
                t_state(i)=0;
            else
                t_state(i)=1;
            end
        end
        %calculate auxiliary variables A and B along with t, w and l states
        for i=1:nnodes
            A(i)=t_state.'*D(:,i);
            B(i)=1+t_state.'*I(:,i);
            w_state(i)=(A(i)+B(i))/sqrt(A(i));
            l_state(i)=(abs(sqrt(A(i))-w_state(i))^2+B(i))^-1;
        end
    end
    lambda=lambda-1;
end
%determine the optimal state of each node based on state of the matched pair
opt_state=false(nnodes,4);%matrix of the state of all nodes
for i=1:floor(nnodes/2)
    if (t_state(i)==1)&&(t_state(i+nnodes/2)==1)
        opt_state(i,3)=1;
        opt_state(i+nnodes/2,3)=1;
    elseif (t_state(i)==0)&&(t_state(i+nnodes/2)==0)
        opt_state(i,4)=1;
        opt_state(i+nnodes/2,4)=1;
    elseif (t_state(i)==1)&&(t_state(i+nnodes/2)==0)
        opt_state(i,2)=1;
        opt_state(i+nnodes/2,1)=1;
    elseif (t_state(i)==0)&&(t_state(i+nnodes/2)==1)
        opt_state(i,1)=1;
        opt_state(i+nnodes/2,2)=1;
    end
end



%Scheme 2: conventional setup where  2 random groups transmit at different time slots
conv_state=false(nnodes,4);
for  i=1:nnodes
   random_number=round(rand);
   if random_number==1
       conv_state(i,3)=1;
   else 
       conv_state(i,2)=1;
   end
end

%Scheme 3: random state allocation
rand_state=false(nnodes,4);
for  i=1:nnodes
   random_number=rand;
   if random_number<=0.25
       rand_state(i,1)=1;
   elseif random_number<=0.5
       rand_state(i,2)=1;
   elseif random_number<=0.75
       rand_state(i,3)=1;
   else 
       rand_state(i,4)=1;
   end
end

%Scheme 4: duplex scheme every node in duplex mode all the time
duplex_state=false(nnodes,4);
duplex_state(:,3)=ones(nnodes,1);
%% caclulate results for each scheduling scheme
%calulate SINR
r1=opt_state(:,1);t1=opt_state(:,2);f1=opt_state(:,3);s1=opt_state(:,4);
r2=conv_state(:,1);t2=conv_state(:,2);f2=conv_state(:,3);s2=conv_state(:,4);
r3=rand_state(:,1);t3=rand_state(:,2);f3=rand_state(:,3);s3=rand_state(:,4);
r4=duplex_state(:,1);t4=duplex_state(:,2);f4=duplex_state(:,3);s4=duplex_state(:,4);
SINR1=zeros(nnodes,1);
SINR2=zeros(nnodes,1);
SINR3=zeros(nnodes,1);
SINR4=zeros(nnodes,1);
%Caluculate SINR for each node
for  i=1:nnodes
    SINR1(i)=(r1(i)+f1(i))*((t1+f1).'*D(:,i))/(1+(t1+f1).'*I(:,i));
    SINR2(i)=(r2(i)+f2(i))*((t2+f2).'*D(:,i))/(1+(t2+f2).'*I(:,i));
    SINR3(i)=(r3(i)+f3(i))*((t3+f3).'*D(:,i))/(1+(t3+f3).'*I(:,i));
    SINR4(i)=(r4(i)+f4(i))*((t4+f4).'*D(:,i))/(1+(t4+f4).'*I(:,i));
end

%calculate (weighted) sum rate
SumRate1=node_weight'*log2(1+SINR1);
SumRate2=node_weight'*log2(1+SINR2);
SumRate3=node_weight'*log2(1+SINR3);
SumRate4=node_weight'*log2(1+SINR4);
