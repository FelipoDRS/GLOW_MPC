%this script is for the development of the NN that estimates the internal
%states
X0=[3.34	8.31	10.08];
par=[5.82e-3	14.45	3.29	1.20];
D_list=[];
Y=[];
sensor_data=[];
load('data_for_param_est.mat');
Ytrue=Y;
%N= 15000;%length(sensor_data);%number of data points used
data_WoS=sensor_data(1:17280,2:end); %data without slugging
data_WS=sensor_data(17281:end,2:end); %data without slugging

Y_mod=zeros(length(sensor_data),3);

ODS=odeset('MaxStep', 100);
X0(1)=X0(1)*1e3;
X0(2)=X0(2)*10;
X0(3)=X0(3)*1e3;

global u Dist

u=sensor_data(1,[end-1,end]);
Dist=sensor_data(1,end-2);
[t,ymodelo]=ode15s(@(t,y) MPC_der(t,y,par,u,Dist),[sensor_data(1,1),sensor_data(2,1)],X0,ODS);
Y_mod(1,:)=ymodelo(end,:);
%cheio de nomear variaveis

for i=2:length(sensor_data)
u=sensor_data(i,[end-1,end]);
Dist=sensor_data(i,end-2);
[~,ymodelo]=ode15s(@(t,y) MPC_der(t,y,par,u,Dist),[sensor_data(i-1,1),sensor_data(i,1)] ,ymodelo(end,:),ODS);
Y_mod(i,:)=ymodelo(end,:);

end

%past 60 sec average
W_av=zeros(length(data_WoS),4);
W_av(1:11,:)=data_WoS(1:11,1:4);

W_av2=zeros(length(data_WS),4);
W_av2(1:11,:)=data_WS(1:11,1:4);

for i=12:length(data_WoS)
W_av(i,:)=mean(data_WoS(i-11:i,1:4));
end
for i=12:length(data_WS)
W_av2(i,:)=mean(data_WS(i-11:i,1:4));
end


data_WoS=[W_av,data_WoS];
data_WS=[W_av2,data_WS];

Y_WoS=Y_mod(1:17280,[2,3]);
Y_WS=Y_mod(17281:end,[2,3]);

MNX=mean(data_WoS);
STDX=std(data_WoS);
X=(data_WoS-ones(size(data_WoS))*diag(MNX))./(ones(size(data_WoS))*diag(STDX));

MNY=mean(Y_WoS);
STDY=std(Y_WoS);
Y=(Y_WoS-ones(size(Y_WoS))*diag(MNY))./(ones(size(Y_WoS))*diag(STDY));

X_WS=(data_WS-ones(size(data_WS))*diag(MNX))./(ones(size(data_WS))*diag(STDX));

Y_WS=(Y_WS-ones(size(Y_WS))*diag(MNY))./(ones(size(Y_WS))*diag(STDY));
%Y_WStrue=Ytrue(17281:end,[2,3]);
%Y_WStrue=(Y_WStrue-ones(size(Y_WStrue))*diag(MNY))./(ones(size(Y_WStrue))*diag(STDY));

x = X';
t = Y';


x2 = X_WS';
t2 = Y_WS';
%t3= Y_WStrue';

%This use placeholding data while the estimation isnt good

X_ph=[x,x2];
T_ph=[t,t2]; %changed to the badly estimated data 

trainFcn = 'trainbfg';  % BFGS backpropagation.

% Create a Fitting Network
hiddenLayerSize = 6;
net = fitnet(hiddenLayerSize,trainFcn);

% Choose Input and Output Pre/Post-Processing Functions
% For a list of all processing functions type: help nnprocess

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = 30/100;
net.divideParam.valRatio = 35/100;
net.divideParam.testRatio = 35/100;

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error


% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
    'plotregression', 'plotfit'};
net.performParam.regularization = 1e-3;
% Train the Network
[net,tr] = train(net,X_ph,T_ph);
%[net,tr] = train(net,X,T);
%[net2,tr] = train(net,x2,t3);

% Test the Network
y = net(x);
e = t2-y2;
%performance = perform(net,t,y)

% Recalculate Training, Validation and Test Performance

f=@(x)( (sim(net,(x-MNX')./STDX')).*STDY'+MNY');

