% Date Created: 7/12/2022
% Created by Anuradha Agarwal

%Time span
T=50;
daily = 1:1:T; %Daily
tspan=daily;

%Initial Guess
Fitted_Parameters=[0.001 0.03]; 

%Initial Conditions
N = 1000;
I0 = 10;
R0 = 0;
S0 = 990;

beta = 0.001;
alpha = 0.03;
init_cond = [S0, I0, R0];

Number_Parameters=length(Fitted_Parameters); %Number of parameters.

[t,y_est]=ode45(@SIR_model_states,daily,init_cond,[],Fitted_Parameters);
Prevalencedata = y_est(:,2) + 0.01*(y_est(:,2));

%% MC Optimization
betavalues = linspace(0.0008, 0.0012, 100);
length_of_beta = length(betavalues);
EstiParam = zeros(length_of_beta, Number_Parameters); %Storage for estimated parameters.
function_values = zeros(1, length_of_beta);
beta_val = zeros(1,length_of_beta);
alpha_val = zeros(1,length_of_beta);
j = 0.05;
for i = 1:length_of_beta
    Initial_Guess = j;
    
    %Bounds for the parameters
    Lowerbounds = [0];
    Upperbounds=[1];
    fixp = betavalues(i);
    options=optimset('Disp','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',15000,'MaxFunEval',15000);
    [EstimatedParameters,fval,exitflag]=fmincon(@(Initial_Guess)err_in_data_SIR(Initial_Guess,Prevalencedata, T,fixp),Initial_Guess,[],[],[],[],Lowerbounds, Upperbounds, [], options);
    EstiParam(i,:) = [fixp,EstimatedParameters]; %Stores parameters for current noise level to computer ARE.
    function_values(1,i) = fval;
    beta_val(1,i) = fixp;
    alpha_val(1,i) = EstimatedParameters(1);
    j = EstimatedParameters;
end

threshold = chi2inv(0.95,2) + min(function_values);
figure
plot(betavalues, function_values)
hold on 
yline(threshold, 'r','Threshold')
hold off 

%%  SEIR model
function [dx]  = SIR_model_states(t,x,z)

%parameters to be estimated
beta = z(1);
alpha = z(2);

S = x(1); I = x(2); R = x(3);
dx = zeros(3,1);
dx(1) = -beta*S*I;
dx(2) = beta*S*I-alpha*I;
dx(3) = alpha*I;
end

%%

function  error_in_data = err_in_data_SIR(par,data,T,fp)
%This is the error function for prevalence.
%Initial Conditions
N = 1000;
I0 = 10;
R0 = 0;
S0 = N - I0 - R0;
init_cond = [S0,I0,R0];
daily = 1:1:T; %Daily
tspan=daily;
par1=[fp,par];
[t,y]=ode45(@SIR_model_states,daily,init_cond,[],par1);

Infected=y(tspan(:),2);
weight=1./(0.1* (1+y(:,2)));
error_in_data = sum((weight.*(Infected-data)).^2);

end



