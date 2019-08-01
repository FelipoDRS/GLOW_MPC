# GLOW_MPC
Files used in the paper Development of an MPC for stabilization and optimization of a gas lift oil well

README

Process model is well_disturbance.m
MPC model is MPC_model.m

*_der files only separate the derivatives

monitor files add noise and separate measuremnets available at the top of the well

MPC objective function is MPC_obj.m

Run first NN_script_MPC.m and later MPC_test_script.m. The first fits the neural network and the second tests the MPC

Data_genCOQ897 generates the data for parameter estimation 

fminsearch3.m is a modification of fminseach found in MatLab´s optimization toolbox

Just replace it for fminsearch and the code will run

O modelo do processo é o well_disturbance.m
o modelo do MPC é o MPC_model.m
os arquivos *_der apenas separam a derivada
Os arquivos monitor separam as variáveis medidas no topo da coluna e no caso do monitor_well_disturbance. adiciona ruído de medição
fminsearch3.m é uma modificação do fminsearch que retorna os resultados intermediários da otimização
MPC_obj.m é a função objetivo do MPC
NN_script_MPC.m cria e testa a rede neuronal, MPC_test_script.m implementa o MPC. 
É necessário rodar o primeiro antes do segundo.
Data_genCOQ897 gera os dados da estimação de parâmetros
