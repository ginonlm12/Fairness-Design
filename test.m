POP_SIZE = 10;
GENERATION = 40;
MAX_FEs = POP_SIZE * GENERATION;
M = 21;
K = 10;

D = 20;
dim = K*(M+1);

CostFunction = @(x)Rate(x,AP,Ter,M,K);
best = run_b_ga(CostFunction, dim, GENERATION, POP_SIZE);

