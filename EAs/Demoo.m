clc,clear
addpath(genpath(pwd))

POP_SIZE = 50;
GENERATION = 20;
MAX_FEs = POP_SIZE * GENERATION;
M = 9;
K = 5;

D = 20;
dim = K*(M+1);

% for DE
P_BEST = 0.11;
MEM_SIZE = 6;  
ARC_RATE = 2.6;

UB = 3;
Results = zeros(UB,6);
GA = zeros(UB,3);
DE = zeros(UB,3);
BGA = zeros(UB,3);
BGA2 = zeros(UB,3);
BGA3 = zeros(UB,3);
BGA_byGA = zeros(UB,3);

cvg1 = zeros(UB,GENERATION);
cvg2 = zeros(UB,GENERATION);
cvg3 = zeros(UB,GENERATION);
cvg4 = zeros(UB,GENERATION);

for IterLarge = 1:UB
    disp(IterLarge);
    AP=unifrnd(-D/2,D/2,M,2);
    Ter=unifrnd(-D/2,D/2,K,2);
    
    obj = @(x)sum(x);
    obj_func = @(x)Rate(x,AP,Ter,M,K,obj);
    obj_func_bga = @(x)-Rate_bga(x,AP,Ter,M,K,obj);
    
    Results(IterLarge,5) = Rate1(AP,Ter,M,K,obj);
    
    %% call algorithm to maximize the objective function
    % Binary GA supported by GA
    %[best_solu, cvg4(IterLarge,:), best, GA_best] = run_b_ga_sp_by_GA(dim, obj_func_bga, obj_func, MAX_FEs, POP_SIZE);
    %Results(IterLarge,4) = best;
        % disp(sprintf('%.2f ', ratio(best_solu,K,M)));
        %BGA_byGA(IterLarge,:) = ratio(best_solu,K,M);
        
    % Binary GA
    [best_solu, cvg2(IterLarge,:), best] = run_b_ga(dim, obj_func_bga, MAX_FEs, POP_SIZE);
    Results(IterLarge,2) = best;
        % disp(sprintf('%.2f ', ratio(best_solu,K,M)));
         BGA(IterLarge,:) = ratio(best_solu,K,M);
         
    % Binary GA one crossover
    %fprintf("one crossover \n");
    %[best_solu, cvg1(IterLarge,:), best] = run_b_ga_1(dim, obj_func_bga, MAX_FEs, POP_SIZE);
    %Results(IterLarge,1) = best;
        % disp(sprintf('%.2f ', ratio(best_solu,K,M)));
     %    BGA2(IterLarge,:) = ratio(best_solu,K,M);
         
    % Binary GA two crossover
    %fprintf("two crossover \n");
    %[best_solu, cvg3(IterLarge,:), best] = run_b_ga_2(dim, obj_func_bga, MAX_FEs, POP_SIZE);
    %Results(IterLarge,3) = best;
        % disp(sprintf('%.2f ', ratio(best_solu,K,M)));
         %BGA3(IterLarge,:) = ratio(best_solu,K,M);
        
    % DE
    [best_solu, cvg1(IterLarge,:), best] = run_ide(dim, obj_func, MAX_FEs, POP_SIZE, P_BEST, MEM_SIZE, ARC_RATE, true);
    Results(IterLarge,1) = best;
         disp(sprintf('%.2f ', ratio(best_solu,K,M)));
        DE(IterLarge,:) = ratio(best_solu,K,M);
        
    % GA
    [best_solu, cvg3(IterLarge,:), best] = run_ga(dim, obj_func, MAX_FEs, POP_SIZE, GENERATION);
    Results(IterLarge,3) = best;
    %if GA_best > best
    %    Results(IterLarge,3) = GA_best;
    %    cvg3(IterLarge,GENERATION) = GA_best;
    %    disp("Best objective = " + GA_best);

         disp(sprintf('%.2f ', ratio(best_solu,K,M)));
        GA(IterLarge,:) = ratio(best_solu,K,M);
    
    %[~, index] = max(Results(IterLarge,1:4));
    %Results(IterLarge,6) = index;
  
    
    %fprintf("    DE = %d, BGA = %d, GA = %d, IBGA = %d\n", sum(Results(:,6) == 1), sum(Results(:,6) == 2), sum(Results(:,6) == 3), sum(Results(:,6) == 4));
    %fprintf('-----------------------------------------------------------------------\n');
    
    %{
    if IterLarge > UB - 2
        figure;
        x = 1 : GENERATION;
        plot(x,cvg1(IterLarge,:), 'r-', 'LineWidth', 2);
        hold on;
        plot(x,cvg2(IterLarge,:), 'g--', 'LineWidth', 2);
        plot(x,cvg3(IterLarge,:), 'b:', 'LineWidth', 2);

        title(sprintf('Đồ thị thứ %d', IterLarge));
        xlabel('Generation');
        ylabel('Value');
        legend('DE', 'BGA', 'GA');
        hold off;
    end
    
    
    POP_SIZE = POP_SIZE + 3;
    MAX_FEs = POP_SIZE * GENERATION;
    %}
    end
        %
        figure;
        x = 1 : GENERATION;
        plot(mean(cvg1(1:UB,x)), 'Color', 'r', 'LineWidth', 2);
        hold on;
        plot(mean(cvg2(1:UB,x)), 'Color', 'g', 'LineWidth', 2);
        plot(mean(cvg3(1:UB,x)), 'Color', 'b', 'LineWidth', 2);
        %plot(mean(cvg4(1:UB,x)), 'Color', 'm', 'LineWidth', 2);

        xlabel('Generation');
        ylabel('Value');
        title(sprintf("Trung bình cộng: M=%d, K=%d, Total generation = %d, Iterlarge = %d", M, K, GENERATION,UB));
        legend('DE', 'BGA', 'GA');
        hold off;
     

        %{
        figure;
        x = 1 : GENERATION;
        plot(geomean(cvg1(1:UB,x)), 'r-', 'LineWidth', 2);
        hold on;
        plot(geomean(cvg2(1:UB,x)), 'g--', 'LineWidth', 2);
        plot(geomean(cvg3(1:UB,x)), 'b:', 'LineWidth', 2);
        xlabel('Generation');
        ylabel('Value');
        title(sprintf("Trung bình nhân: M=%d, K=%d, Total generation = %d, Iterlarge = %d", M, K, GENERATION,UB));
        legend('DE', 'BGA', 'GA');
        hold off;      
        

        fprintf("M=%d, K=%d, Total generation = %d, Iterlarge = %d\n", M, K, MAX_FEs/POP_SIZE,UB);
        fprintf("    DE       BGA       GA        IBGA       alpha = 1\n");
        disp(Results(:,1:5));

fprintf("    DE = %d, BGA = %d, GA = %d, IBGA = %d\n", sum(Results(:,6) == 1), sum(Results(:,6) == 2), sum(Results(:,6) == 3), sum(Results(:,6) == 4));
        %}

%{
fprintf("    DE - Satellite only : %f\n", mean(DE(:,1)));
fprintf("    DE - APs only       : %f\n", mean(DE(:,2)));
fprintf("    DE - both           : %f\n\n", mean(DE(:,3)));

fprintf("    BGA - Satellite only: %f\n", mean(BGA(:,1)));
fprintf("    BGA - APs only      : %f\n", mean(BGA(:,2)));
fprintf("    BGA - both          : %f\n\n", mean(BGA(:,3)));

fprintf("    GA - Satellite only : %f\n", mean(GA(:,1)));
fprintf("    GA - APs only       : %f\n", mean(GA(:,2)));
fprintf("    GA - both           : %f\n\n", mean(GA(:,3)));
%}

    %{
    %% draw graph

    figure; 
    hold on; 
    plot(AP(:, 1), AP(:, 2), 'ro');  
    plot(Ter(:, 1), Ter(:, 2), 'bs'); 
    xlabel('X-axis');  
    ylabel('Y-axis'); 
    title('Users and Access points');
    legend('Access points', 'Users');
    hold off;
    %}