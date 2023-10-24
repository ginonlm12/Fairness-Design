n = 100000;
random1 = 0.3 * ones(1, n);
random2 = 0.7 * ones(1, n);

eta_c_values = [2, 200];

figure;

for i = 1:length(eta_c_values)
    eta_c = eta_c_values(i);

    child = sb_crossover(random1, random2, eta_c);

    offspring1 = child(1, :);
    offspring2 = child(2, :);
    data = [offspring1, offspring2];
    
    [f, xi] = ksdensity(data);

    f_normalized = f / sum(f);

    % Vẽ biểu đồ của ước tính mật độ xác suất (đã được chuẩn hóa)
    plot(xi, f, 'DisplayName', ['eta_c = ', num2str(eta_c)]);
    hold on;
end

legend;
xlabel('Giá trị');
ylabel('Ước tính Mật độ xác suất');
title('So sánh Ước tính Mật độ xác suất với các giá trị eta_c');
hold off;
