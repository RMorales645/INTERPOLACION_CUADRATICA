%RONALD ALEXIS MORALES VARELA
%0901-23-6114
function [p, coeff, error_percent] = newton_interpolation(x_points, y_points, x_eval)
    % Interpolación polinomial de Newton con cálculo de error
    % Inputs:
    %   x_points: vector de puntos x conocidos
    %   y_points: vector de valores f(x) conocidos
    %   x_eval: punto donde evaluar el polinomio
    % Outputs:
    %   p: valor del polinomio en x_eval
    %   coeff: coeficientes del polinomio de Newton [b0, b1, b2,...]
    %   error_percent: error porcentual respecto al valor real
    
    n = length(x_points);
    if length(y_points) ~= n
        error('Los vectores x e y deben tener la misma longitud');
    end
    
    % Tabla de diferencias divididas
    F = zeros(n, n);
    F(:,1) = y_points(:);
    
    for j = 2:n
        for i = 1:n-j+1
            F(i,j) = (F(i+1,j-1) - F(i,j-1)) / (x_points(i+j-1) - x_points(i));
        end
    end
    
    % Coeficientes del polinomio (primera fila de F)
    coeff = F(1,:);
    
    % Evaluación del polinomio usando Horner
    p = coeff(n);
    for i = n-1:-1:1
        p = p .* (x_eval - x_points(i)) + coeff(i);
    end
    
    % Cálculo del error porcentual (si x_eval está en x_points)
    [~, idx] = ismember(x_eval, x_points);
    if idx > 0
        valor_real = y_points(idx);
        error_percent = abs((p - valor_real) / valor_real) * 100;
    else
        error_percent = NaN;
        warning('No se puede calcular error: x_eval no está en x_points');
    end
end

% Datos del problema
x = [1, 2, 3];
y = [0, 0.6931472, 1.098612];
x_eval = 2;  % Punto a evaluar (ln(2))

% Llamar a la función
[p, coeff, error] = newton_interpolation(x, y, x_eval);

% Mostrar resultados
fprintf('=== Interpolación de Newton para ln(2) ===\n');
fprintf('Puntos conocidos:\n');
disp([x' y']);
fprintf('\nCoeficientes del polinomio (b0, b1, b2):\n');
disp(coeff);
fprintf('\nEstimación en x=2: %.8f\n', p);
fprintf('Valor real de ln(2): 0.6931472\n');
fprintf('Error porcentual: %.10f%%\n', error);

% Graficar
xx = linspace(min(x), max(x), 100);
yy = zeros(size(xx));
for i = 1:length(xx)
    [yy(i), ~, ~] = newton_interpolation(x, y, xx(i));
end

figure;
plot(x, y, 'ro', 'MarkerSize', 10, 'LineWidth', 2); hold on;
plot(xx, yy, 'b-', 'LineWidth', 1.5);
plot(x_eval, p, 'k*', 'MarkerSize', 10);
xlabel('x'); ylabel('f(x)');
title('Interpolación de Newton para ln(x)');
legend('Puntos conocidos', 'Polinomio interpolador', 'Estimación en x=2');
grid on;