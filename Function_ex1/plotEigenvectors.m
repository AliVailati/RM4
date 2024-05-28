function plotEigenvectors(V_1year, V_3year, V_5year, V_7year, ii)
    figure
    plot(V_1year(:,ii),'-o','DisplayName','1 year')
    hold on
    plot(V_3year(:,ii),'-square','DisplayName','3 year')
    plot(V_5year(:,ii),'-diamond','DisplayName','5 year')
    plot(V_7year(:,ii), '-^', 'DisplayName', '7 year')
    grid on
    title (sprintf('Eigenvectors # %d', ii))
    legend('1year', '3years', '5years', '7years')
    hold off