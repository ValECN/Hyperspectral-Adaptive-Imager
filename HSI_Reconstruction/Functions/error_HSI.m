function y = error_HSI(HSI_1, HSI_2, mode)

switch mode
    case 'l2'
        y = sum(sum(sum((HSI_1-HSI_2).^2)));
    case 'l1'
        y = sum(sum(sum(abs(HSI_1-HSI_2))));
end