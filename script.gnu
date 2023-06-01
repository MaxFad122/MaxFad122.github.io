set term wxt persist
set title 'Зависимость V(t)'
set xlabel 't'
set ylabel 'V'
plot 'data.txt' using 1:2 with points title 'Реальные данные', 'data.txt' using 1:3 with lines title 'Предсказанные данные', 'data_fit.txt' using 1:2 with lines title 'Подобранная зависимость'
