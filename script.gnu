set term wxt persist
set title '����������� V(t)'
set xlabel 't'
set ylabel 'V'
plot 'data.txt' using 1:2 with points title '�������� ������', 'data.txt' using 1:3 with lines title '������������� ������', 'data_fit.txt' using 1:2 with lines title '����������� �����������'
