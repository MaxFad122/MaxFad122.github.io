set title '����������������� ������'
set xlabel 't'
set ylabel 'Z'
plot 'data.txt' with linespoints title '����������������� ������', \
     -8.09989*x**4 + 20.5123*x**3 + 0*x**2 + 0*x + 72.1822 title '����������������� ������', \
     -8.09989*x**2 + 20.5123*x + 0 title '����� �������'
pause -1 '������� ����� ������� ��� ������'
