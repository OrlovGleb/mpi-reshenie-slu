# mpi-reshenie-slu
Решение слу с использованием mpi.
Входные данные программы:
  1) n – размерность матрицы,
  2) m – количество выводимых значений в матрице,
  3) k – задает номер формулы для инициализации матрицы, должен быть равен 0 при
вводе матрицы из файла
  4) filename – имя файла, откуда надо прочитать матрицу. Этот аргумент отсутствует,
если k! = 0.

Программа решает систему линейных уравнений, используя параллельное программирование на основе mpi.
