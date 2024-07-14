#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <locale.h>
#include <Windows.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>

using namespace std;

// слой
struct layer {
	double z_min, z_max; // координаты
	double sigma; // сигма
};
vector<layer> layers; // слои

// неоднородность (аномалия)
struct anomaly {
	double x_min, x_max, z_min, z_max; // координаты
	double sigma; // сигма
};
vector<anomaly> anomalies; // неоднородности (аномалии)

// узел
struct node {
	double x, z; // координаты узла
};
vector<node> grid; // сетка

// константные значения

const int anomalies_index = 2; // коэффициент для дробления сетки в местах неоднородностей
const int drob_index = 2; // коэффициент дляя дробления сетки

const double big_number = 1E+9; // большое число для учета 1 краевых
const double true_solution_at_the_boundary = 0; // значение на границе (0) 

const double m = 1.2566E-6; // Мю
const double lambda = 1 / m; // параметр для матрицы жесткости (1 / Мю)
const double omega = 1E+4; // частота

const int max_number_of_iterations = 10000; // максимальное число итераций
const double eps = 1E-15; // невязка

vector<vector<int>> elements; // конечные элементы

vector<double> b; // глобальный вектор правой части

int global_nodes_count; // количество глобальных узлов
int elements_count; // количество конечных элементов
double border_for_first_conditions; // граница на которой заданы 1 краевые
double current_source_border; // граница на которой задан ток

// глобальная матрица в разреженном формате
vector<double> di; // диагональные элементы матрицы
vector<double> al; // внедиагональные элементы нижнего треугольника (по строкам)
vector<double> au; // внедиагональные элементы верхнего треугольника (по столбцам)
vector<int> ja; // номера столбцов (строк) внедиагональных элементов
vector<int> ia; // вспомогательный вектор для портрета

vector<vector<double>> G; // матрица жесткости
vector<vector<double>> M; // матрица массы
vector<vector<double>> local_p, local_c; // матрицы для блочной сборки

// для ЛОС
vector<double> L; // нижний треугольник 
vector<double> U; // верхний треугольник
vector<double> D; // диагональ
vector<double> z;
vector<double> r;
vector<double> y;
vector<double> p;
vector<double> t;
vector<double> s;
vector<double> xp; // для погрешности
vector<double> x0; // начальное приближение/полученное решение

// синус компонента тока
double sin_component(double x, double z) {
	if (z == current_source_border)
		return 1;
	else
		return 0;
}

// косинус компонента тока
double cos_component(double x, double z) {
	if (z == current_source_border)
		return 1;
	else
		return 0;
}

// значение сигмы
double sigma(int element_number) {
	int ver1 = elements[element_number][0];
	int ver4 = elements[element_number][3];

	// учет аномалий
	for (int i = 0; i < anomalies.size(); i++) {
		if (grid[ver1].z >= anomalies[i].z_min && grid[ver4].z <= anomalies[i].z_max)
			if (grid[ver1].x >= anomalies[i].x_min && grid[ver4].x <= anomalies[i].x_max)
				return anomalies[i].sigma;
	}

	// учет слоев
	for (int i = 0; i < layers.size(); i++) {
		if (grid[ver1].z >= layers[i].z_min && grid[ver4].z <= layers[i].z_max)
				return layers[i].sigma;
	}

	// иначе воздух
	return 0;
}

// генерация сетки
void grid_generation() {

	vector<double> grid_x; // сетка по x
	vector<double> grid_z; // сетка по z

	double x_min, x_max, z_min, z_max; // координаты границ расчетной области
	double hx_start, hz_start; // начальный шаг
	double kx, kz; // коэффициент разрядки
	double x_start, z_start; // координаты начала разрядки
	int layers_count; // количество слоев
	int anomalies_count; // количество неоднородностей
	int x_segments_count, z_segments_count; // количество отрезков

	ifstream input("input.txt");

	input >> x_min >> x_max >> z_min >> z_max;
	input >> hx_start >> hz_start;
	input >> kx >> kz;
	input >> x_start >> z_start;


	input >> layers_count;
	layers.resize(layers_count);
	for (int i = 0; i < layers_count; i++) {
		input >> layers[i].z_min >> layers[i].z_max >> layers[i].sigma;
	}

	input >> anomalies_count;
	anomalies.resize(anomalies_count);
	for (int i = 0; i < anomalies_count; i++) {
		input >> anomalies[i].x_min >> anomalies[i].x_max >> anomalies[i].z_min >> anomalies[i].z_max >> anomalies[i].sigma;
	}

	double x = x_start;
	double z = z_start;
	double hx = hx_start;
	double hz = hz_start;

	// генерация сетки по х
	while (x < x_max) {
		grid_x.push_back(x);
		x += hx;
		hx *= kx;
	}
	grid_x.push_back(x_max);

	hx = hx_start;
	x = x_start - hx;
	while (x > x_min) {
		grid_x.push_back(x);
		hx *= kx;
		x -= hx;
	}
	grid_x.push_back(x_min);

	// генерация сетки по z
	while (z < z_max) {
		grid_z.push_back(z);
		z += hz;
		hz *= kz;
	}
	grid_z.push_back(z_max);

	hz = hz_start;
	z = z_start - hz;
	while (z > z_min) {
		grid_z.push_back(z);
		hz *= kz;
		z -= hz;
	}
	grid_z.push_back(z_min);

	// учет слоев
	for (int i = 0; i < layers_count; i++) {
		grid_z.push_back(layers[i].z_max);
		grid_z.push_back(layers[i].z_min);
	}

	// учет неоднородностей
	for (int i = 0; i < anomalies_count; i++) {
		x = anomalies[i].x_min;
		while (x < anomalies[i].x_max) {
			grid_x.push_back(x);
			x += abs(anomalies[i].x_max - anomalies[i].x_min) / anomalies_index;
		}
		grid_x.push_back(anomalies[i].x_max);
		//grid_x.push_back(anomalies[i].x_min);


		z = anomalies[i].z_min;
		while (z < anomalies[i].z_max) {
			grid_z.push_back(z);
			z += abs(anomalies[i].z_max - anomalies[i].z_min) / anomalies_index;
		}
		grid_z.push_back(anomalies[i].z_max);
		//grid_z.push_back(anomalies[i].z_min);
	}

	// сортировка
	sort(grid_x.begin(), grid_x.end());
	sort(grid_z.begin(), grid_z.end());

	// удаление дубликатов
	grid_x.erase(unique(grid_x.begin(), grid_x.end()), grid_x.end());
	grid_z.erase(unique(grid_z.begin(), grid_z.end()), grid_z.end());


	// дробление сетки
	if (drob_index == 2) {
		double drob_h;
		int size = grid_x.size() - 1;

		for (int i = 0; i < size; i++) {
			drob_h = (grid_x[i + 1] - grid_x[i]) / drob_index;
			for (int j = 1; j < drob_index; j++) 
				grid_x.push_back(grid_x[i] + drob_h * j);
		}

		size = grid_z.size() - 1;

		for (int i = 0; i < size - 1; i++) {
			drob_h = (grid_z[i + 1] - grid_z[i]) / drob_index;
			for (int j = 1; j < drob_index; j++)
				grid_z.push_back(grid_z[i] + drob_h * j);
		}
	}

	// сортировка
	sort(grid_x.begin(), grid_x.end());
	sort(grid_z.begin(), grid_z.end());

	// удаление дубликатов
	grid_x.erase(unique(grid_x.begin(), grid_x.end()), grid_x.end());
	grid_z.erase(unique(grid_z.begin(), grid_z.end()), grid_z.end());

	x_segments_count = grid_x.size() - 1;
	z_segments_count = grid_z.size() - 1;

	global_nodes_count = grid_x.size() * grid_z.size();
	grid.resize(global_nodes_count);

	// вывод сетки в файл
	ofstream grids("grid.txt");
	grids << global_nodes_count << endl;
	int i = 0;
	for (int iy = 0; iy < z_segments_count + 1; iy++)
		for (int ix = 0; ix < x_segments_count + 1; ix++) {
			grid[i].x = grid_x[ix];
			grid[i].z = grid_z[iy];
			grids << i << " " << grid[i].x << " " << grid[i].z << endl;
			i++;
		}

	// вывод элементов в файл
	elements_count = x_segments_count * z_segments_count;
	elements.resize(elements_count);
	ofstream elem("elem.txt");
	elem << elements_count << endl;
	int number = 0;
	for (int j = 0; j < z_segments_count; j++) 
		for (int i = 0; i < x_segments_count; i++) {
			elements[number].resize(4);
			elements[number][0] = i + j * (x_segments_count + 1);
			elements[number][1] = i + j * (x_segments_count + 1) + 1;
			elements[number][2] = i + (j + 1) * (x_segments_count + 1);
			elements[number][3] = i + (j + 1) * (x_segments_count + 1) + 1;
			elem << number << " " << elements[i][0] << " " << elements[i][1] << " " << elements[i][2] << " " << elements[i][3] << endl;
			number++;
		}

	// для учета 1 краевых
	border_for_first_conditions = z_min;
	current_source_border = z_max;
}

// портрет матрицы
void portrait() {
	vector<vector<int>> help_array(global_nodes_count); // список смежности
	vector<int> current_element(4); // вектор вершин одного элемента

	// формируем список смежности
	for (int i = 0; i < elements_count; i++) { // идем по конечным элементам
		for (int j = 0; j < 4; j++)
			current_element[j] = elements[i][j];
		reverse(current_element.begin(), current_element.end()); //разворачиваем список чтобы элементы в нем располагались по убывaнию
		for (int j = 0; j < 4; j++)
			for (int k = j + 1; k < 4; k++) {
				int flag = 1; // метка присутствия в списке смежности
				for (int p = 0; p < help_array[current_element[j]].size() && flag; p++)
					if (help_array[current_element[j]][p] == current_element[k]) 
						flag = 0; // если элемент уже есть в списке смежности
				if (flag) {
					help_array[current_element[j]].push_back(current_element[k]); // если в списке еще нет такого элемента, то добавляем
				}
			}
	}

	for (int i = 0; i < global_nodes_count; i++)
		sort(help_array[i].begin(), help_array[i].end()); // сортируем по возрастанию 

	ia.resize(2 * global_nodes_count + 1);

	ia[0] = 0;
	ia[1] = 0;
	ia[2] = 1; // за первый диагональный

	for (int i = 1; i < global_nodes_count; i++) {
		ia[2 * i + 1] = ia[2 * i] + help_array[i].size() * 2; // по 2 элемента за 1
		ia[2 * (i + 1)] = ia[2 * i + 1] + help_array[i].size() * 2 + 1; // по 2 элемента за 1 + диагональный
	}

	ja.resize(ia[2 * global_nodes_count]);
	au.resize(ia[2 * global_nodes_count]);
	al.resize(ia[2 * global_nodes_count]);
	di.resize(2 * global_nodes_count);

	for (int i = 1, k = 1; i < global_nodes_count; i++) {
		for (int j = 0; j < help_array[i].size(); j++) { // верхняя строка

			ja[k] = 2 * help_array[i][j];
			ja[k + 1] = 2 * help_array[i][j] + 1;
			k += 2;
		}
		for (int j = 0; j < help_array[i].size(); j++) { // нижняя строка без диагонального элемента

			ja[k] = 2 * help_array[i][j];
			ja[k + 1] = 2 * help_array[i][j] + 1;
			k += 2;
		}
		ja[k] = 2 * i; // диагональный элемент
		k++;
	}
}

// вычисление и сборка глобальной матрицы
void global_matrix_calculation() {
	int ver1, ver4; // вершины конечного элемента

	double x_min, z_min, x_max, z_max; // координаты вершин прямоугольника
	double hx, hz; // шаги для коэффициентов
	double kg1, kg2; // коэффициенты для матрицы жесткости
	double km; // коэффициент для матрицы массы

	int begin_index; // начальный индекс, с которого смотрим номера столбцов
	int index; // текущий индекс
	int required_value; // значение которое мы должны найти (глобальный номер узла)

	G.resize(4);
	M.resize(4);
	local_p.resize(4);
	local_c.resize(4);
	b.resize(2 * global_nodes_count);

	for (int i = 0; i < 4; i++) {
		G[i].resize(4);
		M[i].resize(4);
		local_p[i].resize(4);
		local_c[i].resize(4);
	}

	for (int element_number = 0; element_number < elements_count; element_number++) {
		ver1 = elements[element_number][0];
		ver4 = elements[element_number][3];

		x_min = grid[ver1].x;
		z_min = grid[ver1].z;
		x_max = grid[ver4].x;
		z_max = grid[ver4].z;

		hx = x_max - x_min;
		hz = z_max - z_min;

		kg1 = hx / hz;
		kg2 = hz / hx;
		km = hx * hz / 36;

		// вычисляем локальную матрицу жесткости
		G[0][0] = G[1][1] = G[2][2] = G[3][3] = kg2 / 3 + kg1 / 3;
		G[0][1] = G[1][0] = G[2][3] = G[3][2] = -kg2 / 3 + kg1 / 6;
		G[0][2] = G[2][0] = G[1][3] = G[3][1] = kg2 / 6 - kg1 / 3;
		G[0][3] = G[3][0] = G[1][2] = G[2][1] = -kg2 / 6 - kg1 / 6;

		// вычисляем локальную матрицу массы
		M[0][0] = M[1][1] = M[2][2] = M[3][3] = 4 * km;
		M[1][0] = M[0][1] = M[2][0] = M[0][2] = M[1][3] = M[3][1] = M[2][3] = M[3][2] = 2 * km;
		M[0][3] = M[1][2] = M[2][1] = M[3][0] = km;


		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				local_p[i][j] = lambda * G[i][j];
				local_c[i][j] = omega * sigma(element_number) * M[i][j];
			}
		}

		for (int i = 0; i < 4; i++) {
			// диагональные элементы за данный глобальный узел
			di[elements[element_number][i] * 2] += local_p[i][i];
			di[elements[element_number][i] * 2 + 1] += local_p[i][i];

			// с компонента диагонального элемента
			int c_diag_pos = ia[elements[element_number][i] * 2 + 2] - 1; // позиция в массиве al(au) 
			al[c_diag_pos] += local_c[i][i];
			au[c_diag_pos] -= local_c[i][i];

			
			// учет первой строки(столбца)
			begin_index = ia[elements[element_number][i] * 2];
			for (int j = 0; j < i; j++, begin_index++) {
				index = begin_index; 
				required_value = elements[element_number][j] * 2;
				while (ja[index] != required_value) {
					index++;
				}
				al[index] += local_p[i][j];
				al[index + 1] -= local_c[i][j];

				au[index] += local_p[i][j];
				au[index + 1] += local_c[i][j];
			}

			// учет второй строки(столбца)
			begin_index = ia[elements[element_number][i] * 2 + 1];
			for (int j = 0; j < i; j++, begin_index++) {
				index = begin_index;
				required_value = elements[element_number][j] * 2;
				while (ja[index] != required_value) {
					index++;
				}
				al[index] += local_c[i][j];
				al[index + 1] += local_p[i][j];

				au[index] -= local_c[i][j];
				au[index + 1] += local_p[i][j];
			}
		}
	}
}

// вычисление глобального вектора правой части
void global_b_calculation() {
	for (int i = 0; i < grid.size(); i++) {
		b[2 * i] = sin_component(grid[i].x, grid[i].z);
		b[2 * i + 1] = cos_component(grid[i].x, grid[i].z);
	}
}

// краевые условия первого типа для нижней границы (по условию задачи)
void first_boundary_conditions() {
	vector<int> nodes_with_first_conditions; // вершины на которых задано краевое условие первого типа
	// поиск всех вершин на заданной границе
	for (int i = 0; i < grid.size(); i++) {
		if (grid[i].z == border_for_first_conditions)
			nodes_with_first_conditions.push_back(i);
	}

	// учет 1 краевых
	for (int i = 0; i < nodes_with_first_conditions.size(); i++) {
		b[nodes_with_first_conditions[i] * 2] = big_number * true_solution_at_the_boundary;
		b[nodes_with_first_conditions[i] * 2 + 1] = big_number * true_solution_at_the_boundary;

		di[nodes_with_first_conditions[i] * 2] = big_number;
		di[nodes_with_first_conditions[i] * 2 + 1] = big_number;
	}
}

void LUsq_decomposition() {

	for (int i = 0; i < ia[2 * global_nodes_count]; i++) {
		L[i] = al[i];
		U[i] = au[i];
	}
	for (int i = 0; i < 2 * global_nodes_count; i++)
		D[i] = di[i];
	
	for (int i = 0; i < 2 * global_nodes_count; i++) { // идем по строкам
		int start_index_current_line = ia[i]; // с какого индекса в этой строке 
		int start_index_next_line = ia[i + 1]; // с какого индекса в следующей строке
		double diag_add = 0;
		for (int j = start_index_current_line; j < start_index_next_line; j++) { // по количеству элементов в строке (разность индексов) 
			double l_add = 0;
			double u_add = 0;
			int column_pos = ja[j]; // для верхнего треугольника определяем номер столбца в нижнем
			int start_column_index = ia[column_pos]; // по номеру столбца определяем индексы из ia для нашего столбца в верхнем треугольнике 
			int start_next_column_index = ia[column_pos + 1]; // и для следующего столбца чтобы умножить значение из нижнего(верхнего) треугольника на значение из верхнего (нижнего) треугольника
			int line_index = start_index_current_line; // начальный индекс элементов строки
			int column_index = start_column_index; // начальный индекс элементов столбца другого треугольника
			// для нашего треугольника line_index < j для прохода по всем элементам строки до текущего чтобы учесть все значениея в нашем треугольнике
			// для другого треугольника column_index < start_next_column_index для значения треугольника L(U) для учета всех чисел в столбце(строке) противоположного треугольника U(L)
			while (line_index < j && column_index < start_next_column_index)
				if (ja[line_index] == ja[column_index]) { // для одинаковых индексов в обоих треугольниках
					l_add += L[line_index] * U[column_index];
					u_add += L[column_index] * U[line_index];
					line_index++;
					column_index++;
				}
				else if (ja[line_index] > ja[column_index])
					column_index++;
				else
					line_index++;
			L[j] = L[j] - l_add;
			U[j] = U[j] - u_add;
			L[j] = L[j] / D[ja[j]];
			U[j] = U[j] / D[ja[j]];
			diag_add += U[j] * L[j];
		}
		D[i] -= diag_add;
		D[i] = sqrt(D[i]);
	} 
}

void LOS() {
	double alpha, betta, p_r, p_p, p_t, r_r, r_r0;

	r_r0 = 0;
	// начальное приближение
	for (int i = 0; i < global_nodes_count * 2; i++)
		x0[i] = 0;
	
	// y = A * x0
	for (int i = 0; i < 2 * global_nodes_count; i++)
		y[i] = 0;
	for (int i = 0; i < 2 * global_nodes_count; i++) { // идем по строкам
		for (int j = ia[i]; j < ia[i + 1]; j++) { // учет элементов строки умноженные на столбец вектора принцип как в сборке глобальной матрицы
			y[i] += al[j] * x0[ja[j]];
			y[ja[j]] += au[j] * x0[i];
		}
		y[i] += di[i] * x0[i]; // учет диагонального элемента
	}

	// r0 = f - A * x0 = f - y
	for (int i = 0; i < 2 * global_nodes_count; i++)
		r[i] = b[i] - y[i];
	
	// r0 = L^-1 * (f - A * x0) = L^-1 * r0
	r[0] = r[0] / D[0];
	for (int i = 1; i < 2 * global_nodes_count; i++) { // по строкам
		double s = 0;
		for (int j = ia[i]; j < ia[i + 1]; j++)
			s += L[j] * r[ja[j]];
		r[i] = r[i] - s;
		r[i] = r[i] / D[i];
	}	

	//////////////////////////////
	// (r, r)
	for (int j = 0; j < 2 * global_nodes_count; j++)
		r_r0 += r[j] * r[j];
	//////////////////////////////
	
	// z0 = U^-1 * r0
	for (int i = 0; i < 2 * global_nodes_count; i++)
		z[i] = r[i];
	z[2 * global_nodes_count - 1] = z[2 * global_nodes_count - 1] / D[2 * global_nodes_count - 1];
	for (int i = 2 * global_nodes_count - 1; i >= 1; i--) {
		for (int j = ia[i]; j < ia[i + 1]; j++)
			z[ja[j]] -= U[j] * z[i];
		z[i - 1] = z[i - 1] / D[i - 1];
	}
	
	// y = A * z0
	for (int i = 0; i < 2 * global_nodes_count; i++)
		y[i] = 0;
	for (int i = 0; i < 2 * global_nodes_count; i++) { // идем по строкам
		for (int j = ia[i]; j < ia[i + 1]; j++) { // учет элементов строки умноженные на столбец вектора принцип как в сборке глобальной матрицы
			y[i] += al[j] * z[ja[j]];
			y[ja[j]] += au[j] * z[i];
		}
		y[i] += di[i] * z[i]; // учет диагонального элемента
	}
	
	// p0 = L^-1 * (A * z0) = L^-1 * y
	p[0] = y[0] / D[0];
	for (int i = 1; i < 2 * global_nodes_count; i++) { // по строкам
		double s = 0;
		for (int j = ia[i]; j < ia[i + 1]; j++)
			s += L[j] * p[ja[j]];
		p[i] = y[i] - s;
		p[i] = p[i] / D[i];
	}
	for (int i = 1; i < max_number_of_iterations; i++) {
		p_r = 0;
		p_p = 0;
		p_t = 0;
		r_r = 0;

		// (p, r); (p, p)
		for (int j = 0; j < 2 * global_nodes_count; j++) {
			p_r += p[j] * r[j];
			p_p += p[j] * p[j];
		}
		
		//alpha = p_r / p_p;
		alpha = p_r / p_p;
		
		// x = x + alpha * z
		for (int j = 0; j < 2 * global_nodes_count; j++)
			x0[j] = x0[j] + alpha * z[j];

		// r = r - alpha * p
		for (int j = 0; j < 2 * global_nodes_count; j++)
			r[j] = r[j] - alpha * p[j];

		// t = U^-1 * r
		for (int j = 0; j < 2 * global_nodes_count; j++)
			t[j] = r[j];

		t[2 * global_nodes_count - 1] = t[2 * global_nodes_count - 1] / D[2 * global_nodes_count - 1];
		for (int k = 2 * global_nodes_count - 1; k >= 1; k--) {
			for (int j = ia[k]; j < ia[k + 1]; j++)
				t[ja[j]] -= U[j] * t[k];
			t[k - 1] = t[k - 1] / D[k - 1];
		}

		// y = A * U^-1 * r = A * t
		for (int j = 0; j < 2 * global_nodes_count; j++)
			y[j] = 0;

		for (int k = 0; k < 2 * global_nodes_count; k++) { // идем по строкам
			for (int j = ia[k]; j < ia[k + 1]; j++) { // учет элементов строки умноженные на столбец вектора принцип как в сборке глобальной матрицы
				y[k] += al[j] * t[ja[j]];
				y[ja[j]] += au[j] * t[k];
			}
			y[k] += di[k] * t[k]; // учет диагонального элемента
		}

		// t = L^-1 * A * U^-1 * r = L^-1 * y
		t[0] = y[0] / D[0];
		for (int k = 1; k < 2 * global_nodes_count; k++) { // по строкам
			double s = 0;
			for (int j = ia[k]; j < ia[k + 1]; j++)
				s += L[j] * t[ja[j]];
			t[k] = y[k] - s;
			t[k] = t[k] / D[k];
		}

		// betta = - (p, L^-1 * A * U^-1 * r) / (p, p) = - (p, t) / (p, p)
		for (int j = 0; j < 2 * global_nodes_count; j++)
			p_t += p[j] * t[j];
		betta = - p_t / p_p;

		// y = U^-1 * r
		for (int j = 0; j < 2 * global_nodes_count; j++)
			y[j] = r[j];
		y[2 * global_nodes_count - 1] = y[2 * global_nodes_count - 1] / D[2 * global_nodes_count - 1];
		for (int k = 2 * global_nodes_count - 1; k >= 1; k--) {
			for (int j = ia[k]; j < ia[k + 1]; j++)
				y[ja[j]] -= U[j] * y[k];
			y[k - 1] = y[k - 1] / D[k - 1];
		}

		// z = U^-1 * r + betta * z = y + betta * z
		for (int j = 0; j < 2 * global_nodes_count; j++)
			z[j] = y[j] + betta * z[j];

		// p = L^-1 * A * U ^ -1 * r + betta * p = t + betta * p
		for (int j = 0; j < 2 * global_nodes_count; j++)
			p[j] = t[j] + betta * p[j];

		// (r, r)
		for (int j = 0; j < 2 * global_nodes_count; j++)
			r_r += r[j] * r[j];

		//////////////////////////////////////////////////
		//cout << i << ' ' << r_r << ' ' << r_r0 << ' ' << r_r / r_r0 << endl;

		if (r_r < eps)
			return;
		/////////////////////////////////////////////////

		/*if (i == 1) {
			r_r0 = r_r;
			cout << i << endl;
			cout << r_r0 << endl;
			if (r_r0 < eps)
				return;
		}
		else {
			cout << i << endl;
			cout << r_r << endl;
			cout << r_r / r_r0 << endl;

			if (r_r < eps)
				return;
		}*/
		
	}

	
}

void slay_solution() {
	L.resize(ia[2 * global_nodes_count]);
	U.resize(ia[2 * global_nodes_count]);
	D.resize(2 * global_nodes_count);
	x0.resize(2 * global_nodes_count);
	y.resize(2 * global_nodes_count);
	z.resize(2 * global_nodes_count);
	r.resize(2 * global_nodes_count);
	p.resize(2 * global_nodes_count);
	t.resize(2 * global_nodes_count);
	s.resize(2 * global_nodes_count);

	LUsq_decomposition();
	LOS();
}

vector<double> result_at_point(double x, double z) {
	int ver1, ver4; // вершины конечного элемента
	double x_min, z_min, x_max, z_max; // координаты вершин прямоугольника
	double hx, hz; // шаги
	vector<double> result = { 0, 0 };

	// определить, какому элементу принадлежит
	for (int element_number = 0; element_number < elements_count; element_number++) {
		ver1 = elements[element_number][0];
		ver4 = elements[element_number][3];

		x_min = grid[ver1].x;
		z_min = grid[ver1].z;
		x_max = grid[ver4].x;
		z_max = grid[ver4].z;

		hx = x_max - x_min;
		hz = z_max - z_min;

		if (x >= x_min && x <= x_max && z >= z_min && z <= z_max) {// вычисление соответствующих базисных функций

			vector<double> psi(4);
			psi[0] = (x_max - x) / hx * (z_max - z) / hz;
			psi[1] = (x - x_min) / hx * (z_max - z) / hz;
			psi[2] = (x_max - x) / hx * (z - z_min) / hz;
			psi[3] = (x - x_min) / hx * (z - z_min) / hz;

			result[0] = psi[0] * x0[2 * elements[element_number][0]] + psi[1] * x0[2 * elements[element_number][1]] + psi[2] * x0[2 * elements[element_number][2]] + psi[3] * x0[2 * elements[element_number][3]];
			result[1] = psi[0] * x0[2 * elements[element_number][0] + 1] + psi[1] * x0[2 * elements[element_number][1] + 1] + psi[2] * x0[2 * elements[element_number][2] + 1] + psi[3] * x0[2 * elements[element_number][3] + 1];
			return result;
		}
	}
	return result;
}

vector<double> calculation_Bs(double x, double z) {
	int ver1, ver4; // вершины конечного элемента
	double x_min, z_min, x_max, z_max; // координаты вершин прямоугольника
	double hx, hz; // шаги
	vector<double> Bs = { 0, 0 };

	// определить, какому элементу принадлежит
	for (int element_number = 0; element_number < elements_count; element_number++) {
		ver1 = elements[element_number][0];
		ver4 = elements[element_number][3];

		x_min = grid[ver1].x;
		z_min = grid[ver1].z;
		x_max = grid[ver4].x;
		z_max = grid[ver4].z;

		hx = x_max - x_min;
		hz = z_max - z_min;

		if (x >= x_min && x <= x_max && z >= z_min && z <= z_max) {// вычисление соответствующих базисных функций

			vector<double> psi(4);
			vector<double> dx(2);
			vector<double> dz(2);

			psi[0] = (x_max - x) / hx;
			psi[1] = (x - x_min) / hx;
			psi[2] = (z_max - z) / hz;
			psi[3] = (z - z_min) / hz;

			dx[0] = (double)-1 / hx;
			dx[1] = (double)1 / hx;

			dz[0] = (double)-1 / hz;
			dz[1] = (double)1 / hz;

			Bs[0] = dx[0] * psi[2] * x0[2 * elements[element_number][0]] + dx[1] * psi[2] * x0[2 * elements[element_number][1]] + dx[0] * psi[3] * x0[2 * elements[element_number][2]] + dx[1] * psi[3] * x0[2 * elements[element_number][3]];
			Bs[1] = dz[0] * psi[0] * x0[2 * elements[element_number][0]] + dz[0] * psi[1] * x0[2 * elements[element_number][1]] + dz[1] * psi[0] * x0[2 * elements[element_number][2]] + dz[1] * psi[1] * x0[2 * elements[element_number][3]];
			return Bs;
		}
	}
	return Bs;
}

vector<double> calculation_Bс(double x, double z) {
	int ver1, ver4; // вершины конечного элемента
	double x_min, z_min, x_max, z_max; // координаты вершин прямоугольника
	double hx, hz; // шаги
	vector<double> Bc = { 0, 0 };

	// определить, какому элементу принадлежит
	for (int element_number = 0; element_number < elements_count; element_number++) {
		ver1 = elements[element_number][0];
		ver4 = elements[element_number][3];

		x_min = grid[ver1].x;
		z_min = grid[ver1].z;
		x_max = grid[ver4].x;
		z_max = grid[ver4].z;

		hx = x_max - x_min;
		hz = z_max - z_min;

		if (x >= x_min && x <= x_max && z >= z_min && z <= z_max) {// вычисление соответствующих базисных функций

			vector<double> psi(4);
			vector<double> dx(2);
			vector<double> dz(2);

			psi[0] = (double)((x_max - x) / hx);
			psi[1] = (double)((x - x_min) / hx);
			psi[2] = (double)((z_max - z) / hz);
			psi[3] = (double)((z - z_min) / hz);

			dx[0] = (double)-1 / hx;
			dx[1] = (double)1 / hx;

			dz[0] = (double)-1 / hz;
			dz[1] = (double)1 / hz;

			Bc[0] = dx[0] * psi[2] * x0[2 * elements[element_number][0] + 1] + dx[1] * psi[2] * x0[2 * elements[element_number][1] + 1] + dx[0] * psi[3] * x0[2 * elements[element_number][2] + 1] + dx[1] * psi[3] * x0[2 * elements[element_number][3] + 1];
			Bc[1] = dz[0] * psi[0] * x0[2 * elements[element_number][0] + 1] + dz[0] * psi[1] * x0[2 * elements[element_number][1] + 1] + dz[1] * psi[0] * x0[2 * elements[element_number][2] + 1] + dz[1] * psi[1] * x0[2 * elements[element_number][3] + 1];
			return Bc;
		}
	}
	return Bc;
}

void result_in_receivers() {
	vector<node> receivers = {
		/*{-100, 0},
		{100, 0},
		{300, 0},
		{500, 0},
		{700, 0},
		{900,0},
		{1100, 0},
		{1300, 0},
		{1500, 0},
		{1700, 0},
		{1900,0},
		{1100, 0},
		{1300, 0},
		{1500, 0},
		{1700, 0},
		{1900,0}*/
		//{10,-10},
		{0, 0},
		{10, 0},
		{20, 0},
		{30, 0},
		{40, 0},
		{50, 0},
		{60, 0},
		{70, 0},
		{80, 0},
		{90, 0},
		{100, 0},
		{110, 0},
		{120, 0},
		{130, 0},
		{140, 0},
		{150, 0},
		{160, 0},
		{170, 0},
		{180, 0},
		{190, 0},
		{200, 0},
		{210, 0},
		{220, 0},
		{230, 0},
		{240, 0},
		{250, 0},
		{260, 0},
		{270, 0},
		{280, 0},
		{290, 0},
		{300, 0},
		{310, 0},
		{320, 0},
		{330, 0},
		{340, 0},
		{350, 0},
		{360, 0},
		{370, 0},
		{380, 0},
		{390, 0},
		{400, 0},
		{410, 0},
		{420, 0},
		{430, 0},
		{440, 0},
		{450, 0},
		{460, 0},
		{470, 0},
		{480, 0},
		{490, 0},
		{500, 0},
		{510, 0},
		{520, 0},
		{530, 0},
		{540, 0},
		{550, 0},
		{560, 0},
		{570, 0},
		{580, 0},
		{590, 0},
		{600, 0},
		{610, 0},
		{620, 0},
		{630, 0},
		{640, 0},
		{650, 0},
		{660, 0},
		{670, 0},
		{680, 0},
		{690, 0},
		{700, 0},
		{710, 0},
		{720, 0},
		{730, 0},
		{740, 0},
		{750, 0},
		{760, 0},
		{770, 0},
		{780, 0},
		{790, 0},
		{800, 0},
		{810, 0},
		{820, 0},
		{830, 0},
		{840, 0},
		{850, 0},
		{860, 0},
		{870, 0},
		{880, 0},
		{890, 0},
		{900, 0},
		{910, 0},
		{920, 0},
		{930, 0},
		{940, 0},
		{950, 0},
		{960, 0},
		{970, 0},
		{980, 0},
		{990, 0},
		{1000, 0},
	};
	ofstream file_A("A.txt");
	ofstream file_E("E.txt");

	ofstream file_B("B.txt");

	ofstream file_Z("Z.txt");

	vector<double> result_in_receiver = { 0, 0 }; // результат в приемнике

	double As, Ac; // компоненты вектора потенциалов
	double Es, Ec; // компоненты напряженности поля
	double E; // комплексный модуль напряженности 

	vector<double> Bs = { 0, 0 }; // компоненты индукции
	vector<double> Bc = { 0, 0 }; // компоненты индукции
	double B; // комплексный модуль индукции 

	double Z; // значение импеданса
	double p; // значение кажущегося удельного сопротивления
	double s; // полученное значение проводимости среды

	for (int i = 0; i < receivers.size(); i++) {
		result_in_receiver = result_at_point(receivers[i].x, receivers[i].z);

		As = result_in_receiver[0];
		Ac = result_in_receiver[1];

		Es = Ac * omega;
		Ec = - As * omega;
		E = sqrt(Es * Es + Ec * Ec);

		file_A << As << ' ' << Ac << endl;
		file_E << Es << ' ' << Ec << ' ' << E << endl;

		Bs = calculation_Bs(receivers[i].x, receivers[i].z);
		Bc = calculation_Bс(receivers[i].x, receivers[i].z);
		B = sqrt(Bs[0] * Bs[0] + Bs[1] * Bs[1] + Bc[0] * Bc[0] + Bc[1] * Bc[1]);

		file_B << Bs[0] << ' ' << Bs[1] << ' ' << Bc[0] << ' ' << Bc[1] << ' ' << B << endl;

		Z = m * E / B;
		p = Z * Z / (m * omega);
		s = 1 / p;
		

		file_Z << Z << ' ' << p << ' ' << s << endl;
	}
}

void main() {
	grid_generation();
	portrait();
	global_matrix_calculation();
	global_b_calculation();
	first_boundary_conditions();
	slay_solution();
	result_in_receivers();
}