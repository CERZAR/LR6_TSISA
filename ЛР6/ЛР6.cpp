#include "pch.h"
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <fcntl.h>
#include <io.h>

#define M_PI 3.1415926535897932384626433832795

struct Data
{
	double h,
		distanace,
		omega,
		delta,
		j;
	std::vector<double> alpha;
	std::vector<double> filtered_func;
};

struct Result
{
	std::vector<Data> database;
	std::vector<Data> answer;
	std::vector<int> r;
};

double func(double &x)
{
	return sin(x) + 0.5;
}

double get_random(double min, double max)
{
	return (double)rand() / RAND_MAX * (max - min) + min;
}

double get_tests_number() // N
{
	double p = 0.95,
		epsilon = 0.01,
		x_min = 0,
		x_max = M_PI;
	return 1 + log(1 - p) / log(1 - epsilon/(x_max - x_min));
}

std::vector<double> get_vector_k(int &k_min, int &k_max)
{
	double x_min = 0,
		x_max = M_PI;
	std::vector<double> vector;
	for (int i = k_min; i < k_max; i++)
	{
		vector.push_back((x_min + i * (x_max - x_min)) / k_max);
	}
	return vector;
}

std::vector<double> get_vector_sigma(int &k_min, int &k_max)
{
	double amplitude = 0.25;
	std::vector<double> vector;
	for (int i = k_min; i < k_max; i++)
	{
		vector.push_back(get_random(-1*amplitude, amplitude));
	}
	return vector;
}

void get_vector_func(int &k_min, int &k_max, std::vector<double> x_k)
{
	std::vector<double> vector;
	std::ofstream out("function.txt");
	if (out.is_open())
	{
		for (size_t i = k_min; i < k_max; i++)
		{
			out << "(" << x_k[i] << ";" << func(x_k[i]) << ") ";
		}
	}
	out.close();
}

std::vector<double> get_vector_noisy_func(int &k_min, int &k_max, std::vector<double> &x_k, std::vector<double> &sigma)
{
	std::vector<double> vector1;
	std::ofstream out("noisy_function.txt");
	std::vector<double> vector;
	for (int i = k_min; i < k_max; i++)
	{
		vector1.push_back(func(x_k[i]) + sigma[i]);
	}
	if (out.is_open())
	{
		for (int i = k_min; i < k_max; i++)
		{
			out << "(" << x_k[i] << ";" << vector1[i] << ") ";
		}
	}
	return vector1;
}

std::vector<double> get_vector_alpha(int r)
{
	std::vector<double> vector(r);
	//int up,
	//	down;
	//double current_double = 0;
	//double sum = 0;
	//up = down = r / 2;
	//while (up != r)
	//{
	//	if (up != r - 1)
	//	{
	//		current_double = get_random(0, 1 - sum);
	//		sum += current_double;
	//		vector[up] = vector[down] = current_double;
	//		up++;
	//		down--;
	//	}
	//	else
	//	{
	//		current_double = (1 - sum) / 2;
	//		vector[up] = vector[down] = current_double;
	//		up++;
	//	}
	//}
	if (r == 3)
	{
		vector[1] = get_random(0, 1);
		vector[0] = (1 - vector[1]) * 0.5;
		vector[2] = (1 - vector[1]) * 0.5;
	}
	else
	{
		vector[2] = get_random(0, 1);
		vector[1] = vector[3] = get_random(0, 1 - vector[2]) * 0.5;
		vector[0] = vector[4] = (1 - vector[1] - vector[2] - vector[3]) * 0.5;
	}
	return vector;
}

std::vector<double> get_vector_filtered_func(int &r, int &k_max, std::vector<double> &noisy_func, std::vector<double> &alpha) // arithmetic average
{
	int m = (r - 1) / 2;
	std::vector<double> vector1;
	double sum;
	for (int k = m; k < k_max - m; ++k)
	{
		sum = 0;
		for (int j = k - m; j < k + m + 1; ++j)
		{
			sum += noisy_func[j] * alpha[j + m - k];
		}
		vector1.push_back(sum);
	}
	return vector1;
}

double find_noisy_criterion(int k_max, std::vector<double> &filtered_func) // omega 
{
	double sum = 0;
	for (int k = 1; k < k_max; ++k)
	{
		sum += pow(filtered_func[k] - filtered_func[k - 1], 2);
	}
	return sqrt(sum);
}

double find_difference_criterion(int k_max, std::vector<double> &filtered_func, std::vector<double> &noisy_func) // delta 
{
	double sum = 0;
	for (int k = 0; k < k_max; ++k)
	{
		sum += pow(filtered_func[k] - noisy_func[k], 2);
	}
	return sqrt(sum / k_max);
}

double find_distance(double &omega, double &delta) // Euclidean metric
{
	return(sqrt(omega*omega + delta*delta));
}

Data find_convolution(double h, int &k_min, int &k_max, int &r, std::vector<double> &x_k, std::vector<double> &sigma, std::vector<double> noisy_func)
{
	Data current;
	current.alpha = get_vector_alpha(r);
	current.filtered_func = get_vector_filtered_func(r, k_max, noisy_func, current.alpha);
	current.omega = find_noisy_criterion(current.filtered_func.size(), current.filtered_func);
	current.delta = find_difference_criterion(current.filtered_func.size(), current.filtered_func, noisy_func);
	current.j = h * current.omega + (1 - h)*current.delta;
	return current;
}

void find_answer(Result &res)
{
	Data min;

	//for (size_t i = 0; i < res.r.size(); ++i)
	//{
	//	min.distanace = 9999999.;
	//	for (size_t j = 11*i; j < res.database.size() / res.r.size(); ++j)
	//	{
	//		if (res.database[j].distanace < min.distanace)
	//			min = res.database[j];
	//	}
	//	res.answer.push_back(min);
	//}
	min.distanace = 9999999.;
	for (size_t j = 0; j < 11; ++j)
	{
		if (res.database[j].distanace < min.distanace)
			min = res.database[j];
	}
	res.answer.push_back(min);
	min.distanace = 9999999.;
	for (size_t j = 11; j < res.database.size(); ++j)
	{
		if (res.database[j].distanace < min.distanace)
			min = res.database[j];
	}
	res.answer.push_back(min);
}

void researh(Result &res)
{
	Data data;
	int k_min = 0,
		k_max = 100,
		n = get_tests_number();
	double l_max = 10;
	//std::vector<int> r = { 3, 5 };
	std::vector<double> x_k = get_vector_k(k_min, k_max);
	std::vector<double> sigma = get_vector_sigma(k_min, k_max);
	std::vector<double> noisy_func = get_vector_noisy_func(k_min, k_max, x_k, sigma);

	double h;
	Data min;
	int r;
	for (int i = 0; i < 2; ++i)
	{
		if (i == 0) r = 3;
		else r = 5;
		for (int l = 0; l <= l_max; ++l)
		{
			h = l / l_max;
			min.j = 99999999.;
			for (int i = 0; i < n; ++i)
			{
				data = find_convolution(h, k_min, k_max, r, x_k, sigma, noisy_func);
				if (data.j < min.j)
					min = data;
			}
			min.h = h;
			min.distanace = find_distance(min.omega, min.delta);
			res.database.push_back(min);
		}
		res.r.push_back(r);
	}
	find_answer(res);
}

void output_begin(int r)
{
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x250C");
	if (r == 3)
	for (int i = 0; i < 96; i++)
	{
			if ((i == 13) || (i == 27) || (i == 41) || (i == 70) || (i == 83) || (i == 100)) wprintf(L"\x252C");
			else wprintf(L"\x2500");
	}
	else
	for (int i = 0; i < 112; i++)
	{
		if ((i == 13) || (i == 27) || (i == 41) || (i == 86) || (i == 99) || (i == 115)) wprintf(L"\x252C");
		else wprintf(L"\x2500");
	}
	wprintf(L"\x2510\n\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << "  Lambda(h)  ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << "      J      ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << "  Distance   ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	if (r == 3)std::cout << "		Alpha	       ";
	else 	std::cout << "			Alpha		       ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << "   Noisy(w) ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << "Diff(delta) ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502\n");
	wprintf(L"\x251C");
	if (r == 3)
	for (int i = 0; i < 96; i++)
	{
			if ((i == 13) || (i == 27) || (i == 41) || (i == 70) || (i == 83) || (i == 100)) wprintf(L"\x253C");
			else wprintf(L"\x2500");
	}
	else
	{
		for (int i = 0; i < 112; i++)
		{
			if ((i == 13) || (i == 27) || (i == 41) || (i == 86) || (i == 99) || (i == 114)) wprintf(L"\x253C");
			else wprintf(L"\x2500");
		}
	}
	wprintf(L"\x2524\n");
}



void print_res(Result &res)
{
	//std::cout << "h	dis	alpha	omega	delta\n";
	//for (size_t i = 0; i < res.database.size() / 2; ++i)
	//{
	//	std::cout << std::fixed << std::setprecision(1) << res.database[i].h << "	" << std::fixed << std::setprecision(4) << res.database[i].distanace << "	";
	//	for (int k = 0; k < 3; ++k)
	//		std::cout << res.database[i].alpha[k] << "	";
	//	std::cout << res.database[i].omega << "	" << res.database[i].delta << "\n";
	//}
	//std::cout << "\n";
	//for (size_t i = 11; i < res.database.size(); ++i)
	//{
	//	std::cout << std::fixed << std::setprecision(1) << res.database[i].h << "	" << std::fixed << std::setprecision(4) << res.database[i].distanace << "	";
	//	for (int k = 0; k < 5; ++k)
	//		std::cout << res.database[i].alpha[k] << "	";
	//	std::cout << res.database[i].omega << "	" << res.database[i].delta << "\n";
	//}
	output_begin(3);
	for (size_t i = 0; i < 11; i++)
	{
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(7) << std::right << std::fixed << std::setprecision(1) << res.database[i].h << "      ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(9) << std::right << std::fixed << std::setprecision(4) << res.database[i].j << "    ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(9) << std::right << std::fixed << std::setprecision(5) << res.database[i].distanace << "    ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(3) << std::right << std::fixed << std::setprecision(4) << "[" << res.database[i].alpha[0] << ", " << res.database[i].alpha[1] << ", " << res.database[i].alpha[2] << "]  ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(9) << std::right << res.database[i].omega << "   ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(9) << std::right << res.database[i].delta << "   |\n";
	}
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2514");
	for (int i = 0; i < 96; i++)
	{
		if ((i == 13) || (i == 27) || (i == 41) || (i == 70) || (i == 83) || (i == 100)) wprintf(L"\x2534");
		else wprintf(L"\x2500");
	}
	wprintf(L"\x2518\n");

	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(7) << std::right << std::fixed << std::setprecision(1) << res.answer[0].h << "      ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(9) << std::right << std::fixed << std::setprecision(4) << res.answer[0].j << "    ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(9) << std::right << std::fixed << std::setprecision(5) << res.answer[0].distanace << "    ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(3) << std::right << std::fixed << std::setprecision(4) << "[" << res.answer[0].alpha[0] << ", " << res.answer[0].alpha[1] << ", " << res.answer[0].alpha[2] << "]  ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(9) << std::right << res.answer[0].omega << "   ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(9) << std::right << res.answer[0].delta << "   |  <- best\n\n\n";

	output_begin(5);
	for (size_t i = 12; i < res.database.size(); i++)
	{
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(7) << std::right << std::fixed << std::setprecision(1) << res.database[i].h << "      ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(9) << std::right << std::fixed << std::setprecision(4) << res.database[i].j << "    ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(9) << std::right << std::fixed << std::setprecision(5) << res.database[i].distanace << "    ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(3) << std::right << std::fixed << std::setprecision(4) << "[" << res.database[i].alpha[0] << ", " << res.database[i].alpha[1] << ", " << res.database[i].alpha[2] << ", " << res.database[i].alpha[3] << ", " << res.database[i].alpha[4] << "]  ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(9) << std::right << res.database[i].omega << "   ";
		_setmode(_fileno(stdout), _O_U16TEXT);
		wprintf(L"\x2502");
		_setmode(_fileno(stdout), _O_TEXT);
		std::cout << std::setw(9) << std::right << res.database[i].delta << "   |\n";
	}
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2514");
	for (int i = 0; i < 112; i++)
	{
		if ((i == 13) || (i == 27) || (i == 41) || (i == 86) || (i == 99) || (i == 114)) wprintf(L"\x2534");
		else wprintf(L"\x2500");
	}
	wprintf(L"\x2518\n");

	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(7) << std::right << std::fixed << std::setprecision(1) << res.answer[1].h << "      ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(9) << std::right << std::fixed << std::setprecision(4) << res.answer[1].j << "    ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(9) << std::right << std::fixed << std::setprecision(5) << res.answer[1].distanace << "    ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(3) << std::right << std::fixed << std::setprecision(4) << "[" << res.answer[1].alpha[0] << ", " << res.answer[1].alpha[1] << ", " << res.answer[1].alpha[2] << ", " << res.answer[1].alpha[3] << ", " << res.answer[1].alpha[4] << "]  ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(9) << std::right << res.answer[1].omega << "   ";
	_setmode(_fileno(stdout), _O_U16TEXT);
	wprintf(L"\x2502");
	_setmode(_fileno(stdout), _O_TEXT);
	std::cout << std::setw(9) << std::right << res.answer[1].delta << "   |  <- best\n\n\n";
}

void print_dots(Result res, int seed)
{
	int k_min = 0,
		k_max = 100;
	std::vector<double> x_k = get_vector_k(k_min, k_max);
	get_vector_func(k_min, k_max, x_k);
	std::ofstream out("3_filtered_function.txt");
	if (out.is_open())
	{
		for (size_t i = 0; i < res.answer[0].filtered_func.size(); i++)
		{
			out << "(" << x_k[i] << ";" << res.answer[0].filtered_func[i] << ") ";
		}
	}
	out.close();
	out.open("5_filtered_function.txt");
	if (out.is_open())
	{
		out << seed << "\n";
		for (size_t i = 0; i < res.answer[1].filtered_func.size(); i++)
		{
			out << "(" << x_k[i] << ";" << res.answer[1].filtered_func[i] << ") ";
		}
	}
	out.close();
	out.open("3_filtered_dots.txt");
	if (out.is_open())
	{
		for (size_t i = 0; i < 11; i++)
		{
			out << res.database[i].delta << " ; " << res.database[i].omega << "\n";
		}
	}
	out.close();
	out.open("5_filtered_dots.txt");
	if (out.is_open())
	{
		for (size_t i = 11; i < res.database.size(); i++)
		{
			out << res.database[i].delta << " ; " << res.database[i].omega << "\n";
		}
	}
	out.close();
}

int main()
{
	int seed = time(NULL);
	std::cout << seed << "\n";
	srand(seed);
	Result res;
	researh(res);
	print_res(res);
	print_dots(res, seed);
}