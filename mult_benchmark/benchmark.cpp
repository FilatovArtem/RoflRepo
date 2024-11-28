#include <iostream>
#include <vector>
#include <chrono>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <cstring>

// Функция для вычисления статистических показателей
struct Statistics {
    double mean;
    double median;
    double stddev;
    double min;
    double max;

    Statistics(const std::vector<double>& times) {
        // Среднее значение
        mean = std::accumulate(times.begin(), times.end(), 0.0) / times.size();

        // Минимум и максимум
        auto [min_it, max_it] = std::minmax_element(times.begin(), times.end());
        min = *min_it;
        max = *max_it;

        // Стандартное отклонение
        double sq_sum = std::inner_product(times.begin(), times.end(), times.begin(), 0.0);
        stddev = std::sqrt(sq_sum / times.size() - mean * mean);

        // Медиана
        std::vector<double> sorted_times = times;
        std::sort(sorted_times.begin(), sorted_times.end());
        if (sorted_times.size() % 2 == 0) {
            median = (sorted_times[sorted_times.size()/2 - 1] + sorted_times[sorted_times.size()/2]) / 2;
        } else {
            median = sorted_times[sorted_times.size()/2];
        }
    }
};

// Функция для измерения времени выполнения
template<typename Func, typename... Args>
std::vector<double> benchmark_function(Func func, int num_runs, Args... args) {
    std::vector<double> times;
    times.reserve(num_runs);

    for (int i = 0; i < num_runs; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        func(args...);
        auto end = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double> duration = end - start;
        times.push_back(duration.count());
    }

    return times;
}

// Вывод статистики
void print_statistics(const std::string& name, const Statistics& stats) {
    std::cout << "\nСтатистика для " << name << ":\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Среднее время: " << stats.mean << " сек\n";
    std::cout << "Медиана: " << stats.median << " сек\n";
    std::cout << "Стд. отклонение: " << stats.stddev << " сек\n";
    std::cout << "Мин. время: " << stats.min << " сек\n";
    std::cout << "Макс. время: " << stats.max << " сек\n";
}

#define UNUSED(x) (void)(x)

// Функция для генерации случайной матрицы
void generate_matrix(double* matrix, int size) {
    for (int i = 0; i < size * size; ++i) {
        matrix[i] = (double)rand() / RAND_MAX * 200.0 - 100.0; // значения от -100 до 100
    }
}

// Функция для вычисления нормы матрицы
double calculate_norm(double* matrix, int size) {
    double norm = 0.0;
    for (int i = 0; i < size * size; ++i) {
        double abs_val = std::fabs(matrix[i]);
        if (abs_val > norm) norm = abs_val;
    }
    return norm;
}

// Здесь должны быть ваши две функции для сравнения
// Например:
void ilia_mult(double *a, double *b, double *res, int m1, int m2, int m3, int m, double norm) {
    int count_b = 0;
    UNUSED(m1);
    UNUSED(m2);
    UNUSED(m3);


    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            res[i * m + j] = 0.0;
            double abs_value = fabs(b[i * m + j]);
            if (1e+250 * norm < abs_value || abs_value < 1e-250 * norm) {
                b[i * m + j] = 0.0;
                ++count_b;
            }
        }
    }
    if (count_b == m * m) {
        return;
    }

    for (int k = 0; k < m; ++k) {
        for (int i = 2; i < m; i += 3) {
            double temp = a[i * m + k];
            double temp_1 = a[(i - 1) * m + k];
            double temp_2 = a[(i - 2) * m + k];

            for (int j = 2; j < m; j += 3) {
                double b_j = b[k * m + j];
                double b_jm1 = b[k * m + j - 1];
                double b_jm2 = b[k * m + j - 2];

                res[i * m + j] += b_j * temp;
                res[(i - 1) * m + j] += b_j * temp_1;
                res[(i - 2) * m + j] += b_j * temp_2;

                res[i * m + (j - 1)] += b_jm1 * temp;
                res[(i - 1) * m + (j - 1)] += b_jm1 * temp_1;
                res[(i - 2) * m + (j - 1)] += b_jm1 * temp_2;

                res[i * m + (j - 2)] += b_jm2 * temp;
                res[(i - 1) * m + (j - 2)] += b_jm2 * temp_1;
                res[(i - 2) * m + (j - 2)] += b_jm2 * temp_2;
            }

            for (int j = (m / 3) * 3; j < m; ++j) {
                double b_j = b[k * m + j];
                res[i * m + j] += b_j * temp;
                res[(i - 1) * m + j] += b_j * temp_1;
                res[(i - 2) * m + j] += b_j * temp_2;
            }
        }

        for (int i = (m / 3) * 3; i < m; ++i) {
            double temp = a[i * m + k];
            for (int j = 2; j < m; j += 3) {
                double b_j = b[k * m + j];
                double b_jm1 = b[k * m + j - 1];
                double b_jm2 = b[k * m + j - 2];

                res[i * m + j] += b_j * temp;
                res[i * m + (j - 1)] += b_jm1 * temp;
                res[i * m + (j - 2)] += b_jm2 * temp;
            }

            for (int j = (m / 3) * 3; j < m; ++j) {
                double b_j = b[k * m + j];
                res[i * m + j] += b_j * temp;
            }
        }
    }
}

inline void Artem_matrix_mult(double *a, double *b, double *c, int m, double norm)
{
    int zero_count_a = 0;
    int zero_count_b = 0;

    int i;
    int j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m; j++)
        {
            c[i * m + j] = 0.0;
            if (fabs(a[i * m + j]) < 1e-50 * norm)
            {
                a[i * m + j] = 0.;
                zero_count_a++;
            }
            if (fabs(b[i * m + j]) < 1e-50 * norm)
            {
                b[i * m + j] = 0.;
                zero_count_b++;
            }
        }
    }

    if (zero_count_a == m * m || zero_count_b == m * m) {
        return;
    }

    int t, q, r;
    double s00, s01, s02, s10, s11, s12, s20, s21, s22;
    int v3 = m % 3;
    int h3 = m % 3;
    for (r = 0; r < v3; r++) // номер строки
    {
        for (t = 0; t < h3; t++) // номер столбца
        {
            s00 = 0;
            for (q = 0; q < m; q++)
            {
                s00 += a[r * m + q] * b[q * m + t];
            }
            c[r * m + t] = s00;
        }
        for (; t < m; t += 3)
        {
            s00 = 0;
            s01 = 0;
            s02 = 0;
            for (q = 0; q < m; q++)
            {
                s00 += a[r * m + q] * b[q * m + t];
                s01 += a[r * m + q] * b[q * m + t + 1];
                s02 += a[r * m + q] * b[q * m + t + 2];
            }
            c[r * m + t] = s00;
            c[r * m + t + 1] = s01;
            c[r * m + t + 2] = s02;
        }
    }
    for (; r < m; r += 3)
    {
        for (t = 0; t < h3; t++)
        {
            s00 = 0;
            s10 = 0;
            s20 = 0;
            for (q = 0; q < m; q++)
            {
                s00 += a[r * m + q] * b[q * m + t];
                s10 += a[(r + 1) * m + q] * b[q * m + t];
                s20 += a[(r + 2) * m + q] * b[q * m + t];
            }
            c[r * m + t] = s00;
            c[(r + 1) * m + t] = s10;
            c[(r + 2) * m + t] = s20;
        }
        for (; t < m; t += 3)
        {
            s00 = 0;
            s01 = 0;
            s02 = 0;
            s10 = 0;
            s11 = 0;
            s12 = 0;
            s20 = 0;
            s21 = 0;
            s22 = 0;
            for (q = 0; q < m; q++)
            {
                s00 += a[r * m + q] * b[q * m + t];
                s01 += a[r * m + q] * b[q * m + t + 1];
                s02 += a[r * m + q] * b[q * m + t + 2];
                s10 += a[(r + 1) * m + q] * b[q * m + t];
                s11 += a[(r + 1) * m + q] * b[q * m + t + 1];
                s12 += a[(r + 1) * m + q] * b[q * m + t + 2];
                s20 += a[(r + 2) * m + q] * b[q * m + t];
                s21 += a[(r + 2) * m + q] * b[q * m + t + 1];
                s22 += a[(r + 2) * m + q] * b[q * m + t + 2];
            }
            c[r * m + t] = s00;
            c[r * m + t + 1] = s01;
            c[r * m + t + 2] = s02;
            c[(r + 1) * m + t] = s10;
            c[(r + 1) * m + t + 1] = s11;
            c[(r + 1) * m + t + 2] = s12;
            c[(r + 2) * m + t] = s20;
            c[(r + 2) * m + t + 1] = s21;
            c[(r + 2) * m + t + 2] = s22;
        }
    }
}

int main() {
    const int matrix_size = 300; // размер матрицы
    const int num_runs = 10; // количество прогонов в одной серии
    const int num_series = 100; // количество серий тестов
    
    std::vector<double> series_ratios; // для хранения отношений времени между сериями
    series_ratios.reserve(num_series);
    
    std::vector<double> all_times1, all_times2; // для хранения всех времен выполнения
    all_times1.reserve(num_series * num_runs);
    all_times2.reserve(num_series * num_runs);

    std::cout << "\nСравнение производительности умножения матриц " << matrix_size << "x" << matrix_size << "\n";
    std::cout << "Количество прогонов в серии: " << num_runs << "\n";
    std::cout << "Количество серий: " << num_series << "\n\n";

    // Выделяем память под матрицы
    double *a = new double[matrix_size * matrix_size];
    double *b = new double[matrix_size * matrix_size];
    double *a_copy = new double[matrix_size * matrix_size];  // Копии для второй функции
    double *b_copy = new double[matrix_size * matrix_size];
    double *res1 = new double[matrix_size * matrix_size];
    double *res2 = new double[matrix_size * matrix_size];

    for (int series = 0; series < num_series; ++series) {
        std::cout << "Серия " << series + 1 << ":\n";
        
        std::vector<double> times1, times2;
        times1.reserve(num_runs);
        times2.reserve(num_runs);

        for (int run = 0; run < num_runs; ++run) {
            // Генерируем новые матрицы для каждого прогона
            generate_matrix(a, matrix_size);
            generate_matrix(b, matrix_size);
            
            // Создаем копии матриц для второй функции
            std::memcpy(a_copy, a, matrix_size * matrix_size * sizeof(double));
            std::memcpy(b_copy, b, matrix_size * matrix_size * sizeof(double));
            
            double norm = calculate_norm(b, matrix_size);
            double norm_copy = norm; // Та же норма для второй функции

            // Замер времени для ilia_mult
            auto start1 = std::chrono::high_resolution_clock::now();
            ilia_mult(a, b, res1, matrix_size, matrix_size, matrix_size, matrix_size, norm);
            auto end1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration1 = end1 - start1;
            times1.push_back(duration1.count());
            all_times1.push_back(duration1.count());

            // Замер времени для Artem_matrix_mult
            auto start2 = std::chrono::high_resolution_clock::now();
            Artem_matrix_mult(a_copy, b_copy, res2, matrix_size, norm_copy);
            auto end2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration2 = end2 - start2;
            times2.push_back(duration2.count());
            all_times2.push_back(duration2.count());
            
            // Проверяем, что результаты совпадают
            if (run == 0 && series == 0) {
                bool results_match = true;
                for (int i = 0; i < matrix_size * matrix_size; ++i) {
                    if (std::fabs(res1[i] - res2[i]) > 1e-10) {
                        results_match = false;
                        std::cout << "Внимание: результаты не совпадают на позиции " << i 
                                 << ": " << res1[i] << " vs " << res2[i] << "\n";
                        break;
                    }
                }
                if (results_match) {
                    std::cout << "Проверка пройдена: результаты функций совпадают\n\n";
                }
            }
        }

        // Вычисление и вывод статистики для серии
        Statistics stats1(times1);
        Statistics stats2(times2);

        print_statistics("ilia_mult", stats1);
        print_statistics("Artem_matrix_mult", stats2);

        // Сохраняем отношение средних времен для этой серии
        double ratio = stats2.mean / stats1.mean;
        series_ratios.push_back(ratio);
        std::cout << "Отношение времени выполнения Artem_matrix_mult/ilia_mult: " 
                  << std::fixed << std::setprecision(2) << ratio << "\n\n";
    }

    // Вычисление и вывод общей статистики
    std::cout << "\n=== ОБЩАЯ СТАТИСТИКА ПО ВСЕМ СЕРИЯМ ===\n";
    Statistics total_stats1(all_times1);
    Statistics total_stats2(all_times2);
    Statistics ratios_stats(series_ratios);

    std::cout << "\nСтатистика отношения времен выполнения (Artem_matrix_mult/ilia_mult):\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Среднее отношение: " << ratios_stats.mean << "\n";
    std::cout << "Медиана отношений: " << ratios_stats.median << "\n";
    std::cout << "Стд. отклонение отношений: " << ratios_stats.stddev << "\n";
    std::cout << "Мин. отношение: " << ratios_stats.min << "\n";
    std::cout << "Макс. отношение: " << ratios_stats.max << "\n";

    std::cout << "\nОбщая статистика по всем запускам:\n";
    print_statistics("ilia_mult (все запуски)", total_stats1);
    print_statistics("Artem_matrix_mult (все запуски)", total_stats2);

    // Освобождаем память
    delete[] a;
    delete[] b;
    delete[] a_copy;
    delete[] b_copy;
    delete[] res1;
    delete[] res2;

    return 0;
}
