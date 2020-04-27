import math
import numpy as np
from numpy.linalg import solve
from scipy.stats import f, t
from functools import partial
from random import randint
from prettytable import PrettyTable

while True:
    def cohren(f1, f2, q=0.05):
        q1 = q / f1
        fisher_value = f.ppf(q=1 - q1, dfn=f2, dfd=(f1 - 1) * f2)
        return fisher_value / (fisher_value + f1 - 1)


    fisher = partial(f.ppf, q=1 - 0.05)
    student = partial(t.ppf, q=1 - 0.025)

    X1max,X1min,X2max,X2min,X3max,X3min = 50,10,65,25,65,50

    Xmax_average = (X1max + X2max + X3max) / 3  # Xcp(max)
    Xmin_average = (X1min + X2min + X3min) / 3  # Xcp(min)

    y_max,y_min = 235,212

    # матриця ПФЕ
    x0_factor = [1, 1, 1, 1, 1, 1, 1, 1]
    x1_factor = [-1, -1, 1, 1, -1, -1, 1, 1]
    x2_factor = [-1, 1, -1, 1, -1, 1, -1, 1]
    x3_factor = [-1, 1, 1, -1, 1, -1, -1, 1]
    x1x2_factor = [a * b for a, b in zip(x1_factor, x2_factor)]
    x1x3_factor = [a * b for a, b in zip(x1_factor, x3_factor)]
    x2x3_factor = [a * b for a, b in zip(x2_factor, x3_factor)]
    x1x2x3_factor = [a * b * c for a, b, c in zip(x1_factor, x2_factor, x3_factor)]

    m = 3  # кількість повторень кожної комбінації

    y1, y2, y3 = [], [], []
    for i in range(0, 8):
        y1.append(randint(y_min, y_max))  # заповнення у
        y2.append(randint(y_min, y_max))
        y3.append(randint(y_min, y_max))

    Y_row1 = [y1[0], y2[0], y3[0]]  # перший рядок у
    Y_row2 = [y1[1], y2[1], y3[1]]
    Y_row3 = [y1[2], y2[2], y3[2]]
    Y_row4 = [y1[3], y2[3], y3[3]]
    Y_row5 = [y1[4], y2[4], y3[4]]
    Y_row6 = [y1[5], y2[5], y3[5]]
    Y_row7 = [y1[6], y2[6], y3[6]]
    Y_row8 = [y1[7], y2[7], y3[7]]

    Y_average1 = np.average(Y_row1)  # середній у першого рядку
    Y_average2 = np.average(Y_row2)
    Y_average3 = np.average(Y_row3)
    Y_average4 = np.average(Y_row4)
    Y_average5 = np.average(Y_row5)
    Y_average6 = np.average(Y_row6)
    Y_average7 = np.average(Y_row7)
    Y_average8 = np.average(Y_row8)
    Y_average = [round(Y_average1, 3), round(Y_average2, 3), round(Y_average3, 3), round(Y_average4, 3),
                 # округлення до тисячних
                 round(Y_average5, 3), round(Y_average6, 3), round(Y_average7, 3), round(Y_average8, 3)]

    x0 = [1, 1, 1, 1, 1, 1, 1, 1]
    x1 = [15, 15, 45, 45, 15, 15, 45, 45]  # заміна -1 на x1(min) 1 на x1(max)
    x2 = [-25, 10, -25, 10, -25, 10, -25, 10]  # заміна -1 на x2(min) 1 на x2(max)
    x3 = [-45, 50, 50, -45, 50, -45, -45, 50]  # заміна -1 на x3(min) 1 на x3(max)
    x1x2 = [a * b for a, b in zip(x1, x2)]
    x1x3 = [a * b for a, b in zip(x1, x3)]
    x2x3 = [a * b for a, b in zip(x2, x3)]
    x1x2x3 = [a * b * c for a, b, c in zip(x1, x2, x3)]

    list_for_solve_b = [x0_factor, x1_factor, x2_factor, x3_factor, x1x2_factor, x1x3_factor, x2x3_factor,
                        x1x2x3_factor]
    list_for_solve_a = list(zip(x0, x1, x2, x3, x1x2, x1x3, x2x3, x1x2x3))

    N = 8  # кількість повторення дослідів
    list_bi = []  # b(i)
    for k in range(N):
        S = 0
        for i in range(N):
            S += (list_for_solve_b[k][i] * Y_average[i]) / N
        list_bi.append(round(S, 5))

    Disp1,Disp2,Disp3,Disp4,Disp5,Disp6,Disp7,Disp8 = 0,0,0,0,0,0,0,0  # дисперсії

    for i in range(m):
        Disp1 += ((Y_row1[i] - np.average(Y_row1)) ** 2) / m  # ((y-ycp)^2)/3
        Disp2 += ((Y_row2[i] - np.average(Y_row2)) ** 2) / m
        Disp3 += ((Y_row3[i] - np.average(Y_row3)) ** 2) / m
        Disp4 += ((Y_row4[i] - np.average(Y_row4)) ** 2) / m
        Disp5 += ((Y_row5[i] - np.average(Y_row5)) ** 2) / m
        Disp6 += ((Y_row6[i] - np.average(Y_row6)) ** 2) / m
        Disp7 += ((Y_row7[i] - np.average(Y_row7)) ** 2) / m
        Disp8 += ((Y_row8[i] - np.average(Y_row8)) ** 2) / m
    sum_dispersion = Disp1 + Disp2 + Disp3 + Disp4 + Disp5 + Disp6 + Disp7 + Disp8  # сума дисперсій
    disp_list = [round(Disp1, 3), round(Disp2, 3), round(Disp3, 3), round(Disp4, 3), round(Disp5, 3), round(Disp6, 3),
                 round(Disp7, 3), round(Disp8, 3)]

    pt1 = PrettyTable()  # формування таблиці

    column_names1 = ["X0", "X1", "X2", "X3", "X1X2", "X1X3", "X2X3", "X1X2X3", "Y1", "Y2", "Y3", "Y", "S^2"]  # назви
    pt1.add_column(column_names1[0], x0_factor)  # запис назв
    pt1.add_column(column_names1[1], x1_factor)
    pt1.add_column(column_names1[2], x2_factor)
    pt1.add_column(column_names1[3], x3_factor)
    pt1.add_column(column_names1[4], x1x2_factor)
    pt1.add_column(column_names1[5], x1x3_factor)
    pt1.add_column(column_names1[6], x2x3_factor)
    pt1.add_column(column_names1[7], x1x2x3_factor)
    pt1.add_column(column_names1[8], y1)
    pt1.add_column(column_names1[9], y2)
    pt1.add_column(column_names1[10], y3)
    pt1.add_column(column_names1[11], Y_average)
    pt1.add_column(column_names1[12], disp_list)
    print(pt1, "\n")
    # рівняння регресії з ефектом взаємодії
    print("y = {} + {}*x1 + {}*x2 + {}*x3 + {}*x1x2 + {}*x1x3 + {}*x2x3 + {}*x1x2x3 \n".format(list_bi[0], list_bi[1],
                                                                                               list_bi[2], list_bi[3],
                                                                                               list_bi[4], list_bi[5],
                                                                                               list_bi[6], list_bi[7]))

    pt2 = PrettyTable()  # створення таблиці
    pt2.add_column(column_names1[0], x0)  # запис назв (113 рядок)
    pt2.add_column(column_names1[1], x1)
    pt2.add_column(column_names1[2], x2)
    pt2.add_column(column_names1[3], x3)
    pt2.add_column(column_names1[4], x1x2)
    pt2.add_column(column_names1[5], x1x3)
    pt2.add_column(column_names1[6], x2x3)
    pt2.add_column(column_names1[7], x1x2x3)
    pt2.add_column(column_names1[8], y1)
    pt2.add_column(column_names1[9], y2)
    pt2.add_column(column_names1[10], y3)
    pt2.add_column(column_names1[11], Y_average)
    pt2.add_column(column_names1[12], disp_list)
    print(pt2, '\n')

    list_ai = [round(i, 5) for i in solve(list_for_solve_a, Y_average)]
    print("y = {} + {}*x1 + {}*x2 + {}*x3 + {}*x1x2 + {}*x1x3 + {}*x2x3 + {}*x1x2x3".format(list_ai[0], list_ai[1],
                                                                                            list_ai[2], list_ai[3],
                                                                                            list_ai[4], list_ai[5],
                                                                                            list_ai[6], list_ai[7]))

    Gp = max(Disp1, Disp2, Disp3, Disp4, Disp5, Disp6, Disp7, Disp8) / sum_dispersion  # експерементальне
    F1 = m - 1
    N = len(y1)
    F2 = N
    Gt = cohren(F1, F2)  # теоретичне
    print("\nGp = ", Gp, " Gt = ", Gt)
    if Gp < Gt:
        print("Дисперсія однорідна!\n")

        Dispersion_B = sum_dispersion / N
        Dispersion_beta = Dispersion_B / (m * N)
        S_beta = math.sqrt(abs(Dispersion_beta))

        lst = [x0_factor, x1_factor, x2_factor, x3_factor, x1x2_factor, x1x3_factor,x2x3_factor, x1x2x3_factor]
        beta_list_test, beta_list = [], []
        def beta_calc():
            for j in range(len(x0_factor)):
                for i in range(len(x0_factor)):
                    beta_list_test.append((Y_average[i] * lst[j][i]) / N)
                beta_list.append(sum(beta_list_test))
                beta_list_test.clear()
        beta_calc()
        print(f"beta_list:{beta_list}")

        t_list = []
        def t_calc():
            for i in range(N):
                t_list.append(abs(beta_list[i]) / S_beta)
        t_calc()

        F3 = F1 * F2
        d = 0
        T = student(df=F3)
        print("t табличне = ", T)
        for i in range(len(t_list)):
            if t_list[i] < T:
                beta_list[i] = 0
                print("Гіпотеза підтверджена, beta{} = 0".format(i))
            else:
                print("Гіпотеза не підтверджена.\nbeta{} = {}".format(i, beta_list[i]))
                d += 1

        Y_counted_for_Student=[]
        def y_calc():
            for i in range(N):
                Y_counted_for_Student.append(beta_list[0] + beta_list[1] * x1[i] + beta_list[2] * x2[i] + beta_list[3] * x3[i] + beta_list[4] * x1x2[i] + beta_list[5] * x1x3[i] + beta_list[6] * x2x3[i] + beta_list[7] * x1x2x3[i])
        y_calc()

        F4 = N - d
        Dispersion_ad = 0
        for i in range(len(Y_counted_for_Student)):
            Dispersion_ad += ((Y_counted_for_Student[i] - Y_average[i]) ** 2) * m / (N - d)
        Fp = Dispersion_ad / Dispersion_beta
        Ft = fisher(dfn=F4, dfd=F3)
        if Fp > Ft:
            print("Рівняння регресії неадекватне.")
            print("--------------------------------------5 лаба -----------------------------------------")
            import math
            from _pydecimal import Decimal
            from scipy.stats import f, t
            import random
            from functools import reduce
            from itertools import compress
            import numpy as np

            x1min, x2min, x3min = -8, -6, -10
            x1max, x2max, x3max = 8, 3, 7

            x_avr_min = -3
            x_avr_max = 7

            x0_i = [(x1max + x1min) / 2, (x2max + x2min) / 2, (x3max + x3min) / 2]
            det_x_i = [(x1min - x0_i[0]), (x2min - x0_i[1]), (x3min - x0_i[2])]
            l = 1.215

            raw_naturalized_factors_table = [[x1min, x2min, x3min],
                                             [x1min, x2max, x3max],
                                             [x1max, x2min, x3max],
                                             [x1max, x2max, x3min],

                                             [x1min, x2min, x3max],
                                             [x1min, x2max, x3min],
                                             [x1max, x2min, x3min],
                                             [x1max, x2max, x3max],

                                             [-l * det_x_i[0] + x0_i[0], x0_i[1], x0_i[2]],
                                             [l * det_x_i[0] + x0_i[0], x0_i[1], x0_i[2]],
                                             [x0_i[0], -l * det_x_i[1] + x0_i[1], x0_i[2]],
                                             [x0_i[0], l * det_x_i[1] + x0_i[1], x0_i[2]],
                                             [x0_i[0], x0_i[1], -l * det_x_i[2] + x0_i[2]],
                                             [x0_i[0], x0_i[1], l * det_x_i[2] + x0_i[2]],

                                             [x0_i[0], x0_i[1], x0_i[2]]]

            raw_factors_table = [[-1, -1, -1],
                                 [-1, +1, +1],
                                 [+1, -1, +1],
                                 [+1, +1, -1],

                                 [-1, -1, +1],
                                 [-1, +1, -1],
                                 [+1, -1, -1],
                                 [+1, +1, +1],

                                 [-1.215, 0, 0],
                                 [+1.215, 0, 0],
                                 [0, -1.215, 0],
                                 [0, +1.215, 0],
                                 [0, 0, -1.215],
                                 [0, 0, +1.215],

                                 [0, 0, 0]]


            def generate_factors_table(raw_array):
                return [row + [row[0] * row[1], row[0] * row[2], row[1] * row[2], row[0] * row[1] * row[2]]
                        + list(map(lambda x: round(x ** 2, 5), row))
                        for row in raw_array]


            def x_i(i):
                try:
                    assert i <= 10
                except:
                    raise AssertionError("i must be smaller or equal 10")
                with_null_factor = list(map(lambda x: [1] + x, generate_factors_table(raw_factors_table)))
                res = [row[i] for row in with_null_factor]
                return np.array(res)


            def cochran_criteria(m, N, y_table):
                print()
                print(
                    "---------------------------------------Criterion Cohrena----------------------------------------")
                print()
                y_variations = [np.var(i) for i in y_table]
                max_y_variation = max(y_variations)
                gp = max_y_variation / sum(y_variations)
                f1 = m - 1
                f2 = N
                p = 0.95
                q = 1 - p
                gt = get_cochran_value(f1, f2, q)
                print("Gp = {}; Gt = {}; f1 = {}; f2 = {}; q = {:.2f}".format(gp, gt, f1, f2, q))
                if gp < gt:
                    print("Gp < Gt => дисперсії рівномірні ")
                    return True
                else:
                    print("Gp > Gt => дисперсії нерівномірні ")
                    return False


            def student_criteria(m, N, y_table, beta_coefficients):
                print()
                print(
                    "---------------------------------------Criterion Student`a----------------------------------------")
                print()
                average_variation = np.average(list(map(np.var, y_table)))

                variation_beta_s = average_variation / N / m
                standard_deviation_beta_s = math.sqrt(variation_beta_s)
                t_i = np.array(
                    [abs(beta_coefficients[i]) / standard_deviation_beta_s for i in range(len(beta_coefficients))])
                f3 = (m - 1) * N
                q = 0.05

                t = get_student_value(f3, q)
                importance = [True if el > t else False for el in list(t_i)]

                print("Оцінки коефіцієнтів βs: " + ", ".join(
                    list(map(lambda x: str(round(float(x), 3)), beta_coefficients))))
                print("Коефіцієнти ts:         " + ", ".join(list(map(lambda i: "{:.2f}".format(i), t_i))))
                print("f3 = {}; q = {}; tтабл = {}".format(f3, q, t))
                beta_i = ["β0", "β1", "β2", "β3", "β12", "β13", "β23", "β123", "β11", "β22", "β33"]
                importance_to_print = ["важливий" if i else "неважливий" for i in importance]
                to_print = map(lambda x: x[0] + " " + x[1], zip(beta_i, importance_to_print))
                x_i_names = list(
                    compress(["", "x1", "x2", "x3", "x12", "x13", "x23", "x123", "x1^2", "x2^2", "x3^2"], importance))
                betas_to_print = list(compress(beta_coefficients, importance))
                print(*to_print, sep="; ")
                equation = " ".join(
                    ["".join(i) for i in zip(list(map(lambda x: "{:+.2f}".format(x), betas_to_print)), x_i_names)])
                print("Рівняння регресії без незначимих членів: y = " + equation)
                return importance


            def calculate_theoretical_y(x_table, b_coefficients, importance):
                x_table = [list(compress(row, importance)) for row in x_table]
                b_coefficients = list(compress(b_coefficients, importance))
                y_vals = np.array([sum(map(lambda x, b: x * b, row, b_coefficients)) for row in x_table])
                return y_vals


            def fisher_criteria(m, N, d, naturalized_x_table, y_table, b_coefficients, importance):
                f3 = (m - 1) * N
                f4 = N - d
                q = 0.05

                theoretical_y = calculate_theoretical_y(naturalized_x_table, b_coefficients, importance)
                theoretical_values_to_print = list(zip(map(lambda x: "x1 = {0[1]}, x2 = {0[2]}, x3 = {0[3]}".format(x),
                                                           naturalized_x_table), theoretical_y))

                y_averages = np.array(list(map(np.average, y_table)))
                s_ad = m / (N - d) * (sum((theoretical_y - y_averages) ** 2))
                y_variations = np.array(list(map(np.var, y_table)))
                s_v = np.average(y_variations)
                f_p = float(s_ad / s_v)
                f_t = get_fisher_value(f3, f4, q)

                print()
                print(
                    "---------------------------------------Criterion Fishera----------------------------------------")
                print()
                print("Теоретичні значення y для різних комбінацій факторів:")
                print("\n".join(["{arr[0]}: y = {arr[1]}".format(arr=el) for el in theoretical_values_to_print]))
                print("Fp = {}, Ft = {}".format(f_p, f_t))
                print("Fp < Ft => модель адекватна" if f_p < f_t else "Fp > Ft => модель неадекватна")
                return True if f_p < f_t else False


            def m_ij(*arrays):
                return np.average(reduce(lambda accum, el: accum * el, arrays))


            def get_cochran_value(f1, f2, q):
                partResult1 = q / f2  # (f2 - 1)
                params = [partResult1, f1, (f2 - 1) * f1]
                fisher = f.isf(*params)
                result = fisher / (fisher + (f2 - 1))
                return Decimal(result).quantize(Decimal('.0001'))


            def get_student_value(f3, q):
                return Decimal(abs(t.ppf(q / 2, f3))).quantize(Decimal('.0001'))


            def get_fisher_value(f3, f4, q):
                return Decimal(abs(f.isf(q, f4, f3))).quantize(Decimal('.0001'))


            print()
            print("---------------------------------------start----------------------------------------")
            print()

            factors_table = generate_factors_table(raw_factors_table)
            for row in factors_table:
                print(row)
            naturalized_factors_table = generate_factors_table(raw_naturalized_factors_table)
            with_null_factor = list(map(lambda x: [1] + x, naturalized_factors_table))

            m = 32
            N = 15
            ymin = 196
            ymax = 209
            y_arr = [[random.randint(ymin, ymax) for _ in range(m)] for _ in range(N)]
            while not cochran_criteria(m, N, y_arr):
                m += 1
                y_arr = [[random.randint(ymin, ymax) for _ in range(m)] for _ in range(N)]

            y_i = np.array([np.average(row) for row in y_arr])

            coefficients = [[m_ij(x_i(column) * x_i(row)) for column in range(11)] for row in range(11)]

            free_values = [m_ij(y_i, x_i(i)) for i in range(11)]

            beta_coefficients = np.linalg.solve(coefficients, free_values)
            print(list(map(int, beta_coefficients)))

            importance = student_criteria(m, N, y_arr, beta_coefficients)
            d = len(list(filter(None, importance)))
            fisher_criteria(m, N, d, naturalized_factors_table, y_arr, beta_coefficients, importance)
            break
        else:
            print("Рівняння регресії адекватне!")
            break

    else:
        print("Дисперсія неоднорідна. Спробуйте ще раз.")
        m += 1