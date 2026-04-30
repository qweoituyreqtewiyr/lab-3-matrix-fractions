#include <iostream>
#include <vector>
#include "fraction.h"
#include "matrix.h"

// Вспомогательная функция для вывода вектора в строку
void printVector(const std::vector<Fraction>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << (i < vec.size() - 1 ? "\t" : "");
    }
    std::cout << "\n";
}

int main() {
    try {
        // Подготовка общих данных для тестов умножения
        Matrix A(2, 2);
        A(0, 0) = Fraction(1); A(0, 1) = Fraction(2);
        A(1, 0) = Fraction(3); A(1, 1) = Fraction(4);

        Matrix B(2, 2);
        B(0, 0) = Fraction(5); B(0, 1) = Fraction(6);
        B(1, 0) = Fraction(7); B(1, 1) = Fraction(8);

        std::vector<Fraction> vec = {Fraction(2), Fraction(3)};

        // Подготовка данных для СЛАУ и обратных матриц (система 3x3)
        Matrix Sys(3, 3);
        Sys(0,0)=Fraction(2);  Sys(0,1)=Fraction(1);  Sys(0,2)=Fraction(-1);
        Sys(1,0)=Fraction(-3); Sys(1,1)=Fraction(-1); Sys(1,2)=Fraction(2);
        Sys(2,0)=Fraction(-2); Sys(2,1)=Fraction(1);  Sys(2,2)=Fraction(2);

        std::vector<Fraction> b = {Fraction(8), Fraction(-11), Fraction(-3)};

        // ---------------------------------------------------------
        std::cout << "Тест 1\n";
        std::cout << "Условие:\nДана матрица A:\n";
        A.print();
        std::cout << "Дан вектор v:\n";
        printVector(vec);
        std::cout << "Метод: умножение матрицы на вектор\n";
        std::cout << "Ответ:\n";
        printVector(A * vec);
        std::cout << "\n";

        // ---------------------------------------------------------
        std::cout << "Тест 2\n";
        std::cout << "Условие:\nДана матрица A:\n"; A.print();
        std::cout << "Дана матрица B:\n"; B.print();
        std::cout << "Метод: классический алгоритм умножения матриц\n";
        std::cout << "Ответ:\n";
        (A * B).print();
        std::cout << "\n";

        // ---------------------------------------------------------
        std::cout << "Тест 3\n";
        std::cout << "Условие:\nДана матрица A:\n"; A.print();
        std::cout << "Дана матрица B:\n"; B.print();
        std::cout << "Метод: блочный алгоритм умножения матриц\n";
        std::cout << "Ответ:\n";
        A.multiplyBlock(B, 2).print(); // Размер блока 2 для теста
        std::cout << "\n";

        // ---------------------------------------------------------
        std::cout << "Тест 4\n";
        std::cout << "Условие:\nДана матрица A:\n"; A.print();
        std::cout << "Дана матрица B:\n"; B.print();
        std::cout << "Метод: алгоритм Штрассена для умножения матриц\n";
        std::cout << "Ответ:\n";
        A.multiplyStrassen(B).print();
        std::cout << "\n";

        // ---------------------------------------------------------
        std::cout << "Тест 5\n";
        std::cout << "Условие:\nДана матрица коэффициентов:\n"; Sys.print();
        std::cout << "Вектор свободных членов:\n"; printVector(b);
        std::cout << "Метод: решение систем линейных алгебраических уравнений методом Гаусса\n";
        std::cout << "Ответ:\n";
        printVector(Sys.solveGauss(b, false));
        std::cout << "\n";

        // ---------------------------------------------------------
        std::cout << "Тест 6\n";
        std::cout << "Условие:\nДана матрица коэффициентов:\n"; Sys.print();
        std::cout << "Вектор свободных членов:\n"; printVector(b);
        std::cout << "Метод: решение систем линейных алгебраических уравнений методом Жордана-Гаусса\n";
        std::cout << "Ответ:\n";
        printVector(Sys.solveJordanGauss(b, false));
        std::cout << "\n";

        // ---------------------------------------------------------
        std::cout << "Тест 7\n";
        std::cout << "Условие:\nДана квадратная матрица:\n"; Sys.print();
        std::cout << "Метод: нахождение обратной матрицы методом Гаусса\n";
        std::cout << "Ответ:\n";
        Sys.inverseGauss(false).print();
        std::cout << "\n";

        // ---------------------------------------------------------
        std::cout << "Тест 8\n";
        std::cout << "Условие:\nДана квадратная матрица:\n"; Sys.print();
        std::cout << "Метод: нахождение обратной матрицы методом Жордана-Гаусса\n";
        std::cout << "Ответ:\n";
        Sys.inverseJordanGauss(false).print();
        std::cout << "\n";

    } catch (const std::exception& e) {
        std::cerr << "КРИТИЧЕСКАЯ ОШИБКА: " << e.what() << "\n";
    }

    return 0;
}
