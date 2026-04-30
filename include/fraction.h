#pragma once
#include <iostream>
#include <string>
#include <gmp.h>
#include <stdexcept>

class Fraction {
private:
    mpz_t num; // Числитель
    mpz_t den; // Знаменатель
    void simplify();
public:
    // конструкторы и деструктор
    Fraction(long long n = 0, long long d = 1);
    Fraction(const Fraction& other);      // Копирование
    Fraction(Fraction&& other) noexcept;  // Перемещение
    ~Fraction();
    // операторы присваивания
    Fraction& operator=(const Fraction& other);
    Fraction& operator=(Fraction&& other) noexcept;
    // проверк
    bool isZero() const;
    Fraction operator+(const Fraction& other) const;
    Fraction operator-(const Fraction& other) const;
    Fraction operator*(const Fraction& other) const;
    Fraction operator/(const Fraction& other) const;
    // составное присваивание (изменяют текущий обьект)
    Fraction& operator+=(const Fraction& other);
    Fraction& operator-=(const Fraction& other);
    Fraction& operator*=(const Fraction& other);
    Fraction& operator/=(const Fraction& other);
    // операторы сравнения (нужны для метода гаусса)
    bool operator==(const Fraction& other) const;
    bool operator!=(const Fraction& other) const;
    bool operator<(const Fraction& other) const;
    bool operator>(const Fraction& other) const;
    // разное
    std::string toDecimal(int precision = 10) const;
    friend std::ostream& operator<<(std::ostream& os, const Fraction& frac);
};
// функция ввода с клавиатуры
Fraction inputFraction(const std::string& fractionName);
