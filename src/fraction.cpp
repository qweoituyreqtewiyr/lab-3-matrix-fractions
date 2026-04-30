#include "fraction.h"
#include <cstdlib>
void Fraction::simplify() {
    if (mpz_sgn(den) < 0) {
        mpz_neg(num, num);
        mpz_neg(den, den);
    }
    if (mpz_cmp_ui(num, 0) == 0) {
        mpz_set_ui(den, 1);
        return;
    }
    mpz_t gcd;
    mpz_init(gcd);
    mpz_gcd(gcd, num, den);
    mpz_divexact(num, num, gcd);
    mpz_divexact(den, den, gcd);
    mpz_clear(gcd);
}
Fraction::Fraction(long long n, long long d) {
    if (d == 0) throw std::invalid_argument("Знаменатель не может быть нулем!");
    mpz_init_set_si(num, n);
    mpz_init_set_si(den, d);
    simplify();
}
Fraction::Fraction(const Fraction& other) {
    mpz_init_set(num, other.num);
    mpz_init_set(den, other.den);
}
//перемещение указателей GMP
Fraction::Fraction(Fraction&& other) noexcept {
    mpz_init(num);
    mpz_init(den);
    mpz_swap(num, other.num);
    mpz_swap(den, other.den);
}
Fraction& Fraction::operator=(const Fraction& other) {
    if (this != &other) {
        mpz_set(num, other.num);
        mpz_set(den, other.den);
    }
    return *this;
}
Fraction& Fraction::operator=(Fraction&& other) noexcept {
    if (this != &other) {
        mpz_swap(num, other.num);
        mpz_swap(den, other.den);
    }
    return *this;
}
Fraction::~Fraction() {
    mpz_clear(num);
    mpz_clear(den);
}
bool Fraction::isZero() const {
    return mpz_cmp_ui(num, 0) == 0;
}
Fraction Fraction::operator+(const Fraction& other) const {
    Fraction result(*this);
    result += other;
    return result;
}
Fraction Fraction::operator-(const Fraction& other) const {
    Fraction result(*this);
    result -= other;
    return result;
}
Fraction Fraction::operator*(const Fraction& other) const {
    Fraction result(*this);
    result *= other;
    return result;
}
Fraction Fraction::operator/(const Fraction& other) const {
    Fraction result(*this);
    result /= other;
    return result;
}
Fraction& Fraction::operator+=(const Fraction& other) {
    mpz_t temp1, temp2;
    mpz_init(temp1); mpz_init(temp2);
    mpz_mul(temp1, num, other.den);
    mpz_mul(temp2, other.num, den);
    mpz_add(num, temp1, temp2);
    mpz_mul(den, den, other.den);
    simplify();
    mpz_clear(temp1); mpz_clear(temp2);
    return *this;
}
Fraction& Fraction::operator-=(const Fraction& other) {
    mpz_t temp1, temp2;
    mpz_init(temp1); mpz_init(temp2);
    mpz_mul(temp1, num, other.den);
    mpz_mul(temp2, other.num, den);
    mpz_sub(num, temp1, temp2);
    mpz_mul(den, den, other.den);
    simplify();
    mpz_clear(temp1); mpz_clear(temp2);
    return *this;
}
Fraction& Fraction::operator*=(const Fraction& other) {
    mpz_mul(num, num, other.num);
    mpz_mul(den, den, other.den);
    simplify();
    return *this;
}
Fraction& Fraction::operator/=(const Fraction& other) {
    if (other.isZero()) throw std::runtime_error("Деление на ноль!");
    mpz_mul(num, num, other.den);
    mpz_mul(den, den, other.num);
    simplify();
    return *this;
}
bool Fraction::operator==(const Fraction& other) const {
    return mpz_cmp(num, other.num) == 0 && mpz_cmp(den, other.den) == 0;
}
bool Fraction::operator!=(const Fraction& other) const { return !(*this == other); }
bool Fraction::operator<(const Fraction& other) const {
    mpz_t left, right;
    mpz_init(left); mpz_init(right);
    mpz_mul(left, num, other.den);
    mpz_mul(right, other.num, den);
    bool result = mpz_cmp(left, right) < 0;
    mpz_clear(left); mpz_clear(right);
    return result;
}
bool Fraction::operator>(const Fraction& other) const { return other < *this; }
std::string Fraction::toDecimal(int precision) const {
    std::string result;
    mpz_t q, r;
    mpz_init(q); mpz_init(r);
    mpz_tdiv_qr(q, r, num, den);
    if (mpz_sgn(num) < 0 && mpz_cmp_ui(q, 0) == 0) result += "-";
    char* q_str = mpz_get_str(NULL, 10, q);
    result += q_str;
    free(q_str);
    if (precision > 0 && mpz_cmp_ui(r, 0) != 0) {
        result += ".";
        mpz_abs(r, r);
        for (int i = 0; i < precision; ++i) {
            if (mpz_cmp_ui(r, 0) == 0) break;
            mpz_mul_ui(r, r, 10);
            mpz_tdiv_qr(q, r, r, den);
            char* digit = mpz_get_str(NULL, 10, q);
            result += digit;
            free(digit);
        }
    }
    mpz_clear(q); mpz_clear(r);
    return result;
}
std::ostream& operator<<(std::ostream& os, const Fraction& frac) {
    char* num_str = mpz_get_str(NULL, 10, frac.num);
    char* den_str = mpz_get_str(NULL, 10, frac.den);
    os << num_str;
    if (mpz_cmp_ui(frac.den, 1) != 0) os << "/" << den_str;
    free(num_str); free(den_str);
    return os;
}
