#include <iostream>
#include <cstdlib>
#include <sstream>

namespace numbers {

    unsigned long long gcd(unsigned long long a, unsigned long long b) {
        while (b != 0) {
            a %= b;
            std::swap(a, b);
        }
        return a;
    }

//    template<typename T1, typename T2>
//    [[maybe_unused]] unsigned long long gcd(const T1& a, const T2& b) = delete;

    unsigned long long lcm(unsigned long long a, unsigned long long b) {
        unsigned long long g = gcd(a, b);
        a /= g;
        return a * b;
    }

    class RatioSignIn {
    public:
        long long numerator;
        unsigned long long denominator;

        // если ты ещё раз захочешь написать тут explicit,
        // то не делай этого! Будут проблемы с "= 0" в FSS
        RatioSignIn(long long n = 0, unsigned long long d = 1) {
            reduce(n, d);
            numerator = n, denominator = d;
        }

        bool operator==(int b) const {
            return (numerator == b && denominator == 1);
        }

        RatioSignIn& operator+=(const RatioSignIn other) {
            long long new_numerator;
            unsigned long long new_denominator;
            new_denominator = lcm(denominator, other.denominator);
            new_numerator = numerator * (long long) (new_denominator / denominator) + other.numerator * (long long) (new_denominator / other.denominator);
            reduce(new_numerator, new_denominator);
            numerator = new_numerator, denominator = new_denominator;
            return *this;
        }

        RatioSignIn& operator-=(const RatioSignIn other) {
            *this += -other;
            return *this;
        }

        RatioSignIn& operator*=(const RatioSignIn other) {
            long long n1 = numerator, n2 = other.numerator;
            unsigned long long d1 = denominator, d2 = other.denominator;
            reduce(n1, d2);
            reduce(n2, d1);
            numerator = n1 * n2, denominator = d1 * d2;
            return *this;
        }

        RatioSignIn& operator/=(const RatioSignIn other) {
            *this *= ~other;
            return *this;
        }

        RatioSignIn operator-() const {
            return {-numerator, denominator};
        }

        RatioSignIn operator~() const {
            if (numerator == 0) {
                return {0, 1};
            } else {
                long long new_numerator;
                unsigned long long new_denominator;
                if (numerator > 0) {
                    new_numerator = (long long) denominator;
                    new_denominator = (unsigned long long) numerator;
                } else {
                    new_numerator = (long long) -denominator;
                    new_denominator = (unsigned long long) -numerator;
                }
                return {new_numerator, new_denominator};
            }
        }

        static void reduce(long long& n, unsigned long long& d) {
            unsigned long long g = gcd((unsigned long long) std::abs(n), d);
            n /= (long long) g, d /= g;
        }

        void print() const {
            std::cout << numerator;
            if (denominator != 1) {
                std::cout << "/" << denominator;
            }
        }
    };

    class RatioSignOut {
        bool sign;
        unsigned long long numerator, denominator;

        explicit RatioSignOut(unsigned long long n = 0, unsigned long long d = 1, bool s = false) : sign(s) {
            reduce(n, d);
            numerator = n, denominator = d;
        }

        bool operator==(int b) const {
            return (numerator == b && denominator == 1);
        }

        static void reduce(unsigned long long& n, unsigned long long& d) {
            unsigned long long g = gcd(n, d);
            n /= g, d /= g;
        }
    };

    class Residue {
    public:
        static unsigned long long modulus;
        unsigned long long remainder;

        Residue(unsigned long long n = 0) : remainder(n) {}

        bool operator==(unsigned long long b) const {
            return (remainder == b % modulus);
        }

        Residue& operator+=(const Residue& other) {
            remainder += other.remainder;
            if (remainder >= modulus) {
                remainder -= modulus;
            }
            return *this;
        }

        Residue& operator-=(const Residue& other) {
            *this += -other;
            return *this;
        }

        Residue& operator*=(const Residue& other) {
            remainder = (remainder * other.remainder) % modulus;
            return *this;
        }

        Residue& operator/=(const Residue& other) {
            *this *= ~other;
            return *this;
        }

        Residue operator-() const {
            if (remainder == 0) {
                return 0;
            } else {
                return modulus - remainder;
            }
        }

        Residue operator~() const {
            if ( gcd(remainder, modulus) != 1 ) {
                return 0;
            } else {
                unsigned long long k = 1;
                while (remainder * k % modulus != 1) {
                    ++k;
                }
                return k;
            }
        }

        void print() const {
            std::cout << remainder;
        }
    };

    typedef RatioSignIn Ratio;

    Ratio operator+(const Ratio& r1, const Ratio& r2) {
        Ratio r = r1;
        r += r2;
        return r;
    }

    Ratio operator-(const Ratio& r1, const Ratio& r2) {
        Ratio r = r1;
        r -= r2;
        return r;
    }

    Ratio operator*(const Ratio& r1, const Ratio& r2) {
        Ratio r = r1;
        r *= r2;
        return r;
    }

    Ratio operator/(const Ratio& r1, const Ratio& r2) {
        Ratio r = r1;
        r /= r2;
        return r;
    }

    class RatioVector {
    public:
        unsigned length;
        Ratio* memory;

        explicit RatioVector(unsigned l = 0) : length(l) {
            memory = new Ratio[length];
        }

        RatioVector(const RatioVector& other) : length(other.length) {
            memory = new Ratio[length];
            for (unsigned i = 0; i < length; ++i) {
                memory[i] = other[i];
            }
        }

        RatioVector& operator=(const RatioVector& other) {
            if (this != &other) {
                delete[] memory;
                length = other.length;
                memory = new Ratio[length];
                for (unsigned i = 0; i < length; ++i) {
                    memory[i] = other[i];
                }
            }
            return *this;
        }

        ~RatioVector() {
            delete[] memory;
        }

        Ratio& operator[](unsigned k) const {
            return memory[k];
        }

        const RatioVector& operator-=(const RatioVector& other) const {
            for (unsigned i = 0; i < length; ++i) {
                memory[i] -= other[i];
            }
            return *this;
        }

        void simplify() const {
            unsigned k = 0;
            while (k < length && memory[k] == 0) {
                ++k;
            }

            if (k < length) {
                Ratio r = memory[k];
                for (unsigned i = k + 1; i < length; ++i) {
                    r.numerator = (long long) gcd((unsigned long long) std::abs(r.numerator), (unsigned long long) std::abs(memory[i].numerator));
                    r.denominator = gcd(r.denominator, memory[i].denominator);
                }
                for (unsigned i = k; i < length; ++i) {
                    memory[i] /= r;
                }
            }
        }

        void print() const {
            std::cout << "(";
            memory[0].print();
            for (unsigned i = 1; i < length; ++i) {
                std::cout << ", ";
                memory[i].print();
            }
            std::cout << ")" << std::endl;
        }
    };

    RatioVector operator*(Ratio coefficient, const RatioVector& vector) {
        RatioVector res = vector;
        for (unsigned i = 0; i < res.length; ++i) {
            res.memory[i] *= coefficient;
        }
        return res;
    }

    class RatioMatrix {
    public:
        const unsigned height, width;
        RatioVector* memory;

        explicit RatioMatrix(unsigned h = 0, unsigned w = 0) : height(h), width(w) {
            memory = new RatioVector[height];
            for (unsigned i = 0; i < height; ++i) {
                memory[i] = RatioVector(width);
            }
        }

        RatioMatrix& operator=(const RatioMatrix& other) {
            if (this != &other && other.height == height && other.width == width) {
                for (unsigned i = 0; i < height; ++i) {
                    memory[i] = other.memory[i];
                }
            }
            return *this;
        }


        RatioVector& operator[](unsigned k) const {
            return memory[k];
        }

        ~RatioMatrix() {
            delete[] memory;
        }

        void print() const {
            for (unsigned i = 0; i < height; ++i) {
                memory[i].print();
            }
        }
    };

    class SLAE {
    public:
        const unsigned equations, variables;
        unsigned rank = 0;
        RatioMatrix matrix;
        RatioVector result;
        unsigned* bf;

        explicit SLAE(unsigned e, unsigned v) : equations(e), variables(v), matrix(e, v), result(e) {
            bf = new unsigned[variables];
        }

        SLAE& operator=(const SLAE& other) {
            if (this != &other && other.equations == equations && other.variables == variables) {
                matrix = other.matrix;
                result = other.result;
            }
            return *this;
        }

        ~SLAE() {
            delete[] bf;
        }

        void input() const {
            std::string str;
            char c;
            Ratio* p;

            for (unsigned i = 0; i < equations; ++i) {
                for (unsigned j = 0; j < variables + 1; ++j) {
                    std::stringstream container;
                    std::cin >> str;
                    container << str;

                    if (j < variables) {
                        p = &matrix[i][j];
                    } else {
                        p = &result[i];
                    }
                    container >> (*p).numerator;
                    if ( str.find('/') != std::string::npos ) {
                        container >> c;
                        container >> (*p).denominator;
                    } else {
                        (*p).denominator = 1;
                    }
                }
            }
        }

        void diag() {
            Ratio coefficient;
            unsigned h;

            for (unsigned i = 0; i < variables; ++i) {
                h = rank;
                while (h < equations && matrix[h][i] == 0) {
                    ++h;
                }

                if (h == equations) {
                    bf[variables - (i - rank) - 1] = i;
                } else {
                    if (h > rank) {
                        std::swap(matrix[h], matrix[rank]);
                        std::swap(result[h], result[rank]);
                    }
                    for (unsigned j = 0; j < equations; ++j) {
                        if (j != rank) {
                            coefficient = matrix[j][i] / matrix[rank][i];
                            matrix[j] -= coefficient * matrix[rank];
                            result[j] -= coefficient * result[rank];
                        }
                    }

                    bf[rank] = i;
                    ++rank;
                }
            }
        }

        bool solvable() const {
            for (unsigned i = rank; i < equations; ++i) {
                std::cout << "lol: ";
                result[i].print();
                std::cout << std::endl;
                if (result[i] != 0) {
                    return false;
                }
            }
            return true;
        }

        RatioMatrix FSS() const {
            RatioMatrix s(variables - rank + 1, variables);

            for (unsigned j = 0; j < rank; ++j) {
                s[0][bf[j]] = result[j] / matrix[j][bf[j]];
            }
            for (unsigned j = rank; j < variables; ++j) {
                s[0][bf[j]] = 0;
            }

            for (unsigned i = 1; i < 1 + variables - rank; ++i) {
                for (unsigned j = 0; j < rank; ++j) {
                    s[i][bf[j]] = -matrix[j][bf[rank + i - 1]] / matrix[j][bf[j]];
                }
                for (unsigned j = rank; j < variables; ++j) {
                    s[i][bf[j]] = 0;
                }
                s[i][bf[rank + i - 1]] = 1;
                s[i].simplify();
            }

            /*
            for (unsigned i = 0; i < rank; ++i) {
                s[0][bf[i]] = vector[i] / matrix[i][bf[i]];
                for (unsigned j = 1; j < variables - rank + 1; ++j) {
                    s[j][bf[i]] = -matrix[i][bf[rank + j - 1]] / matrix[i][bf[i]];
                }
            }
            for (unsigned i = rank; i < variables; ++i) {
                s[0][bf[i]] = 0;
                for (unsigned j = 1; j < variables - rank + 1; ++j) {
                    s[j][bf[i]] = 0;
                }
                s[i - rank + 1][bf[i]] = 1;
            }
            */

            return s;
        }

        [[maybe_unused]] void print() const {
            std::cout << "matrix:" << std::endl;
            matrix.print();
            std::cout << "vector:" << std::endl;
            result.print();
            std::cout << std::endl;
        }
    };
}

int main() {
    /*
     * Hello! I can help you with solving SLAE!
     * Please, give me a SLAE like that:
     *
     * eqs vars
     * a_{1, 1} a_{1, 2} ... a_{1, vars} b_1
     * a_{2, 1} a_{2, 2} ... a_{2, vars} b_2
     * ...
     * a_{eqs, 1} a_{eqs, 2} ... a_{eqs, vars} b_eqs
     */

    unsigned equations = 0, variables = 0;
    std::cin >> equations >> variables;

    numbers::SLAE system(equations, variables);
    system.input();
    system.diag();

    if ( system.solvable() ) {
        std::cout << "YES" << std::endl;
        numbers::RatioMatrix gs = system.FSS();
        gs[0].print();
        for (unsigned i = 1; i < system.variables - system.rank + 1; ++i) {
            std::cout << "+ c" << i << " * ";
            gs[i].print();
        }
    } else {
        std::cout << "NO" << std::endl;
    }

    return 0;
}