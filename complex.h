
typedef struct complex_t {
    double re;
    double im;
} complex;

complex conv_from_polar(double r, double radians);
complex add(complex left, complex right);
complex multiply(complex left, complex right);