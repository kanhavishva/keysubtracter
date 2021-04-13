#ifndef GMPECC_H
#define GMPECC_H

typedef struct Point {
    mpz_t x;
    mpz_t y;
} Point;

typedef struct Elliptic_Curve {
    mpz_t p;
    mpz_t n;
} Elliptic_Curve;

void Point_Doubling(const Point *P, Point *R);
void Point_Addition(Point *P, Point *Q, Point *R);
void Scalar_Multiplication(const Point *P, Point *R, const mpz_t *m);
void Point_Negation(const Point *A, Point *S);
void init_doublingG(const Point *P);

Elliptic_Curve EC;
Point G;
Point DoublingG[256];


#endif
