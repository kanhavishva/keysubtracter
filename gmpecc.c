#include <gmp.h>
#include "gmpecc.h"


void Point_Doubling(const Point *P, Point *R)
{
    mpz_t slope, temp;
    mpz_init(temp);
    mpz_init(slope);
    if (mpz_cmp_ui(P->y, 0) != 0) {
        mpz_mul_ui(temp, P->y, 2);
        mpz_invert(temp, temp, EC.p);
        mpz_mul(slope, P->x, P->x);
        mpz_mul_ui(slope, slope, 3);
        mpz_mul(slope, slope, temp);
        mpz_mod(slope, slope, EC.p);
        mpz_mul(R->x, slope, slope);
        mpz_sub(R->x, R->x, P->x);
        mpz_sub(R->x, R->x, P->x);
        mpz_mod(R->x, R->x, EC.p);
        mpz_sub(temp, P->x, R->x);
        mpz_mul(R->y, slope, temp);
        mpz_sub(R->y, R->y, P->y);
        mpz_mod(R->y, R->y, EC.p);
    } else {
        mpz_set_ui(R->x, 0);
        mpz_set_ui(R->y, 0);
    }
    mpz_clear(temp);
    mpz_clear(slope);
}

void Point_Addition(Point *P, Point *Q, Point *R)
{
    mpz_t T, S;
    mpz_init(T);
    mpz_init(S);
    mpz_mod(Q->x, Q->x, EC.p);
    mpz_mod(Q->y, Q->y, EC.p);
    mpz_mod(P->x, P->x, EC.p);
    mpz_mod(P->y, P->y, EC.p);
    if (mpz_cmp_ui(P->x, 0) == 0 && mpz_cmp_ui(P->y, 0) == 0) {
        mpz_set(R->x, Q->x);
        mpz_set(R->y, Q->y);
    } else {
        if (mpz_cmp_ui(Q->y, 0) != 0) {
            mpz_sub(T, EC.p, Q->y);
            mpz_mod(T, T, EC.p);
        } else {
            mpz_set_ui(T, 0);
        }
        if (mpz_cmp(P->y, T) == 0 && mpz_cmp(P->x, Q->x) == 0) {
            mpz_set_ui(R->x, 0);
            mpz_set_ui(R->y, 0);
        } else {
            if (mpz_cmp(P->x, Q->x) == 0 && mpz_cmp(P->y, Q->y) == 0) {
                Point_Doubling(P, R);
            } else {
                mpz_set_ui(S, 0);
                mpz_sub(T, P->x, Q->x);
                mpz_mod(T, T, EC.p);
                mpz_invert(T, T, EC.p);
                mpz_sub(S, P->y, Q->y);
                mpz_mul(S, S, T);
                mpz_mod(S, S, EC.p);
                mpz_mul(R->x, S, S);
                mpz_sub(R->x, R->x, P->x);
                mpz_sub(R->x, R->x, Q->x);
                mpz_mod(R->x, R->x, EC.p);
                mpz_sub(T, P->x, R->x);
                mpz_mul(R->y, S, T);
                mpz_sub(R->y, R->y, P->y);
                mpz_mod(R->y, R->y, EC.p);
            }
        }
    }
    mpz_clear(T);
    mpz_clear(S);
}

void Scalar_Multiplication(const Point *P, Point *R, const mpz_t *m)
{
    Point T, Q;
    long no_of_bits, i;
    no_of_bits = mpz_sizeinbase(*m, 2);
    mpz_init_set_ui(Q.x, 0);
    mpz_init_set_ui(Q.y, 0);
    mpz_init_set_ui(T.x, 0);
    mpz_init_set_ui(T.y, 0);
    mpz_set_ui(R->x, 0);
    mpz_set_ui(R->y, 0);
    if (mpz_cmp_ui(*m, 0) != 0) {
        mpz_set(Q.x, P->x);
        mpz_set(Q.y, P->y);
        for (i = 0; i < no_of_bits; i++) {
            if (mpz_tstbit(*m, i)) {
                mpz_set(T.x, R->x);
                mpz_set(T.y, R->y);
                mpz_set(Q.x, DoublingG[i].x);
                mpz_set(Q.y, DoublingG[i].y);
                Point_Addition(&T, &Q, R);
            }
        }
    }
    mpz_clear(T.x);
    mpz_clear(T.y);
    mpz_clear(Q.x);
    mpz_clear(Q.y);
}

void Point_Negation(const Point *A, Point *S)
{
    Point Q;
    mpz_t T;
    mpz_init(T);
    mpz_init(Q.x);
    mpz_init(Q.y);
    mpz_set(Q.x, A->x);
    mpz_set(Q.y, A->y);
    mpz_sub(T, EC.p, Q.y);
    mpz_set(S->x, Q.x);
    mpz_set(S->y, T);
    mpz_clear(T);
    mpz_clear(Q.x);
    mpz_clear(Q.y);
}

/*
        Precalculate G Doublings for Scalar_Multiplication
*/
void init_doublingG(const Point *P)
{
    int i = 0;
    mpz_init(DoublingG[i].x);
    mpz_init(DoublingG[i].y);
    mpz_set(DoublingG[i].x, P->x);
    mpz_set(DoublingG[i].y, P->y);
    i = 1;
    while (i < 256) {
        mpz_init(DoublingG[i].x);
        mpz_init(DoublingG[i].y);
        Point_Doubling(&DoublingG[i - 1], &DoublingG[i]);
        mpz_mod(DoublingG[i].x, DoublingG[i].x, EC.p);
        mpz_mod(DoublingG[i].y, DoublingG[i].y, EC.p);
        i++;
    }
}
