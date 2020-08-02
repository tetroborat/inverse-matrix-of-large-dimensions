#ifndef NORMAL_SYSTEM_EQUATIONS_H
#define NORMAL_SYSTEM_EQUATIONS_H
#include <QVector>

class normal_system_equations
{
public:
    normal_system_equations();
    static QVector<double> multiplication_matrix_and_number(double number, QVector<double> matrix);
    static QVector<double> projection_count(QVector<double> a, QVector<double> b);
    static double scalar_multiplication_string_matrix(QVector<double> matrix1, QVector<double> matrix2);
    static QVector<QVector<double> > multiplication_matrix(QVector<QVector<double> > matrix1, QVector<QVector<double> > matrix2);
    static QVector<QVector<double> > transposition_matrix(QVector<QVector<double> > matrix);
    static QVector<double> difference_string_matrix(QVector<double> matrix1, QVector<double> matrix2);
    static QVector<QVector<double> > ortogonalization_and_ortonormal_matrix(QVector<QVector<double> > matrix);
    static QVector<double> count_b(int i, QVector<QVector<double> > A);
    static QVector<QVector<double> > reverse_matrix_3x3_count(QVector<QVector<double> > matrix);
    static double determenant_matrix_3x3_count(QVector<QVector<double> > matrix);
    static double determenant_matrix_2x2_count(int i,int j,QVector<QVector<double> > matrix_3x3);
    static QVector<QVector<double> > difference__matrix(QVector<QVector<double> > matrix1, QVector<QVector<double> > matrix2);
    static double evclidova_norma(QVector<QVector<double> > matrix);
    static double evclidova_norma(QVector<double> matrix);
    static QVector<QVector<double> > create_up_tiangle_matrix(QVector<QVector<double> > matrix_A,QVector<QVector<double> > matrix_V);

    static QVector<QVector<double> > reverse_matrix_count(QVector<QVector<double> > matrix);
    static QVector<QVector<double> > count_q_Khol (QVector <QVector<double> > matrix_A, QVector <QVector<double> > matrix_Ph, QVector <QVector<double> > matrix_H);
    static QVector<QVector<double> > count_q_Gr_Sh (QVector <QVector<double> > matrix_A, QVector <QVector<double> > matrix_Ph, QVector <QVector<double> > matrix_H);
};

#endif // NORMAL_SYSTEM_EQUATIONS_H
