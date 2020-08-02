#include "normal_system_equations.h"
#include <math.h>

normal_system_equations::normal_system_equations()
{

}

QVector<double> normal_system_equations::multiplication_matrix_and_number(double number, QVector<double> matrix)
{
    for (int i = 0;i<matrix.size();i++)
        matrix[i]*=number;
    return matrix;
}

QVector<double> normal_system_equations::projection_count(QVector<double> a, QVector<double> b)
{
    return multiplication_matrix_and_number(scalar_multiplication_string_matrix(a,b)/scalar_multiplication_string_matrix(b,b),b);
}

double normal_system_equations::scalar_multiplication_string_matrix(QVector<double> matrix1, QVector<double> matrix2)
{
    double result = 0;
    for (int i = 0;i<matrix1.size();i++)
        result += matrix1[i]*matrix2[i];
    return result;
}

QVector<QVector<double> > normal_system_equations::multiplication_matrix(QVector<QVector<double> > matrix1, QVector<QVector<double> > matrix2)
{
    QVector<QVector<double> > result(matrix2.size(),QVector<double> (matrix1[0].size(),0));
    for (int i = 0;i<matrix1[0].size();i++)//i - количество строк внутри столбца первой матрицы
        for (int j = 0;j<matrix2.size();j++)//j - количество столбцов второй матрицы
        {
            for (int k = 0;k<matrix1.size();k++)//k - количество столбцов первой матрицы
                 result[j][i]+=matrix1[k][i]*matrix2[j][k];
        }
    return result;
}

QVector<QVector<double> > normal_system_equations::transposition_matrix(QVector<QVector<double> > matrix)
{
    QVector<double> non;
    QVector<QVector<double> > result(matrix[0].size(),non);
    for (int i = 0;i<matrix.size();i++)
        for (int j = 0;j<matrix[i].size();j++)
            result[j].push_back(matrix[i][j]);
    return result;
}

QVector<QVector<double> > normal_system_equations::reverse_matrix_3x3_count(QVector<QVector<double> > matrix)
{
    QVector<QVector<double> > result(matrix.size(),QVector<double> (matrix[0].size(),1/determenant_matrix_3x3_count(matrix)));
    for (int i = 0;i<matrix.size();i++)
        for (int j = 0;j<matrix[0].size();j++)
            result[i][j] *= pow(-1,i+j)*determenant_matrix_2x2_count(i,j,matrix);
    return normal_system_equations::transposition_matrix(result);
}

double normal_system_equations::determenant_matrix_3x3_count(QVector<QVector<double> > matrix)
{
    return matrix[0][0]*matrix[1][1]*matrix[2][2]
            +matrix[1][0]*matrix[2][1]*matrix[0][2]
            +matrix[0][1]*matrix[1][2]*matrix[2][0]
            -matrix[2][0]*matrix[1][1]*matrix[0][2]
            -matrix[0][0]*matrix[1][2]*matrix[2][1]
            -matrix[0][1]*matrix[1][0]*matrix[2][2];
}

double normal_system_equations::determenant_matrix_2x2_count(int i,int j,QVector<QVector<double> > matrix_3x3)
{
    QVector<QVector<double> > result;
    for (int ii = 0;ii<matrix_3x3.size();ii++)
        if (ii!=i)
        {
            result.push_back(QVector<double> (0));
            for (int jj = 0;jj<matrix_3x3[0].size();jj++)
                if (jj!=j)
                {
                    result[result.size()-1].push_back(matrix_3x3[ii][jj]);
                }
        }
    return result[0][0]*result[1][1]-result[0][1]*result[1][0];
}

QVector<QVector<double> > normal_system_equations::difference__matrix(QVector<QVector<double> > matrix1, QVector<QVector<double> > matrix2)
{
    QVector<QVector<double> > result(matrix1.size(),QVector<double> (matrix1[0].size(),0));
    for (int i = 0;i<matrix1.size();i++)
        for (int j = 0;j<matrix1[0].size();j++)
            result[i][j]=matrix1[i][j]-matrix2[i][j];
    return result;
}

double normal_system_equations::evclidova_norma(QVector<QVector<double> > matrix)
{
    double result = 0;
    for (int i = 0; i < matrix[0].size(); ++i) {
        for (int j = 0; j < matrix.size(); ++j) {
            result += matrix[j][i]*matrix[j][i];
        }
    }
    return sqrt(result);
}

double normal_system_equations::evclidova_norma(QVector<double> matrix)
{
    double result = 0;
    for (int j = 0; j < matrix.size(); ++j) {
        result += matrix[j]*matrix[j];
        }
    return sqrt(result);
}

QVector<QVector<double> > normal_system_equations::create_up_tiangle_matrix(QVector<QVector<double> > matrix_A, QVector<QVector<double> > matrix_V)
{
    QVector<QVector<double> > matrix_R (matrix_A.size(),QVector <double> (matrix_A.size(),0));
    for (int i = 0; i < matrix_R.size(); ++i)
    {
        for (int j = 0; j < i+1; ++j)
        {
            matrix_R[i][j] = scalar_multiplication_string_matrix(matrix_A[i],matrix_V[j]);
        }
    }
    return matrix_R;
}

QVector<double> normal_system_equations::difference_string_matrix(QVector<double> matrix1, QVector<double> matrix2)
{
    for (int i = 0;i<matrix1.size();i++)
        matrix1[i]-=matrix2[i];
    return matrix1;
}

QVector<QVector<double> > normal_system_equations::ortogonalization_and_ortonormal_matrix(QVector<QVector<double> > matrix)
{
    QVector<QVector<double> > B(matrix.size(),QVector<double> ());
    double norma;
    for (int i = 0;i<matrix.size();i++)
        B[i] = count_b(i+1,matrix);
    for (int i = 0;i<B.size();i++){
        norma = evclidova_norma(B[i]);
        for (int j = 0;j<B[0].size();j++)
            B[i][j] /= norma;
    }
    return B;
}

QVector<double> normal_system_equations::count_b(int i, QVector<QVector<double> > A)
{
    if (i-1==0)
        return A[0];
    else
    {
        for (int j = 0;j<i-1;j++)
            A[i-1] = difference_string_matrix(A[i-1],projection_count(A[i-1],count_b(j+1,A)));
        return A[i-1];
    }
}



QVector<QVector<double> > normal_system_equations::reverse_matrix_count(QVector<QVector<double> > matrix)
{
    QVector<QVector<double> > A11(matrix.size()/2,QVector<double>(matrix.size()/2,0)),A12(matrix.size()/2,QVector<double>(matrix.size()/2,0)),A21(matrix.size()/2,QVector<double>(matrix.size()/2,0)),A22(matrix.size()/2,QVector<double>(matrix.size()/2,0)),rA22(matrix.size()/2,QVector<double>(matrix.size()/2,0)),rA11(matrix.size()/2,QVector<double>(matrix.size()/2,0)),C11(matrix.size()/2,QVector<double>(matrix.size()/2,0)),C12(matrix.size()/2,QVector<double>(matrix.size()/2,0)),C21(matrix.size()/2,QVector<double>(matrix.size()/2,0)),C22 (matrix.size()/2,QVector<double>(matrix.size()/2,0));
    QVector<QVector<double> > nul_matrix (3,QVector<double> (3,0));
    for (int i = 0;i<matrix.size()/2;i++)
        for (int j = 0;j<matrix.size()/2;j++)
        {
            A11[i][j] = matrix[i][j];
            A12[i][j] = matrix[i+matrix.size()/2][j];
            A21[i][j] = matrix[i][j+matrix.size()/2];
            A22[i][j] = matrix[i+matrix.size()/2][j+matrix.size()/2];
        }
    rA11 = reverse_matrix_3x3_count(A11);
    rA22 = reverse_matrix_3x3_count(A22);
    C11 = reverse_matrix_3x3_count(
                difference__matrix(A11,multiplication_matrix(multiplication_matrix(A12,rA22),A21)));
    C21 = difference__matrix(
                nul_matrix,multiplication_matrix(
                    multiplication_matrix(rA22,A21),C11));
    C22 = reverse_matrix_3x3_count(
                difference__matrix(A22,multiplication_matrix(
                                       multiplication_matrix(A21,rA11),A12)));
    C12 = difference__matrix(
                nul_matrix,multiplication_matrix(
                    multiplication_matrix(rA11,A12),C22));
    for (int i = 0;i<matrix.size()/2;i++)
        for (int j = 0;j<matrix.size()/2;j++)
        {
            matrix[i][j] = C11[i][j];
            matrix[i+matrix.size()/2][j] = C12[i][j];
            matrix[i][j+matrix.size()/2] = C21[i][j];
            matrix[i+matrix.size()/2][j+matrix.size()/2] = C22[i][j];
        }
    return matrix;
}

QVector<QVector<double> > normal_system_equations::count_q_Khol(QVector<QVector<double> > matrix_A, QVector<QVector<double> > matrix_Ph, QVector<QVector<double> > matrix_H)
{
    for (int i = 0;i<matrix_Ph.size();i++)
        matrix_Ph[i][i] = pow(matrix_Ph[i][i],0.5);

    matrix_A = multiplication_matrix(matrix_Ph,matrix_A);
    QVector<QVector<double> > matrix_B = multiplication_matrix(matrix_Ph,matrix_H);
    QVector <double> result (matrix_A.size(),0);
    QVector<QVector<double> > matrix_C = multiplication_matrix(transposition_matrix(matrix_A),matrix_A);
    QVector<double> matrix_dH = multiplication_matrix(transposition_matrix(matrix_A),matrix_B)[0];

    QVector<QVector<double> > matrix_U (matrix_A.size(),result);
    for (int j = 0; j < matrix_U[0].size(); ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            matrix_U[j][j] += pow(matrix_U[j][k],2);
        }
        matrix_U[j][j] = pow(matrix_C[j][j]-matrix_U[j][j],0.5);
        for (int nu = j+1; nu < matrix_U.size(); ++nu)
        {
            for (int k = 0; k < j; ++k)
            {
                matrix_U[nu][j] += matrix_U[nu][k]*matrix_U[j][k];
            }
            matrix_U[nu][j] = (matrix_C[nu][j] - matrix_U[nu][j])/matrix_U[j][j];
        }
    }
    QVector<QVector<double> > matrix_Utr = transposition_matrix(matrix_U);
    QVector<double> matrix_Y (matrix_A.size(),0);

    for (int i = 0; i < matrix_Y.size(); ++i)
    {
        matrix_Y[i] += matrix_dH[i];
        for (int j = 0; j < i; ++j)
        {
            matrix_Y[i] -= matrix_Y[j]*matrix_Utr[j][i];
        }
        matrix_Y[i] /= matrix_Utr[i][i];
    }


    for (int i = 0; i < result.size(); ++i)
    {
        result[result.size()-1-i] += matrix_Y[result.size()-1-i];
        for (int j = 0; j < i; ++j)
        {
            result[result.size()-1-i] -= result[result.size()-1-j]*matrix_U[result.size()-1-j][result.size()-1-i];
        }
        result[result.size()-1-i] /= matrix_U[result.size()-1-i][result.size()-1-i];
    }

    return QVector<QVector<double> > (1,result);
}

QVector<QVector<double> > normal_system_equations::count_q_Gr_Sh(QVector<QVector<double> > matrix_A, QVector<QVector<double> > matrix_Ph, QVector<QVector<double> > matrix_H)
{
    QVector <double> result (matrix_A.size(),0);
    for (int i = 0;i<matrix_Ph.size();i++)
        matrix_Ph[i][i] = sqrt(matrix_Ph[i][i]);

    matrix_A = multiplication_matrix(matrix_Ph,matrix_A);

    QVector<QVector<double> > matrix_V = ortogonalization_and_ortonormal_matrix(matrix_A);
    QVector<QVector<double> > Vt_multi_B = multiplication_matrix(transposition_matrix(matrix_V),multiplication_matrix(matrix_Ph,matrix_H));
    QVector<QVector<double> > matrix_R = create_up_tiangle_matrix(matrix_A,matrix_V);

    for (int i = 0; i < matrix_R.size(); ++i)
    {
        result[result.size()-1-i] += Vt_multi_B[0][result.size()-1-i];
        for (int j = 0; j < i; ++j)
        {
            result[result.size()-1-i] -= result[result.size()-1-j]*matrix_R[result.size()-1-j][result.size()-1-i];
        }
        result[result.size()-1-i] /= matrix_R[result.size()-1-i][result.size()-1-i];
    }
    return QVector<QVector<double> > (1,result);
}


