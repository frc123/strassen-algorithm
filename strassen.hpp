#include "matrix.hpp"

template <typename T>
UniqueMatrix<T> StrassenMultiplicationAux
    (const Matrix<T>& a, const Matrix<T>& b)
{
    size_t size, middle;
    size = a.Rows();
    if (size == 1)
    {
        return UniqueMatrix<T>(typename Matrix<T>::Data(1, 
            std::vector<T>(1, a.At(0, 0) * b.At(0, 0))));
    }
    middle = size >> 1;
    UniqueMatrix<T> c(size, size);
    const SharedMatrix<T> a_11 = a.ShareMatrix(0, middle, 0, middle);
    const SharedMatrix<T> a_12 = a.ShareMatrix(0, middle, middle, middle);
    const SharedMatrix<T> a_21 = a.ShareMatrix(middle, middle, 0, middle);
    const SharedMatrix<T> a_22 = a.ShareMatrix(middle, middle, middle, middle);
    const SharedMatrix<T> b_11 = b.ShareMatrix(0, middle, 0, middle);
    const SharedMatrix<T> b_12 = b.ShareMatrix(0, middle, middle, middle);
    const SharedMatrix<T> b_21 = b.ShareMatrix(middle, middle, 0, middle);
    const SharedMatrix<T> b_22 = b.ShareMatrix(middle, middle, middle, middle);
    UniqueMatrix<T> s_1 = b_12 - b_22;
    UniqueMatrix<T> s_2 = a_11 + a_12;
    UniqueMatrix<T> s_3 = a_21 + a_22;
    UniqueMatrix<T> s_4 = b_21 - b_11;
    UniqueMatrix<T> s_5 = a_11 + a_22;
    UniqueMatrix<T> s_6 = b_11 + b_22;
    UniqueMatrix<T> s_7 = a_12 - a_22;
    UniqueMatrix<T> s_8 = b_21 + b_22;
    UniqueMatrix<T> s_9 = a_11 - a_21;
    UniqueMatrix<T> s_10 = b_11 + b_12;
    UniqueMatrix<T> p_1 = StrassenMultiplicationAux(a_11, s_1); 
    UniqueMatrix<T> p_2 = StrassenMultiplicationAux(s_2, b_22); 
    UniqueMatrix<T> p_3 = StrassenMultiplicationAux(s_3, b_11); 
    UniqueMatrix<T> p_4 = StrassenMultiplicationAux(a_22, s_4); 
    UniqueMatrix<T> p_5 = StrassenMultiplicationAux(s_5, s_6); 
    UniqueMatrix<T> p_6 = StrassenMultiplicationAux(s_7, s_8); 
    UniqueMatrix<T> p_7 = StrassenMultiplicationAux(s_9, s_10); 
    SharedMatrix<T> c_11 = c.ShareMatrix(0, middle, 0, middle);
    SharedMatrix<T> c_12 = c.ShareMatrix(0, middle, middle, middle);
    SharedMatrix<T> c_21 = c.ShareMatrix(middle, middle, 0, middle);
    SharedMatrix<T> c_22 = c.ShareMatrix(middle, middle, middle, middle);
    c_11 += p_5;
    c_11 += p_4;
    c_11 -= p_2;
    c_11 += p_6;
    c_12 += p_1;
    c_12 += p_2;
    c_21 += p_3;
    c_21 += p_4;
    c_22 += p_5;
    c_22 += p_1;
    c_22 -= p_3;
    c_22 -= p_7;
    return c;
}

template <typename T>
UniqueMatrix<T> StrassenMultiplication
    (const Matrix<T>& a, const Matrix<T>& b)
{
    size_t size;
    size = a.Rows();
    if (size != a.Cols() || b.Rows() != b.Cols() 
        || size != b.Rows())
        throw std::runtime_error("Strassen algorithm on matrices with invalid sizes");
    if ((size != 0) && (!(size & (size - 1))) == false)
        throw std::runtime_error("Strassen algorithm on matrices without power of two");
    return StrassenMultiplicationAux(a, b);
}
