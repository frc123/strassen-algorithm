#ifndef _MATRIX_HPP
#define _MATRIX_HPP

#include <vector>
#include <list>
#include <stdexcept>

enum MatrixType { kUnique, kShared };

template <typename T>
class UniqueMatrix;

template <typename T>
class SharedMatrix;

template <typename T>
class Matrix
{
public:
    class Row
    {
    public:
        friend Matrix;
        
        Row() = delete;
        Row(const Row& copy) = delete;
        Row(Row&& move);
        T& operator[](size_t col);
    private:
        Row(Matrix& matrix, size_t row);

        Matrix& matrix_;
        size_t row_;
    };

    class ConstRow
    {
    public:
        friend Matrix;
        
        ConstRow() = delete;
        ConstRow(const ConstRow& copy) = delete;
        ConstRow(ConstRow&& move);
        const T& operator[](size_t col) const;
    private:
        ConstRow(const Matrix& matrix, size_t row);

        const Matrix& matrix_;
        size_t row_;
    };
    
    friend Row;
    friend ConstRow;

    using Data = std::vector< std::vector<T> >;

    virtual ~Matrix() = 0;
    // const T& GetElement(size_t row, size_t col) const;
    // void SetElement(size_t row, size_t col, const T& element);
    // void SetElement(size_t row, size_t col, T&& element);
    Row operator[](size_t row);
    const ConstRow operator[](size_t row) const;
    T& At(size_t row, size_t col);
    const T& At(size_t row, size_t col) const;
    UniqueMatrix<T> operator+(const Matrix<T>& other) const;
    UniqueMatrix<T> operator-(const Matrix<T>& other) const;
    UniqueMatrix<T> operator*(const Matrix<T>& other) const;
    Matrix<T>& operator+=(const Matrix<T>& other);
    Matrix<T>& operator-=(const Matrix<T>& other);
    bool operator==(const Matrix<T>& other) const;
    virtual SharedMatrix<T> ShareMatrix
        (size_t row_base, size_t row_num, size_t col_base, size_t col_num) = 0;
    virtual const SharedMatrix<T> ShareMatrix
        (size_t row_base, size_t row_num, size_t col_base, size_t col_num) const = 0;
    virtual bool Valid() const = 0;
    virtual MatrixType Type() const = 0;
    size_t Rows() const;
    size_t Cols() const;

protected:
    using SharedMatrixPool = std::list< SharedMatrix<T>** >;

    virtual T& Element(size_t row, size_t col) = 0;
    virtual const T& Element(size_t row, size_t col) const = 0;

private:
    Data Addition(const Matrix<T>& other) const;
    Data Subtraction(const Matrix<T>& other) const;
    Data Multiplication(const Matrix<T>& other) const;
    void AdditionAssign(const Matrix<T>& other);
    void SubtractionAssign(const Matrix<T>& other);

protected:
    size_t row_num_;
    size_t col_num_;
};

template <typename T>
class UniqueMatrix : public Matrix<T>
{
public:
    friend class SharedMatrix<T>;

    UniqueMatrix() = delete;
    UniqueMatrix(const UniqueMatrix& copy) = delete;
    UniqueMatrix(UniqueMatrix&& move);
    UniqueMatrix(const typename Matrix<T>::Data& data);
    UniqueMatrix(typename Matrix<T>::Data&& data);
    UniqueMatrix(size_t row_num, size_t col_num);
    UniqueMatrix(std::initializer_list< std::vector<T> > init);
    ~UniqueMatrix();
    SharedMatrix<T> ShareMatrix
        (size_t row_base, size_t row_num, size_t col_base, size_t col_num);
    const SharedMatrix<T> ShareMatrix
        (size_t row_base, size_t row_num, size_t col_base, size_t col_num) const;
    bool Valid() const;
    MatrixType Type() const;
    void ResizeRow(size_t num);
    void ResizeCol(size_t num);

protected:
    T& Element(size_t row, size_t col);
    const T& Element(size_t row, size_t col) const;

private:
    UniqueMatrix** p_p_u_matrix_;
    typename Matrix<T>::Data data_;
    typename Matrix<T>::SharedMatrixPool s_matrix_pool_;
};

template <typename T>
class SharedMatrix : public Matrix<T>
{
public:
    friend class UniqueMatrix<T>;

    SharedMatrix() = delete;
    SharedMatrix(const SharedMatrix& copy) = delete;
    SharedMatrix(SharedMatrix&& move);
    ~SharedMatrix();
    SharedMatrix<T> ShareMatrix
        (size_t row_base, size_t row_num, size_t col_base, size_t col_num);
    const SharedMatrix<T> ShareMatrix
        (size_t row_base, size_t row_num, size_t col_base, size_t col_num) const;
    bool Valid() const;
    MatrixType Type() const;

protected:
    T& Element(size_t row, size_t col);
    const T& Element(size_t row, size_t col) const;

private:
    SharedMatrix(UniqueMatrix<T>** p_p_u_matrix,
        size_t row_base, size_t row_num, size_t col_base, size_t col_num);

    UniqueMatrix<T>** p_p_u_matrix_;
    SharedMatrix<T>** p_p_s_matrix_;
    typename Matrix<T>::SharedMatrixPool::iterator s_matrix_pool_it_;
    size_t row_base_;
    size_t col_base_;
};

template <typename T>
Matrix<T>::Row::Row(Row&& move) : matrix_(move.matrix_), row_(move.row_) {}

template <typename T>
Matrix<T>::Row::Row(Matrix& matrix, size_t row) : matrix_(matrix), row_(row) {}

template <typename T>
T& Matrix<T>::Row::operator[](size_t col)
{
    if (col >= matrix_.col_num_ || col < 0)
        throw std::out_of_range("invalid column index");
    return matrix_.Element(row_, col);
}

template <typename T>
Matrix<T>::ConstRow::ConstRow(ConstRow&& move) 
    : matrix_(move.matrix_), row_(move.row_) {}

template <typename T>
Matrix<T>::ConstRow::ConstRow(const Matrix& matrix, size_t row) 
    : matrix_(matrix), row_(row) {}

template <typename T>
const T& Matrix<T>::ConstRow::operator[](size_t col) const
{
    if (col >= matrix_.col_num_ || col < 0)
        throw std::out_of_range("invalid column index");
    return matrix_.Element(row_, col);
}

template <typename T>
typename Matrix<T>::Row Matrix<T>::operator[](size_t row)
{
    if (this->Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (row >= this->row_num_ || row < 0)
        throw std::out_of_range("invalid row index");
    return Row(*this, row);
}

template <typename T>
const typename Matrix<T>::ConstRow Matrix<T>::operator[](size_t row) const
{
    if (this->Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (row >= this->row_num_ || row < 0)
        throw std::out_of_range("invalid row index");
    return ConstRow(*this, row);
}

template <typename T>
T& Matrix<T>::At(size_t row, size_t col)
{
    if (this->Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (row >= this->row_num_ || row < 0)
        throw std::out_of_range("invalid row index");
    if (col >= this->col_num_ || col < 0)
        throw std::out_of_range("invalid column index");
    return Element(row, col);
}

template <typename T>
const T& Matrix<T>::At(size_t row, size_t col) const
{
    if (this->Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (row >= this->row_num_ || row < 0)
        throw std::out_of_range("invalid row index");
    if (col >= this->col_num_ || col < 0)
        throw std::out_of_range("invalid column index");
    return Element(row, col);
}

template <typename T>
Matrix<T>::~Matrix() {}

template <typename T>
UniqueMatrix<T> Matrix<T>::operator+(const Matrix<T>& other) const
{
    if (this->Valid() == false || other.Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (this->row_num_ != other.row_num_ || this->col_num_ != other.col_num_)
        throw std::runtime_error("addition on matrices with different sizes");
    return UniqueMatrix<T>(Addition(other));
}

template <typename T>
UniqueMatrix<T> Matrix<T>::operator-(const Matrix<T>& other) const
{
    if (this->Valid() == false || other.Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (this->row_num_ != other.row_num_ || this->col_num_ != other.col_num_)
        throw std::runtime_error("subtraction on matrices with different sizes");
    return UniqueMatrix<T>(Subtraction(other));
}

template <typename T>
UniqueMatrix<T> Matrix<T>::operator*(const Matrix<T>& other) const
{
    if (this->Valid() == false || other.Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (this->col_num_ != other.row_num_)
        throw std::runtime_error("multiplication on matrices with invalid sizes");
    return UniqueMatrix<T>(Multiplication(other));
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other)
{
    if (this->Valid() == false || other.Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (this->row_num_ != other.row_num_ || this->col_num_ != other.col_num_)
        throw std::runtime_error("addition on matrices with different sizes");
    AdditionAssign(other);
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other)
{
    if (this->Valid() == false || other.Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (this->row_num_ != other.row_num_ || this->col_num_ != other.col_num_)
        throw std::runtime_error("subtraction on matrices with different sizes");
    SubtractionAssign(other);
    return *this;
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& other) const
{
    size_t i, j;
    if (this->row_num_ != other.row_num_ || this->col_num_ != other.col_num_)
        return false;
    for (i = 0; i < this->row_num_; ++i)
    {
        for (j = 0; j < this->col_num_; ++j)
        {
            if (this->Element(i, j) != other.Element(i, j))
                return false;
        }
    }
    return true;
}

template <typename T>
size_t Matrix<T>::Rows() const
{
    return row_num_;
}

template <typename T>
size_t Matrix<T>::Cols() const
{
    return col_num_;
}

template <typename T>
typename Matrix<T>::Data Matrix<T>::Addition(const Matrix<T>& other) const
{
    size_t i, j;
    Data result(this->row_num_, std::vector<T>(this->col_num_));
    for (i = 0; i < this->row_num_; ++i)
    {
        for (j = 0; j < this->col_num_; ++j)
        {
            result[i][j] = this->Element(i, j) + other.Element(i, j);
        }
    }
    return result;
}

template <typename T>
typename Matrix<T>::Data Matrix<T>::Subtraction(const Matrix<T>& other) const
{
    size_t i, j;
    Data result(this->row_num_, std::vector<T>(this->col_num_));
    for (i = 0; i < this->row_num_; ++i)
    {
        for (j = 0; j < this->col_num_; ++j)
        {
            result[i][j] = this->Element(i, j) - other.Element(i, j);
        }
    }
    return result;
}

template <typename T>
typename Matrix<T>::Data Matrix<T>::Multiplication(const Matrix<T>& other) const
{
    size_t i, k, j;
    Data result(this->row_num_, std::vector<T>(other.col_num_));
    for (i = 0; i < this->row_num_; ++i)
    {
        for (j = 0; j < other.col_num_; ++j)
        {
            result[i][j] = 0;
            for (k = 0; k < this->col_num_; ++k)
            {
                result[i][j] += this->Element(i, k) * other.Element(k, j);
            }
        }
    }
    return result;
}

template <typename T>
void Matrix<T>::AdditionAssign(const Matrix<T>& other)
{
    size_t i, j;
    Data result(this->row_num_, std::vector<T>(this->col_num_));
    for (i = 0; i < this->row_num_; ++i)
    {
        for (j = 0; j < this->col_num_; ++j)
        {
            this->Element(i, j) += other.Element(i, j);
        }
    }
}

template <typename T>
void Matrix<T>::SubtractionAssign(const Matrix<T>& other)
{
    size_t i, j;
    Data result(this->row_num_, std::vector<T>(this->col_num_));
    for (i = 0; i < this->row_num_; ++i)
    {
        for (j = 0; j < this->col_num_; ++j)
        {
            this->Element(i, j) -= other.Element(i, j);
        }
    }
}

template <typename T>
UniqueMatrix<T>::UniqueMatrix(UniqueMatrix&& move)
    : data_(std::move(move.data_)), s_matrix_pool_(std::move(move.s_matrix_pool_))
{
    this->p_p_u_matrix_ = move.p_p_u_matrix_;
    move.p_p_u_matrix_ = new UniqueMatrix*();
    *(move.p_p_u_matrix_) = &move;
    *(this->p_p_u_matrix_) = this;
    this->row_num_ = move.row_num_;
    move.row_num_ = 0;
    this->col_num_ = move.col_num_;
    move.col_num_ = 0;
}

template <typename T>
UniqueMatrix<T>::UniqueMatrix(const typename Matrix<T>::Data& data)
    : data_(data)
{
    p_p_u_matrix_ = new UniqueMatrix*();
    *p_p_u_matrix_ = this;
    this->row_num_ = data_.size();
    this->col_num_ = (this->row_num_ == 0) ? 0 : data_[0].size();
}

template <typename T>
UniqueMatrix<T>::UniqueMatrix(typename Matrix<T>::Data&& data)
    : data_(std::move(data))
{
    p_p_u_matrix_ = new UniqueMatrix*();
    *p_p_u_matrix_ = this;
    this->row_num_ = data_.size();
    this->col_num_ = (this->row_num_ == 0) ? 0 : data_[0].size();
}

template <typename T>
UniqueMatrix<T>::UniqueMatrix(size_t row_num, size_t col_num)
    : data_(typename Matrix<T>::Data(row_num, std::vector<T>(col_num)))
{
    p_p_u_matrix_ = new UniqueMatrix*();
    *p_p_u_matrix_ = this;
    this->row_num_ = row_num;
    this->col_num_ = col_num;
}

template <typename T>
UniqueMatrix<T>::UniqueMatrix(std::initializer_list< std::vector<T> > init)
    : data_(init)
{
    p_p_u_matrix_ = new UniqueMatrix*();
    *p_p_u_matrix_ = this;
    this->row_num_ = data_.size();
    this->col_num_ = (this->row_num_ == 0) ? 0 : data_[0].size();
}

template <typename T>
UniqueMatrix<T>::~UniqueMatrix()
{
    delete p_p_u_matrix_;
    for (SharedMatrix<T>** p_p_s_matrix : s_matrix_pool_)
    {
        (*p_p_s_matrix)->p_p_u_matrix_ = nullptr;
    }
}

template <typename T>
SharedMatrix<T> UniqueMatrix<T>::ShareMatrix
    (size_t row_base, size_t row_num, size_t col_base, size_t col_num)
{
    if (row_num < 0 || row_base < 0 || row_base + row_num > this->row_num_)
        throw std::runtime_error("invalid row_base or row_num");
    if (col_num < 0 || col_base < 0 || col_base + col_num > this->col_num_)
        throw std::runtime_error("invalid col_base or col_num");
    SharedMatrix<T> s_matrix(p_p_u_matrix_, row_base, row_num, col_base, col_num);
    s_matrix_pool_.push_front(s_matrix.p_p_s_matrix_);
    s_matrix.s_matrix_pool_it_ = s_matrix_pool_.begin();
    return s_matrix;
}

template <typename T>
const SharedMatrix<T> UniqueMatrix<T>::ShareMatrix
    (size_t row_base, size_t row_num, size_t col_base, size_t col_num) const
{
    if (row_num < 0 || row_base < 0 || row_base + row_num > this->row_num_)
        throw std::runtime_error("invalid row_base or row_num");
    if (col_num < 0 || col_base < 0 || col_base + col_num > this->col_num_)
        throw std::runtime_error("invalid col_base or col_num");
    SharedMatrix<T> s_matrix(const_cast<UniqueMatrix<T>**>(p_p_u_matrix_), 
        row_base, row_num, col_base, col_num);
    const_cast<UniqueMatrix<T>*>(this)->
        s_matrix_pool_.push_front(const_cast<SharedMatrix<T>**>(s_matrix.p_p_s_matrix_));
    s_matrix.s_matrix_pool_it_ = const_cast<UniqueMatrix<T>*>(this)->s_matrix_pool_.begin();
    return s_matrix;
}

template <typename T>
bool UniqueMatrix<T>::Valid() const
{
    return true;
}

template <typename T>
MatrixType UniqueMatrix<T>::Type() const
{
    return MatrixType::kUnique;
}

template <typename T>
T& UniqueMatrix<T>::Element(size_t row, size_t col)
{
    return data_[row][col];
}

template <typename T>
const T& UniqueMatrix<T>::Element(size_t row, size_t col) const
{
    return data_[row][col];
}

template <typename T>
void UniqueMatrix<T>::ResizeRow(size_t num)
{
    this->row_num_ = num;
    this->data_.resize(this->row_num_);
}

template <typename T>
void UniqueMatrix<T>::ResizeCol(size_t num)
{
    this->col_num_ = num;
    for (std::vector<T>& row : this->data_)
    {
        row.resize(this->col_num_);
    }
}

template <typename T>
SharedMatrix<T>::SharedMatrix(SharedMatrix&& move)
    : s_matrix_pool_it_(std::move(move.s_matrix_pool_it_))
{
    this->p_p_s_matrix_ = move.p_p_s_matrix_;
    move.p_p_s_matrix_ = new SharedMatrix*();
    *(move.p_p_s_matrix_) = &move;
    *(this->p_p_s_matrix_) = this;
    this->row_num_ = move.row_num_;
    move.row_num_ = 0;
    this->col_num_ = move.col_num_;
    move.col_num_ = 0;
    row_base_ = move.row_base_;
    move.row_base_ = 0;
    col_base_ = move.col_base_;
    move.col_base_ = 0;
    this->p_p_u_matrix_ = move.p_p_u_matrix_;
    move.p_p_u_matrix_ = nullptr;
}

template <typename T>
SharedMatrix<T>::~SharedMatrix()
{
    delete p_p_s_matrix_;
    if (Valid())
    {
        (*p_p_u_matrix_)->s_matrix_pool_.erase(s_matrix_pool_it_);
    }
}

template <typename T>
SharedMatrix<T> SharedMatrix<T>::ShareMatrix
    (size_t row_base, size_t row_num, size_t col_base, size_t col_num)
{
    if (this->Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (row_num < 0 || row_base < 0 || row_base + row_num > this->row_num_)
        throw std::runtime_error("invalid row_base or row_num");
    if (col_num < 0 || col_base < 0 || col_base + col_num > this->col_num_)
        throw std::runtime_error("invalid col_base or col_num");
    SharedMatrix s_matrix(p_p_u_matrix_, this->row_base_ + row_base, 
        row_num, this->col_base_ + col_base, col_num);
    (*p_p_u_matrix_)->s_matrix_pool_.push_front(s_matrix.p_p_s_matrix_);
    s_matrix.s_matrix_pool_it_ = (*p_p_u_matrix_)->s_matrix_pool_.begin();
    return s_matrix;
}

template <typename T>
const SharedMatrix<T> SharedMatrix<T>::ShareMatrix
    (size_t row_base, size_t row_num, size_t col_base, size_t col_num) const
{
    if (this->Valid() == false)
        throw std::runtime_error("invalid matrix");
    if (row_num < 0 || row_base < 0 || row_base + row_num > this->row_num_)
        throw std::runtime_error("invalid row_base or row_num");
    if (col_num < 0 || col_base < 0 || col_base + col_num > this->col_num_)
        throw std::runtime_error("invalid col_base or col_num");
    SharedMatrix<T> s_matrix(const_cast<UniqueMatrix<T>**>(p_p_u_matrix_), this->row_base_ + row_base, 
        row_num, this->col_base_ + col_base, col_num);
    (*p_p_u_matrix_)->s_matrix_pool_.push_front(s_matrix.p_p_s_matrix_);
    s_matrix.s_matrix_pool_it_ = (*p_p_u_matrix_)->s_matrix_pool_.begin();
    return s_matrix;
}

template <typename T>
bool SharedMatrix<T>::Valid() const
{
    return p_p_u_matrix_ != nullptr;
}

template <typename T>
MatrixType SharedMatrix<T>::Type() const
{
    return MatrixType::kShared;
}

template <typename T>
T& SharedMatrix<T>::Element(size_t row, size_t col)
{
    return ((*p_p_u_matrix_)->data_)[row_base_ + row][col_base_ + col];
}

template <typename T>
const T& SharedMatrix<T>::Element(size_t row, size_t col) const
{
    return ((*p_p_u_matrix_)->data_)[row_base_ + row][col_base_ + col];
}

template <typename T>
SharedMatrix<T>::SharedMatrix(UniqueMatrix<T>** p_p_u_matrix,
        size_t row_base, size_t row_num, size_t col_base, size_t col_num)
{
    this->p_p_s_matrix_ = new SharedMatrix*();
    *(this->p_p_s_matrix_) = this;
    this->row_num_ = row_num;
    this->col_num_ = col_num;
    this->row_base_ = row_base;
    this->col_base_ = col_base;
    this->p_p_u_matrix_ = p_p_u_matrix;
}

#endif