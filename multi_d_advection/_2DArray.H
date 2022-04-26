#ifndef _2DARRAY_H
#define _2DARRAY_H

#include <vector>
#include <cassert>
#include <ostream>

///
/// Creates a 2D array using 1D vector that takes number of rows, cols and the initial value.
///
class _2DArray {

private:
  std::size_t _rows;
  std::size_t _cols;
  std::vector<double> _data;

public:

  ///
  /// Parametric constructor that takes rows, cols and initial value.
  ///
  _2DArray(std::size_t rows, std::size_t cols, double val=0.0)
    : _rows{rows},
      _cols{cols},
      _data(rows * cols, val)
  {
    assert (rows > 0);
    assert (cols > 0);
  }

  ///
  /// Return number of columns
  ///
  inline std::size_t ncols() const {return _cols;}

  ///
  /// Return number of rows
  ///
  inline std::size_t nrows() const {return _rows;}

  ///
  /// Return the size of the 2D array
  ///
  inline std::size_t size() const {return _cols+_rows;}

  ///
  /// Return array number using array(row,col);
  /// 
  inline double& operator() (int row, int col){
    assert ( row >= 0 && row < _rows);
    assert ( col >= 0 && row < _cols);
    return _data[row*_cols + col];
  }

  ///
  /// make const function return const array value
  ///
  inline const double& operator() (int row, int col) const{
    assert ( row >= 0 && row < _rows);
    assert ( col >= 0 && row < _cols);
    return _data[row*_cols + col];
  }
  

  ///
  /// Prints out the array
  ///
  friend inline std::ostream& operator << (std::ostream& os, const _2DArray & a){
    for (std::size_t row=0; row < a.nrows(); ++row){
      for (std::size_t col=0; col < a.ncols(); ++col){
	os << a(row,col) << " ";
      }
      os << std::endl;
    }
    return os;
  }

};




#endif
