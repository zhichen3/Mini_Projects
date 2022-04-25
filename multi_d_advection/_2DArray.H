#ifndef _2DARRAY_H
#define _2DARRAY_H

#include <vector>
#include <cassert>
#include <ostream>

class _2DArray {

private:
  std::size_t _rows;
  std::size_t _cols;
  std::vector<double> _data;

public:
  _2DArray(std::size_t rows, std::size_t cols, double val=0.0)
    : _rows{rows},
      _cols{cols},
      _data(rows * cols, val)
  {
    assert (rows > 0);
    assert (cols > 0);
  }

  // Return number of columns and rows
  inline std::size_t ncols() const {return _cols;}
  inline std::size_t nrows() const {return _rows;}
  inline std::size_t size() const {return _cols+_rows;}

  // return array numberusing array(row,col);
  inline double& operator() (int row, int col){
    assert ( row >= 0 && row < _rows);
    assert ( col >= 0 && row < _cols);
    return _data[row*_cols + col];
  }
  
  // make const function return const array value
  inline const double& operator() (int row, int col) const{
    assert ( row >= 0 && row < _rows);
    assert ( col >= 0 && row < _cols);
    return _data[row*_cols + col];
  }
  

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
