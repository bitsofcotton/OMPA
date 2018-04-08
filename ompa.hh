#if !defined(_OMPA_)

#include <Eigen/Core>
#include <vector>
#include <cmath>
#include "../konbu/konbu_init.h"

using std::vector;
using std::sqrt;

template <typename T> class OMPA {
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
  typedef struct {
    ;
  } moleculars_t;
  OMPA();
  ~OMPA();
  
  vector<vector<T> > initDicts(const vector<moleculars_t>& moleculars, const int& length) const;
  vector<T> calibrate(const vector<vector<T> >& input, const vector<vector<T> >& zero, const vector<vector<T> >& calibrate) const;
  vector<T> matchMixed(const vector<T>& input, const vector<vector<T> >& dicts) const;
  vector<T> rawmatch(const vector<T>& input, const vector<vector<T> >& dicts) const;
  
private:
  T thresh;
  T cut;
};

template <typename T> OMPA<T>::OMPA() {
  thresh = T(.001);
  cut    = T(.2);
}

template <typename T> OMPA<T>::~OMPA() {
  ;
}

template <typename T> vector<vector<T> > OMPA<T>::initDicts(const vector<moleculars_t>& moleculars, const int& length) const {
  // result is expectrd to raw radiowave per frequency from moleculars_t.
  // but moleculars_t is not so simple in general.
  vector<vector<T> > result;
  for(int i = 0; i < moleculars.size(); i ++) {
    vector<T> work;
    work.resize(length);
    // now stub around molecular_t.
    // with integrate each cell, we produce raw matrix from functions,
    // then, add and invert it with LP.
    // then, we get line from the matrix, then, we can get one solvee from PDE.
    LP<T> konbu;
    result.push_back(work);
  }
  return result;
}

template <typename T> vector<T> OMPA<T>::calibrate(const vector<vector<T> >& input, const vector<vector<T> >& zero, const vector<vector<T> >& calibrate) const {
  // calibrate input. input is expected to the pixel of photos with r\theta rotation.
  // result is expected to be raw radio wave with frequency from them.
  // now stub because we don't have such rotating mechanics...
  vector<T> result;
  return result;
}

template <typename T> vector<T> OMPA<T>::matchMixed(const vector<T>& input, const vector<vector<T> >& dicts) const {
  // result is expected to be the ratio of that seems to be the matter.
  // now stub.
  vector<T> result;
  return result;
}

template <typename T> vector<T> OMPA<T>::rawmatch(const vector<T>& input, const vector<vector<T> >& dicts) const {
  vector<T> result;
  LP<T> konbu;
  Mat A(input.size() * 2, dicts.size());
  Vec b(input.size() * 2), c;
  for(int i = 0; i < input.size(); i ++)
    for(int j = 0; j < dicts.size(); j ++) {
      A(2 * i + 0, j) =   dicts[i];
      A(2 * i + 1, j) = - dicts[i];
      b[2 * i + 0]    =   input[i];
      b[2 * i + 1]    = - input[i];
    }
  const T bnorm(sqrt(b.dot(b)));
  for(int i = 0; i < b.size(); i ++)
    b[i] += bnorm * thresh;
  konbu.inner(c, A, b);
  for(int i = 0; i < c.size(); i ++)
    result.push_back(c[i]);
  return result;
}

#define _OMPA_
#endif

