#if !defined(_OMPA_)

using std::vector;
using std::pair;
using std::make_pair;

extern vector<SimpleVector<myfloat> > filterM;

template <typename T> SimpleVector<T> particle(const vector<pair<SimpleMatrix<T>, SimpleVector<int> > >& conn, const int& sps = 80) {
  assert(0 < conn.size() && 0 < sps);
  SimpleVector<T> spector(sps);
  spector.O();
  if(conn.size() == 1) {
    assert(0 < conn[0].second);
    return spector;
  }
  for(int i = 0; i < conn.size(); i ++)
    spector += particle(conn[i], sps);
  return spector /= T(conn.size());
}

template <typename T> SimpleVector<T> calibrate(const SimpleVector<T>& in, const vector<pair<SimpleMatrix<T>, SimpleVector<int> > >& cal) {
  const auto c(particle(cal, in.size()));
  SimpleVector<T> res(in.size());
  for(int i = 0; i < c.size(); i ++)
    res[i] = c[i] / in[i];
  return res;
}

template <typename T> SimpleVector<T> filter(const T& speed, const SimpleVector<T>& forig) {
  ; // can we do this??
  return forig;
}

template <typename T> SimpleVector<T> mph(const vector<pair<SimpleVector<T>, T> >& in, const vector<vector<T> > col, const T& thresh = T(1) / T(1000)) {
  assert(col.size() == in.size());
  assert(filterM.size() == 3);
  // N.B. integrate (filterM) * light(freq) d(freq) == intensity for each.
  SimpleMatrix<T> test(in.size() * 3, in[0].first.size());
  SimpleVector<T> testb(test.rows());
  for(int i = 0; i < test.rows(); i ++) {
    assert(in[0].first.size() == in[i].first.size());
    assert(col[i].size() == 3);
    for(int j = 0; j < 3; j ++) {
      assert(filterM[j].size() == test.cols());
      test.row(i * 3 + j) = filter(in[i].second, filterM[j]);
      testb[i * 3 + j] = col[i][j];
    }
  }
  auto one(test.col(0));
  one.I(thresh);
  return test.inner(testb - one, testb + one);
}

template <typename T> SimpleVector<T> matchParticle(SimpleVector<T> in, const vector<pair<SimpleMatrix<T>, SimpleVector<int> > >& list, const SimpleVector<T>& cal, const T& thresh = T(1) / T(1000)) {
  assert(in.size() == cal.size());
  for(int i = 0; i < in.size(); i ++)
    in[i] *= cal[i];
  auto one(in);
  one.I(thresh);
  SimpleMatrix<T> test(in.size(), cal.size());
  for(int i = 0; i < cal.size(); i ++) {
    vector<pair<SimpleMatrix<T>, SimpleVector<int> > > work;
    work.resize(1, cal[i]);
    test.setCol(i, particle(work, in.size()));
  }
  return test.inner(in - one, in + one);
}

#define _OMPA_
#endif

