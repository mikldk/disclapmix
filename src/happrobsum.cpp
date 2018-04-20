#include <Rcpp.h>
#include <unordered_map>
#include <functional>

#include <sys/resource.h> // mikl 2016-11-04

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

// [[Rcpp::plugins(cpp11)]]

/*
http://kevinushey.github.io/blog/2015/04/05/debugging-with-valgrind/
*/

using namespace Rcpp;


namespace{
  // a little helper that should IMHO be standardized
  template<typename T>
  std::size_t make_hash(const T& v){
    return std::hash<T>()(v);
  }
  
  // adapted from boost::hash_combine
  void hash_combine(std::size_t& h, const std::size_t& v){
    h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
  }
  
  // hash any container
  template<typename T>
  struct hash_container {
    // FIXME: mikl 2016-11-05: Is this working for only one locus?
    
    size_t operator()(const T& v) const{
      size_t h = 0;
      for( const auto& e : v ) {
        hash_combine(h, make_hash(e));
      }
      return h;
    }
  };
}

namespace std {
  // support for vector<T> if T is hashable
  // (the T... is a required trick if the vector has a non-standard allocator)
  template<>
  struct hash<IntegerVector> : hash_container<IntegerVector> {};
  
  struct equal_to_intvec : binary_function<IntegerVector, IntegerVector, bool> {
    bool operator() (const IntegerVector& x, const IntegerVector& y) const{
      if (x.size() != y.size()){
        return false;
      }
      
      int n = x.size();
      
      for (int i = 0; i < n; ++i) {
        if (x[i] != y[i]) {
          return false;
        }
      }
      
      return true;
    }
  };
  /*
  // the same for map<T,U> if T and U are hashable
  template<typename... T>
  struct hash<map<T...>> : hash_container<map<T...>> {};
  */
  // simply add more containers as needed
}

// FIXME: mikl 2016-11-05: Is this working for only one locus?
class MatchProbSums {
  private:
    int* m_counters;
    int* m_start_alleles;
    int* m_stop_alleles; 

    std::vector<List> m_fits;

    int n_fits;
    
    int loci;
    int* clusters;

    double* hap_sum;
    double* match_within;
    double** match_between;

    double* hap_probs(int* h) {
      double* ph = new double[n_fits];

      //Rcout << n_fits << std::endl;

      for (int f = 0; f < n_fits; ++f) {
        List fit = m_fits[f];

        IntegerMatrix y_f = fit[0];
        NumericMatrix p_f = fit[1];
        NumericVector tau_f = fit[2];

        double hprob = 0.0;

        int clus = clusters[f];

        for (int c = 0; c < clus; c++) {
          IntegerVector yhap = y_f(c, _);
          double component_prob = tau_f(c);
          
          for (int l = 0; l < loci; l++) {
            double pcl = p_f(c, l);
            component_prob *= pow(pcl, abs(h[l] - yhap(l)))*((1-pcl)/(1+pcl));
          }
          
          hprob += component_prob;
        }

        ph[f] = hprob;
      }

      return ph;
    }

    void do_work_with_h(int* h) {
      double* ph = hap_probs(h);

      for (int f = 0; f < n_fits; ++f) {
        hap_sum[f] += ph[f];
        match_within[f] += ph[f] * ph[f];
      }

      for (int f1 = 0; f1 < (n_fits-1); ++f1) {
        for (int f2 = (f1+1); f2 < n_fits; ++f2) {
          match_between[f1][f2] += ph[f1] * ph[f2];
        }
      }
    }
    
    void calc(int level) {
      if (level == loci) {
        do_work_with_h(m_counters);
      } else {
        for (int allele = m_start_alleles[level]; allele <= m_stop_alleles[level]; allele++) {
          m_counters[level] = allele;
          calc(level + 1);
        }
      }
    }

  public:
    MatchProbSums(List fits, IntegerMatrix allele_range) { 
      loci = allele_range.ncol();

      if (allele_range.nrow() != 2) {
        throw std::range_error("Exactly two rows in allele_range required");
      }

      n_fits = fits.size();

      clusters = new int[n_fits];
      hap_sum = new double[n_fits];
      match_within = new double[n_fits];
      match_between = new double*[n_fits];
      
      for (int f = 0; f < n_fits; ++f) {
        match_between[f] = new double[n_fits];
      }

      //print(fits);

      int f = 0;
      for (List::iterator it = fits.begin(); it != fits.end(); ++it ) {
        List l = *it;
        IntegerMatrix y = l["y"];
        NumericMatrix p = l["p"];
        NumericVector t = l["tau"];

        if (y.ncol() != loci) {
          throw std::range_error("Different number of loci (columns) in allele_range and y");
        }

        clusters[f] = y.nrow();
        List fit = List::create(y, p, t);
        m_fits.push_back(fit);
        f += 1;
      }

      for (int f = 0; f < n_fits; ++f) {
        hap_sum[f] = 0.0;
        match_within[f] = 0.0;
      }

      for (int f1 = 0; f1 < (n_fits-1); ++f1) {
        for (int f2 = (f1+1); f2 < n_fits; ++f2) {
          match_between[f1][f2] = 0.0;
        }
      }

      m_counters = new int[loci];
      m_start_alleles = new int[loci];
      m_stop_alleles = new int[loci];
      
      for (int l = 0; l < loci; l++) {
        m_start_alleles[l] = allele_range(0, l);
        m_stop_alleles[l] = allele_range(1, l);
      }

      //print(allele_range);
    }

    List calc() {
      calc(0);

      NumericVector res_hsum(n_fits);
      NumericVector res_within(n_fits);
      NumericMatrix res_between(n_fits, n_fits);

      for (int f = 0; f < n_fits; ++f) {
        res_hsum(f) = hap_sum[f];
        res_within(f) = match_within[f];
      }

      for (int f1 = 0; f1 < (n_fits-1); ++f1) {
        for (int f2 = (f1+1); f2 < n_fits; ++f2) {
          res_between(f1, f2) = match_between[f1][f2];
        }
      }

      List ret = List::create(
        _["hap_sum"] = res_hsum,
        _["match_within"] = res_within,
        _["match_between"] = res_between);

      return ret;
    }    
};

// FIXME: mikl 2016-11-05: Is this working for only one locus?
class MatchProbSumsCache {
  private:
    static const int CACHE_HAPLOTYPES_SUM_SIZE = 10000;

    std::vector<int> m_counters;
    std::vector<int> m_start_alleles;
    std::vector<int> m_stop_alleles; 

    std::vector< std::vector<int> > haplotypes_queue;
    int haplotypes_queue_counter;

    int loci;

    int n_fits;
    std::vector<IntegerMatrix> ys;
    std::vector<NumericMatrix> ps;
    std::vector<NumericVector> taus;
    std::vector<int> clusters;

    std::vector<double> hap_sum;
    std::vector<double> match_within;
    std::vector< std::vector<double> > match_between;

    void do_work_with_haplotype_queue() {
      if (haplotypes_queue_counter <= 0) {
        return;
      }

      std::vector< std::vector<double> > haplotypes_queue_ps(n_fits);

      for (int f = 0; f < n_fits; ++f) {
        haplotypes_queue_ps[f].resize(haplotypes_queue_counter);

        IntegerMatrix y_f = ys[f];
        NumericMatrix p_f = ps[f];
        NumericVector tau_f = taus[f];

        int clus = clusters[f];

        for (int h = 0; h < haplotypes_queue_counter; ++h) {
          std::vector<int> htype = haplotypes_queue[h];
          double hprob = 0.0;

          for (int c = 0; c < clus; c++) {
            IntegerVector yhap = y_f(c, _);
            double component_prob = tau_f(c);
            
            for (int l = 0; l < loci; l++) {
              double pcl = p_f(c, l);
              component_prob *= pow(pcl, abs(htype[l] - yhap(l)))*((1-pcl)/(1+pcl));
            }
            
            hprob += component_prob;
          }

          haplotypes_queue_ps[f][h] = hprob;
        }
      }

      for (int h = 0; h < haplotypes_queue_counter; ++h) {
        std::vector<double> ph(n_fits);

        for (int f = 0; f < n_fits; ++f) {
          ph[f] = haplotypes_queue_ps[f][h];
        }

        for (int f = 0; f < n_fits; ++f) {
          hap_sum[f] += ph[f];
          match_within[f] += ph[f] * ph[f];
        }

        for (int f1 = 0; f1 < (n_fits-1); ++f1) {
          for (int f2 = (f1+1); f2 < n_fits; ++f2) {
            match_between[f1][f2] += ph[f1] * ph[f2];
          }
        }
      }
    }
    
    void calc(int level) {
      if (level == loci) {
        if (haplotypes_queue_counter >= CACHE_HAPLOTYPES_SUM_SIZE) {
          do_work_with_haplotype_queue();

          // Not necessary because the elements are just reused
          //haplotypes_queue.clear();
          //haplotypes_queue.size(CACHE_HAPLOTYPES_SUM_SIZE);
          haplotypes_queue_counter = 0;
        }

        haplotypes_queue[haplotypes_queue_counter++] = m_counters;
      } else {
        for (int allele = m_start_alleles[level]; allele <= m_stop_alleles[level]; allele++) {
          m_counters[level] = allele;
          calc(level + 1);
        }
      }
    }

  public:
    MatchProbSumsCache(const List& fits, const IntegerMatrix& allele_range) { 
      loci = allele_range.ncol();

      if (allele_range.nrow() != 2) {
        throw std::range_error("Exactly two rows in allele_range required");
      }

      n_fits = fits.size();

      clusters.resize(n_fits);

      ys.resize(n_fits);
      ps.resize(n_fits);
      taus.resize(n_fits);

      hap_sum.resize(n_fits);
      match_within.resize(n_fits);
      match_between.resize(n_fits);
      
      for (int f = 0; f < n_fits; ++f) {
        match_between[f].resize(n_fits);
      }

      //print(fits);

      int f = 0;
      for (List::const_iterator it = fits.begin(); it != fits.end(); ++it ) {
        List l = *it;
        IntegerMatrix y = l["y"];
        NumericMatrix p = l["p"];
        NumericVector t = l["tau"];

        if (y.ncol() != loci) {
          throw std::range_error("Different number of loci (columns) in allele_range and y");
        }

        clusters[f] = y.nrow();
        ys[f] = y;
        ps[f] = p;
        taus[f] = t;

        f++;
      }

      for (int f = 0; f < n_fits; ++f) {
        hap_sum[f] = 0.0;
        match_within[f] = 0.0;
      }

      for (int f1 = 0; f1 < (n_fits-1); ++f1) {
        for (int f2 = (f1+1); f2 < n_fits; ++f2) {
          match_between[f1][f2] = 0.0;
        }
      }

      m_counters.resize(loci);
      m_start_alleles.resize(loci);
      m_stop_alleles.resize(loci);
      
      for (int l = 0; l < loci; l++) {
        m_start_alleles[l] = allele_range(0, l);
        m_stop_alleles[l] = allele_range(1, l);
      }

      //print(allele_range);
    }

    List calc() {
      haplotypes_queue.resize(CACHE_HAPLOTYPES_SUM_SIZE);
      haplotypes_queue_counter = 0;

      calc(0);
      do_work_with_haplotype_queue(); // empty queue if anything is left

      NumericVector res_hsum(n_fits);
      NumericVector res_within(n_fits);
      NumericMatrix res_between(n_fits, n_fits);

      for (int f = 0; f < n_fits; ++f) {
        res_hsum(f) = hap_sum[f];
        res_within(f) = match_within[f];
      }

      for (int f1 = 0; f1 < (n_fits-1); ++f1) {
        for (int f2 = (f1+1); f2 < n_fits; ++f2) {
          res_between(f1, f2) = match_between[f1][f2];
        }
      }

      List ret = List::create(
        _["hap_sum"] = res_hsum,
        _["match_within"] = res_within,
        _["match_between"] = res_between);

      return ret;
    }    
};

/*
 Normalise such that sum of haplotypes in db is 1
 */
// FIXME: mikl 2016-11-05: Is this working for only one locus?
class MatchProbSumsCacheOnlyDB {
  private:
    IntegerMatrix m_haplotypes;
    int n_haplotypes;

    int loci;

    int n_fits;
    std::vector<IntegerMatrix> ys;
    std::vector<NumericMatrix> ps;
    std::vector<NumericVector> taus;
    std::vector<int> clusters;

    std::vector<double> hap_sum;
    std::vector<double> match_within;
    std::vector< std::vector<double> > match_between;
    
    void calc_haplotypes() {
      std::vector< std::vector<double> > haplotypes_ps(n_fits);

      for (int f = 0; f < n_fits; ++f) {
        haplotypes_ps[f].resize(n_haplotypes);

        IntegerMatrix y_f = ys[f];
        NumericMatrix p_f = ps[f];
        NumericVector tau_f = taus[f];

        int clus = clusters[f];

        for (int h = 0; h < n_haplotypes; ++h) {
          IntegerVector htype_vec = m_haplotypes(h, Rcpp::_);
          std::vector<int> htype = as< std::vector<int> >(htype_vec);
          
          double hprob = 0.0;

          for (int c = 0; c < clus; c++) {
            IntegerVector yhap = y_f(c, _);
            double component_prob = tau_f(c);
            
            for (int l = 0; l < loci; l++) {
              double pcl = p_f(c, l);
              component_prob *= pow(pcl, abs(htype[l] - yhap(l)))*((1-pcl)/(1+pcl));
            }
            
            hprob += component_prob;
          }

          haplotypes_ps[f][h] = hprob;
          hap_sum[f] += hprob;
        }
      }      

      for (int h = 0; h < n_haplotypes; ++h) {
        std::vector<double> ph(n_fits);

        for (int f = 0; f < n_fits; ++f) {
          ph[f] = haplotypes_ps[f][h];
        }

        for (int f = 0; f < n_fits; ++f) {
          double norm_prob = ph[f]/hap_sum[f];
          match_within[f] += norm_prob * norm_prob;
        }

        for (int f1 = 0; f1 < (n_fits-1); ++f1) {
          double norm_prob1 = ph[f1]/hap_sum[f1];          
          for (int f2 = (f1+1); f2 < n_fits; ++f2) {
            double norm_prob2 = ph[f2]/hap_sum[f2];            
            match_between[f1][f2] += norm_prob1 * norm_prob2;
          }
        }
      }
    }

  public:
    MatchProbSumsCacheOnlyDB(const List& fits, const IntegerMatrix& haplotypes) { 
      m_haplotypes = haplotypes;
      loci = haplotypes.ncol();
      n_haplotypes = haplotypes.nrow();

      n_fits = fits.size();

      clusters.resize(n_fits);

      ys.resize(n_fits);
      ps.resize(n_fits);
      taus.resize(n_fits);

      hap_sum.resize(n_fits);
      match_within.resize(n_fits);
      match_between.resize(n_fits);
      
      for (int f = 0; f < n_fits; ++f) {
        match_between[f].resize(n_fits);
      }
      
      //print(fits);

      int f = 0;
      for (List::const_iterator it = fits.begin(); it != fits.end(); ++it ) {
        List l = *it;
        IntegerMatrix y = l["y"];
        NumericMatrix p = l["p"];
        NumericVector t = l["tau"];

        if (y.ncol() != loci) {
          throw std::range_error("Different number of loci (columns) in haplotypes and y");
        }

        clusters[f] = y.nrow();
        ys[f] = y;
        ps[f] = p;
        taus[f] = t;

        f++;
      }

      for (int f = 0; f < n_fits; ++f) {
        hap_sum[f] = 0.0;
        match_within[f] = 0.0;
      }      

      for (int f1 = 0; f1 < (n_fits-1); ++f1) {
        for (int f2 = (f1+1); f2 < n_fits; ++f2) {
          match_between[f1][f2] = 0.0;
        }
      }
    }

    List calc() {
      calc_haplotypes();
      
      NumericVector res_hsum(n_fits);
      NumericVector res_within(n_fits);
      NumericMatrix res_between(n_fits, n_fits);

      for (int f = 0; f < n_fits; ++f) {
        res_hsum(f) = hap_sum[f];
        res_within(f) = match_within[f];
      }

      for (int f1 = 0; f1 < (n_fits-1); ++f1) {
        for (int f2 = (f1+1); f2 < n_fits; ++f2) {
          res_between(f1, f2) = match_between[f1][f2];
        }
      }

      List ret = List::create(
        _["hap_sum"] = res_hsum,
        _["match_within"] = res_within,
        _["match_between"] = res_between);

      return ret;
    }    
};


// FIXME: mikl 2016-11-05: Is this working for only one locus?
class BinomialHapSums {
  private:  
    std::unordered_map< IntegerVector, std::vector<int>, std::hash<IntegerVector>, std::equal_to_intvec > hap_counts;

    int r;
    std::vector<int> n_db; 

    std::vector<double> hap_sum;
    std::vector<double> match_within;
    std::vector< std::vector<double> > match_between;
    
  public:
    BinomialHapSums(const List& dbs) {
      r = dbs.size();
      n_db.resize(r);
      
      if (r < 2) {
        stop("Expected at least two subpopulations");
      }
      
      hap_sum.resize(r);
      match_within.resize(r);
      match_between.resize(r);      
      for (int i = 0; i < r; ++i) {
        match_between[i].resize(r);
      }

      int i = 0;
      for( List::const_iterator it = dbs.begin(); it != dbs.end(); ++it ) {
        IntegerMatrix db = *it;
        
        int n = db.nrow();
        n_db[i] = n;
        
        for (int j = 0; j < n; ++j) {
          IntegerVector h = db(j, Rcpp::_);

          auto iter = hap_counts.find(h);
          
          if (iter != hap_counts.end()){
            iter->second[i] += 1;
            continue;
          } 

          // not found, add it
          std::vector<int> ns(r);
          ns[i] += 1;
          hap_counts[h] = ns;
        }

        ++i;    
      }
    }
    
    BinomialHapSums(const List& compact_dbs, const List& counts) {
      r = compact_dbs.size();
      n_db.resize(r);
      
      if (r < 2) {
        stop("Expected at least two subpopulations");
      }
      
      if (counts.size() != r) {
        stop("Expected as many counts as compact_dbs");
      }
      
      hap_sum.resize(r);
      match_within.resize(r);
      match_between.resize(r);      
      for (int i = 0; i < r; ++i) {
        match_between[i].resize(r);
      }
      
      int i = 0;
      for( List::const_iterator it = compact_dbs.begin(); it != compact_dbs.end(); ++it ) {
        IntegerMatrix db = *it;
        IntegerVector Ndb = counts[i];
        
        int n = db.nrow();
        
        if (Ndb.length() != n) {
          stop("Expected as many Ns as haplotypes");
        }
        
        n_db[i] = sum(Ndb);
        
        for (int j = 0; j < n; ++j) {
          IntegerVector h = db(j, Rcpp::_);
          
          auto iter = hap_counts.find(h);
          
          if (iter != hap_counts.end()){
            iter->second[i] += Ndb[j];
            continue;
          } 
          
          // not found, add it
          std::vector<int> ns(r);
          ns[i] += Ndb[j]; // the rest is 0
          hap_counts[h] = ns;
        }
        
        ++i;    
      }
    }
    
    List calc() {
      auto ns = hap_counts.begin();
      
      while (ns != hap_counts.end()){
        std::vector<double> ph(r);
        std::vector<int> ns_vec = ns->second;

        for (int i = 0; i < r; ++i) {          
          ph[i] = (double)(ns_vec[i]) / (double)(n_db[i]);
        }

        for (int i = 0; i < r; ++i) {
          hap_sum[i] += ph[i];

          int n_tmp = ns_vec[i];

          if (n_tmp <= 1) { // in this case (n_tmp - 1) / (n_db - 1) is zero and would not contribute to the sum
            continue;
          }

          match_within[i] += ph[i] * ( (double)(n_tmp - 1) / (double)(n_db[i] - 1) );
        }

        for (int i1 = 0; i1 < (r-1); ++i1) {
          for (int i2 = (i1+1); i2 < r; ++i2) {
            match_between[i1][i2] += ph[i1] * ph[i2];
          }
        }

        ++ns;
      }

      NumericVector res_hsum(r);
      NumericVector res_within(r);
      NumericMatrix res_between(r, r);

      for (int i = 0; i < r; ++i) {
        res_hsum(i) = hap_sum[i];
        res_within(i) = match_within[i];
      }

      for (int i1 = 0; i1 < (r-1); ++i1) {
        for (int i2 = (i1+1); i2 < r; ++i2) {
          res_between(i1, i2) = match_between[i1][i2];
        }
      }

      List ret = List::create(
        _["hap_sum"] = res_hsum,
        _["match_within"] = res_within,
        _["match_between"] = res_between);

      return ret;
    }    
};


// FIXME: mikl 2016-11-05: This works for one locus!
class PedigreeSums {
  private:  
    std::unordered_map< int, std::vector<int> > hap_counts;

    int r;
    std::vector<int> n_db; 

    std::vector<double> hap_sum;
    std::vector<double> match_within;
    std::vector< std::vector<double> > match_between;
    
  public:
    PedigreeSums(const List& compact_dbs, const List& counts) {
      // mikl 2016-11-04 ->
      Rcpp::Rcout << "stack size being increased..." << std::endl;
      
      const rlim_t kStackSize = 256 * 1024 * 1024;   // min stack size = 256 MB
      struct rlimit rl;
      int result;

      result = getrlimit(RLIMIT_STACK, &rl);
      if (result == 0) {
        if (rl.rlim_cur < kStackSize)
        {
          rl.rlim_cur = kStackSize;
          result = setrlimit(RLIMIT_STACK, &rl);
          if (result != 0) {
            Rcpp::Rcout << "setrlimit returned result = " << result << std::endl;
          }
        }
      }      
      // <- mikl 2016-11-04
      
      r = compact_dbs.size();
      n_db.resize(r);
      
      if (r < 2) {
        stop("Expected at least two subpopulations");
      }
      
      if (counts.size() != r) {
        stop("Expected as many counts as compact_dbs");
      }
      
      // mikl 2016-11-04
      Progress progress(compact_dbs.size(), true);
      
      hap_sum.resize(r);
      match_within.resize(r);
      match_between.resize(r);      
      for (int i = 0; i < r; ++i) {
        match_between[i].resize(r);
      }
      
      int i = 0;
      for(List::const_iterator it = compact_dbs.begin(); it != compact_dbs.end(); ++it ) {
        IntegerMatrix db = *it;
        
        if (db.ncol() != 1) {
          stop("Only expected one column (the pedigree id)");
        }
        
        IntegerVector Ndb = counts[i];
        
        int n = db.nrow();
        
        if (Ndb.length() != n) {
          stop("Expected as many Ns as haplotypes");
        }
        
        n_db[i] = sum(Ndb);
        
        for (int j = 0; j < n; ++j) {
          int pedid = db(j, 0);
          
          auto iter = hap_counts.find(pedid);
          
          if (iter != hap_counts.end()){
            iter->second[i] += Ndb[j];
            continue;
          } 
          
          // not found, add it
          std::vector<int> ns(r);
          ns[i] += Ndb[j]; // the rest is 0
          hap_counts[pedid] = ns;
        }
        
        ++i;    
        
        // mikl 2016-11-04 ->
        if (Progress::check_abort() ) {
          stop("Aborted...");
        }
        progress.increment();
        // <- mikl 2016-11-04
      }
    }
    
    List calc() {
      // mikl 2016-11-04 ->
      Rcpp::Rcout << "hap_counts size = " << hap_counts.size() << std::endl;
      Progress progress(hap_counts.size(), true);
      // <- mikl 2016-11-04
      
      auto ns = hap_counts.begin();
      
      while (ns != hap_counts.end()) {
        std::vector<double> ph(r);
        std::vector<int> ns_vec = ns->second;

        for (int i = 0; i < r; ++i) {          
          ph[i] = (double)(ns_vec[i]) / (double)(n_db[i]);
        }

        for (int i = 0; i < r; ++i) {
          hap_sum[i] += ph[i];

          int n_tmp = ns_vec[i];

          if (n_tmp <= 1) { // in this case (n_tmp - 1) / (n_db - 1) is zero and would not contribute to the sum
            continue;
          }

          match_within[i] += ph[i] * ( (double)(n_tmp - 1) / (double)(n_db[i] - 1) );
        }

        for (int i1 = 0; i1 < (r-1); ++i1) {
          for (int i2 = (i1+1); i2 < r; ++i2) {
            match_between[i1][i2] += ph[i1] * ph[i2];
          }
        }

        ++ns;
        
        // mikl 2016-11-04 ->
        if (Progress::check_abort() ) {
          stop("Aborted...");
        }
        progress.increment();
        // <- mikl 2016-11-04
      }

      NumericVector res_hsum(r);
      NumericVector res_within(r);
      NumericMatrix res_between(r, r);

      for (int i = 0; i < r; ++i) {
        res_hsum(i) = hap_sum[i];
        res_within(i) = match_within[i];
      }

      for (int i1 = 0; i1 < (r-1); ++i1) {
        for (int i2 = (i1+1); i2 < r; ++i2) {
          res_between(i1, i2) = match_between[i1][i2];
        }
      }

      List ret = List::create(
        _["hap_sum"] = res_hsum,
        _["match_within"] = res_within,
        _["match_between"] = res_between);

      return ret;
    }    
};


void nested_loop_operation_haplotype_probabilities_sum(double* res[], int counters[], int start_alleles[], int stop_alleles[], int level, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int loci = y.ncol();
  
  if (level == loci) {
    int clusters = y.nrow();

    /*
    Rcout << "(";
    for (int l = 0; l < loci; l++) {
      Rcout << counters[l] << ", ";
    }
    Rcout << ")" << std::endl;
    */
    
    double hprob = 0.0;
        
    for (int c = 0; c < clusters; c++) {
      IntegerVector yhap = y(c, _);
      double component_prob = tau(c);
      
      for (int l = 0; l < loci; l++) {
        double pcl = p(c, l);
        component_prob *= pow(pcl, abs(counters[l] - yhap(l)))*((1-pcl)/(1+pcl));
      }
      
      hprob += component_prob;
    }
    
    //Rcout << (*res)[0] << std::endl;
    
    (*res)[0] += hprob;
    (*res)[1] += hprob * hprob;
  } else {
    for (int allele = start_alleles[level]; allele <= stop_alleles[level]; allele++) {
      counters[level] = allele;
      nested_loop_operation_haplotype_probabilities_sum(res, counters, start_alleles, stop_alleles, level + 1, y, p, tau);
    }
  }
}

// [[Rcpp::export]]
List rcpp_calculate_haplotype_probabilities_sum_CLASS(List fits, IntegerMatrix allele_range) {
  MatchProbSums s = MatchProbSums(fits, allele_range);
  return s.calc();
}

// [[Rcpp::export]]
List rcpp_calculate_haplotype_probabilities_sum_CLASS_Cache(List fits, IntegerMatrix allele_range) {
  MatchProbSumsCache s = MatchProbSumsCache(fits, allele_range);
  return s.calc();
}

// [[Rcpp::export]]
List rcpp_calculate_haplotype_probabilities_sum_binomial(List dbs) {
  BinomialHapSums s = BinomialHapSums(dbs);
  return s.calc();
}

// [[Rcpp::export]]
List rcpp_calculate_haplotype_probabilities_sum_binomial_compact_dbs(List compact_dbs, List counts) {
  BinomialHapSums s = BinomialHapSums(compact_dbs, counts);
  return s.calc();
}

// [[Rcpp::export]]
List rcpp_calculate_haplotype_probabilities_sum_pedigrees_compact_dbs(List compact_dbs, List counts) {
  PedigreeSums s = PedigreeSums(compact_dbs, counts);
  return s.calc();
}


// [[Rcpp::export]]
List rcpp_hapsums_disclap_normalised(List fits, IntegerMatrix haplotypes) {
  MatchProbSumsCacheOnlyDB s = MatchProbSumsCacheOnlyDB(fits, haplotypes);
  return s.calc();
}


// [[Rcpp::export]]
NumericVector rcpp_calculate_haplotype_probabilities_sum(IntegerMatrix allele_range, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int loci = allele_range.ncol();

  if (y.ncol() != loci) {
    throw std::range_error("Different number of loci (columns) in allele_range and y");
  }

  if (allele_range.nrow() != 2) {
    throw std::range_error("Exactly two rows in allele_range required");
  }
  
  double* res = new double[2];
  int* counters = new int[loci];
  int* start_alleles = new int[loci];
  int* stop_alleles = new int[loci];
  
  res[0] = 0.0;
  res[1] = 0.0;
  
  for (int l = 0; l < loci; l++) {
    start_alleles[l] = allele_range(0, l);
    stop_alleles[l] = allele_range(1, l);
  }
  
  nested_loop_operation_haplotype_probabilities_sum(&res, counters, start_alleles, stop_alleles, 0, y, p, tau);
  
  // Need to convert to vector, or else the return is not OK. Weird.
  std::vector<double> ans(2);
  ans[0] = res[0];
  ans[1] = res[1];
  
  delete[] res;
  delete[] counters;
  delete[] start_alleles;
  delete[] stop_alleles;
  
  NumericVector ret = NumericVector::create(
    _["sum_prob"] = ans[0],
    _["sum_sq_prob"] = ans[1]);
  
  return ret;
}


// [[Rcpp::export]]
NumericVector rcpp_calculate_haplotype_probabilities_sum_between(IntegerMatrix allele_range, 
  IntegerMatrix y1, IntegerMatrix y2, 
  NumericMatrix p1, NumericMatrix p2, 
  NumericVector tau1, NumericVector tau2) {
  

  NumericVector ret = NumericVector::create(
    _["sum_prob_1"] = 0.0,
    _["sum_prob_2"] = 0.0,
    _["sum_prob_prod"] = 0.0);
  
  return ret;
}

/****************************/

void nested_loop_operation_match_quantities(NumericMatrix res, int counters[], int start_alleles[], int stop_alleles[], int level, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int loci = y.ncol();
  
  if (level == loci) {
    int clusters = y.nrow();
    
    NumericVector happrobs(clusters);
    
    for (int c = 0; c < clusters; c++) {
      IntegerVector yhap = y(c, _);
      double component_prob = 1.0;
      
      for (int l = 0; l < loci; l++) {
        double pcl = p(c, l);
        component_prob *= pow(pcl, abs(counters[l] - yhap(l)))*((1-pcl)/(1+pcl));
      }
      
      happrobs(c) = component_prob;
    }
    
    for (int c1 = 0; c1 < (clusters - 1); c1++) {
      for (int c2 = c1+1; c2 < clusters; c2++) {
        res(c1, c2) += happrobs(c1) * happrobs(c2);
      }
    }
  } else {
    for (int allele = start_alleles[level]; allele <= stop_alleles[level]; allele++) {
      counters[level] = allele;
      nested_loop_operation_match_quantities(res, counters, start_alleles, stop_alleles, level + 1, y, p, tau);
    }
  }
}

// [[Rcpp::export]]
NumericVector rcpp_match_quantities(IntegerMatrix allele_range, IntegerMatrix y, NumericMatrix p, NumericVector tau) {
  int clusters = y.nrow();
  int loci = allele_range.ncol();

  if (y.ncol() != loci) {
    throw std::range_error("Different number of loci (columns) in x and y");
  }

  if (allele_range.nrow() != 2) {
    throw std::range_error("Exactly two rows in allele_range required");
  }
  
  NumericMatrix res(clusters, clusters);
  int* counters = new int[loci];
  int* start_alleles = new int[loci];
  int* stop_alleles = new int[loci];
  
  for (int l = 0; l < loci; l++) {
    start_alleles[l] = allele_range(0, l);
    stop_alleles[l] = allele_range(1, l);
  }
  
  nested_loop_operation_match_quantities(res, counters, start_alleles, stop_alleles, 0, y, p, tau);
  
  delete[] counters;
  delete[] start_alleles;
  delete[] stop_alleles;
  
  return res;
}

