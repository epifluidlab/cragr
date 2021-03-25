#include <Rcpp.h>
#include <string.h>
#include <assert.h>
#include <math.h>
using namespace Rcpp;
using namespace std;


NumericVector rollapply(NumericVector x, int k, double (*func)(NumericVector, int, int, bool), bool na_pad = false, bool na_rm = false, std::string align = "left") {
  int mode;
  if (strcmp(align.c_str(), "left") == 0)
    mode = 1;
  else if (strcmp(align.c_str(), "center") == 0)
    mode = 2;
  else if (strcmp(align.c_str(), "right") == 0)
    mode = 3;
  else
    mode = 0;
  assert(mode > 0);

  // Initialize output vector
  int n = x.size() - k + 1;
  if (n <= 0) {
    return NULL;
  }
  if (k >= x.size()) return NULL;

  int out_len = n;
  if (na_pad) {
    out_len = x.size();
  }

  NumericVector out(out_len);

  // Based on mode, determine the start and end index
  int start = -1;
  int end = -1;
  int shift_l = -1;
  int shift_r = -1;

  if (mode == 1) {
    start = 0;
    end = n;
    shift_l = 0;
    shift_r = k - 1;
  } else if (mode == 2) {
    start = (int)floor((k - 1) / 2);
    end = start + n;
    shift_l = start;
    shift_r = shift_l + k - 1;
  } else if (mode == 3) {
    start = k - 1;
    end = start + n;
    shift_l = start;
    shift_r = shift_l + k - 1;
  }

  for (int i = start; i < end; ++i) {
    double v = func(x, i - shift_l, k, na_rm);
    // double total = 0;
    // for (int j = i - shift_l; j < i - shift_l + k; ++j) {
    //   if (NumericVector::is_na(x[j])) {
    //     if (na_rm) {
    //       continue;
    //     } else {
    //       total = NA_REAL;
    //       break;
    //     }
    //   }
    //   total += x[j];
    // }
    out[i] = v;
  }

  if (na_pad) {
    for (int i = 0; i < start; ++i) {
      out[i] = NA_REAL;
    }
    for (int i = end; i < out_len; ++i) {
      out[i] = NA_REAL;
    }
  }

  return out;
}


double calc_mean(NumericVector x, int start, int length, bool na_rm) {
  double total = 0;
  int count = 0;
  for (int i = start; i < start + length; ++i) {
    if (NumericVector::is_na(x[i])) {
      if (na_rm) {
        continue;
      } else {
        total = NA_REAL;
        break;
      }
    }

    count += 1;
    total += x[i];
  }

  if (NumericVector::is_na(total)) {
    return NA_REAL;
  } else{
    return total / count;
  }
}


double calc_sum(NumericVector x, int start, int length, bool na_rm) {
  double total = 0;
  for (int i = start; i < start + length; ++i) {
    if (NumericVector::is_na(x[i])) {
      if (na_rm) {
        continue;
      } else {
        total = NA_REAL;
        break;
      }
    }
    total += x[i];
  }
  return total;
}


// [[Rcpp::export]]
NumericVector rollmean(NumericVector x, int k, bool na_pad = false, bool na_rm = false, std::string align = "left") {
  return rollapply(x, k, &calc_mean, na_pad, na_rm, align);
}


// [[Rcpp::export]]
NumericVector rollsum(NumericVector x, int k, bool na_pad = false, bool na_rm = false, std::string align = "left") {
  return rollapply(x, k, &calc_sum, na_pad, na_rm, align);
}

