/* Copyright 2021 Lorin Crawford.
 */
#include "mapit/util.h"

context("num_combinations_with_replacement") {
  test_that("num_combinations_with_replacement(2, 2) = 3") {
    // given
    int num_available = 2;
    int num_selected = 2;
    int correct_answer = 3;
    // when
    int result = num_combinations_with_replacement(num_available, num_selected);
    // then
    expect_true(result == correct_answer);
  }
  test_that("num_combinations_with_replacement(3, 2) = 6") {
    // given
    int num_available = 3;
    int num_selected = 2;
    int correct_answer = 6;
    // when
    int result = num_combinations_with_replacement(num_available, num_selected);
    // then
    expect_true(result == correct_answer);
  }
  test_that("num_combinations_with_replacement(1, 2) = 0") {
    // given
    int num_available = 1;
    int num_selected = 2;
    int correct_answer = 1;
    // when
    int result = num_combinations_with_replacement(num_available, num_selected);
    // then
    expect_true(result == correct_answer);
  }
  test_that("num_combinations_with_replacement(1, 1) = 1") {
    // given
    int num_available = 1;
    int num_selected = 1;
    int correct_answer = 1;
    // when
    int result = num_combinations_with_replacement(num_available, num_selected);
    // then
    expect_true(result == correct_answer);
  }
}

context("factorial") {
  test_that("factorial(3) = 6") {
    // given
    int n = 3;
    int correct_answer = 6;
    // when
    int result = factorial(n);
    // then
    expect_true(result == correct_answer);
  }
  test_that("factorial(4) = 24") {
    // given
    int n = 4;
    int correct_answer = 24;
    // when
    int result = factorial(n);
    // then
    expect_true(result == correct_answer);
  }
  test_that("factorial(5) = 120") {
    // given
    int n = 5;
    int correct_answer = 120;
    // when
    int result = factorial(n);
    // then
    expect_true(result == correct_answer);
  }
}

context("matrix_to_vector_of_rows") {
  test_that("matrix_to_vector_of_rows returns correct vectors") {
    // given
    arma::mat matrix(3, 3);
    matrix.eye();
    arma::vec v1(3, arma::fill::zeros);
    v1(0) = 1;
    arma::vec v2(3, arma::fill::zeros);
    v2(1) = 1;
    arma::vec v3(3, arma::fill::zeros);
    v3(2) = 1;
    // when
    std::vector<arma::vec> vectors = matrix_to_vector_of_rows(matrix);
    // then
    expect_true(arma::approx_equal(vectors[0], v1, "absdiff", 0.01));
    expect_true(arma::approx_equal(vectors[1], v2, "absdiff", 0.01));
    expect_true(arma::approx_equal(vectors[2], v3, "absdiff", 0.01));
  }
  test_that("matrix_to_vector_of_rows works with nx1 matrix") {
    // given
    arma::mat matrix(1, 3);
    matrix.zeros();
    matrix(0, 0) = 1;
    arma::vec correct_answer(3, arma::fill::zeros);
    correct_answer(0) = 1;
    // when
    std::vector<arma::vec> result = matrix_to_vector_of_rows(matrix);
    // then
    expect_true(result[0].size() == 3);
    expect_true(typeid(result[0]).name() == typeid(correct_answer).name());
    expect_true(arma::approx_equal(result[0], correct_answer, "absdiff", 0.01));
  }
}

context("vectorise_to_matrix") {
  test_that("vectorise_to_matrix returns correct vectors") {
    // given
    arma::mat matrix(3, 3);
    matrix.eye();
    arma::mat correct_answer(9, 1);
    correct_answer.zeros();
    correct_answer(0, 0) = 1;
    correct_answer(4, 0) = 1;
    correct_answer(8, 0) = 1;
    // when
    arma::mat result = vectorise_to_matrix(matrix);
    // then
    expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.01));
  }
}

context("skip_variant") {
  test_that("test_variant does not skip when empty") {
    // given
    int i = 1;
    arma::vec ind;
    // when
    bool skip = skip_variant(ind, i);
    // then
    expect_true(skip == false);
  }
  test_that("test_variant skips when not in list") {
    // given
    int i = 10;
    arma::vec ind(1);
    ind = 2;
    // when
    bool skip = skip_variant(ind, i);
    // then
    expect_true(skip == true);
  }
  test_that("test_variant does not skip when in list") {
    // given
    int i = 1;
    arma::vec ind(1);
    ind = i + 1;
    // when
    bool skip = skip_variant(ind, i);
    // then
    expect_true(skip == false);
  }
}

context("compute_principal_components") {
  // Problem from:
  // https://mysite.science.uottawa.ca/phofstra/MAT2342/SVDproblems.pdf However
  // it has an error, U has wrong signs
  test_that("compute_principal_components returns right values") {
    // given
    arma::mat X = {{0, 1, 1}, {sqrt(2), 2, 0}, {0, 1, 1}};
    arma::mat answer = {{-2 / sqrt(3), sqrt(2) / sqrt(3), 0},
                        {-4 / sqrt(3), -sqrt(2) / sqrt(3), 0},
                        {-2 / sqrt(3), sqrt(2) / sqrt(3), 0}};
    // when
    arma::mat pc_cols = compute_principal_components(X, 3);
    // then
    expect_true(pc_cols.n_cols == 3);
    expect_true(arma::approx_equal(pc_cols, answer, "absdiff", 0.01));
  }
}

context("index_combinations") {
  test_that("index_combinations(2) = 00, 01, 11") {
    // given
    int num_available = 2;
    arma::mat correct_answer = {{0, 0}, {0, 1}, {1, 1}};
    // when
    arma::mat result = index_combinations(num_available);
    // then
    expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.01));
  }
  test_that("index_combinations(3) = 00, 01, 11, 02, 12, 22") {
    // given
    int num_available = 3;
    arma::mat correct_answer = {{0, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 2}, {2, 2}};
    // when
    arma::mat result = index_combinations(num_available);
    // then
    expect_true(arma::approx_equal(result, correct_answer, "absdiff", 0.01));
  }
}

context("find_row_vector") {
  test_that("find_row_vector((x,y), mat) = i") {
    // given
    arma::rowvec v1 = {0, 0};
    arma::rowvec v2 = {1, 1};
    arma::rowvec v3 = {0, 2};
    arma::rowvec v4 = {2, 2};
    arma::mat indices = {{0, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 2}, {2, 2}};
    int correct_answer_v1 = 0;
    int correct_answer_v2 = 2;
    int correct_answer_v3 = 3;
    int correct_answer_v4 = 5;
    // when
    int result_v1 = find_row_vector(v1, indices);
    int result_v2 = find_row_vector(v2, indices);
    int result_v3 = find_row_vector(v3, indices);
    int result_v4 = find_row_vector(v4, indices);
    // then
    expect_true(result_v1 == correct_answer_v1);
    expect_true(result_v2 == correct_answer_v2);
    expect_true(result_v3 == correct_answer_v3);
    expect_true(result_v4 == correct_answer_v4);
  }
  test_that("find_row_vector throws exception when vector not exist") {
    // given
    arma::rowvec v = {2, 3};
    arma::mat indices = {{0, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 2}, {2, 2}};
    // when
    try {
      int result = find_row_vector(v, indices);
    } catch (const char *msg) {
      // then
      std::string str(msg);
      expect_true(str.compare("Row vector not found.") == 0);
    }
  }
  test_that("if multiple find first") {
    // given
    arma::rowvec v = {1, 1};
    arma::mat indices = {{0, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 1}, {2, 2}};
    int correct_answer = 2;
    // when
    int result = find_row_vector(v, indices);
    // then
    expect_true(result == correct_answer);
  }
}
