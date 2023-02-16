#include <iostream>
#include <cmath>
#include <unordered_set>
#include <string.h>

/*
 * test harness adapted from
 * https://github.com/noporpoise/BitArray/blob/master/dev/bit_array_test.c
 */
const char *suite_name;
char suite_pass;
int suites_run = 0, suites_failed = 0, suites_empty = 0;
int tests_in_suite = 0, tests_run = 0, tests_failed = 0;

#define QUOTE(str) #str
#define ASSERT(x) {tests_run++; tests_in_suite++; if(!(x)) \
    { fprintf(stderr, "failed assert [%s:%i] %s\n", __FILE__, __LINE__, QUOTE(x)); \
      suite_pass = 0; tests_failed++; }}

inline void SUITE_START(const char *name) {
  suite_pass = 1;
  suite_name = name;
  suites_run++;
  tests_in_suite = 0;
}

inline void SUITE_END() {
  printf("Testing %s ", suite_name);
  size_t suite_i;
  for(suite_i = strlen(suite_name); suite_i < 80-8-5; suite_i++) printf(".");
  printf("%s\n", suite_pass ? " pass" : " fail");
  if(!suite_pass) suites_failed++;
  if(!tests_in_suite) suites_empty++;
}
/*
 *  End adapted code
 */

inline std::vector<bool> _bool_array_to_vector(bool *arr, unsigned int n) {
    std::vector<bool> vec;

    for(unsigned int i = 0; i < n; i++)
        vec.push_back(arr[i]);

    return vec;
}

inline std::vector<uint32_t> _uint32_array_to_vector(uint32_t *arr, unsigned int n) {
    std::vector<uint32_t> vec;

    for(unsigned int i = 0; i < n; i++)
        vec.push_back(arr[i]);

    return vec;
}

inline std::vector<double> _double_array_to_vector(double *arr, unsigned int n) {
    std::vector<double> vec;

    for(unsigned int i = 0; i < n; i++)
        vec.push_back(arr[i]);

    return vec;
}

inline std::vector<std::string> _string_array_to_vector(std::string *arr, unsigned int n) {
    std::vector<std::string> vec;

    for(unsigned int i = 0; i < n; i++)
        vec.push_back(arr[i]);

    return vec;
}

inline bool vec_almost_equal(std::vector<double> a, std::vector<double> b) {
    if(a.size() != b.size()) {
        return false;
    }
    for(unsigned int i = 0; i < a.size(); i++) {
        if(!(fabs(a[i] - b[i]) < 0.000001)) {  // sufficient given the tests
            return false;
        }
    }
    return true;
}