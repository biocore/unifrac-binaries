#include "api.hpp"

#include <signal.h>
#include <thread>
#include <vector>

#ifndef API_ONLY
#include "tree.hpp"
#include "biom.hpp"
#include "unifrac.hpp"
#include "unifrac_internal.hpp"
#endif

#include "test_helper.hpp"

// copy of internal function in api... repeated here for testing
uint64_t _testv_comb_2(uint64_t N) {
            // based off of _comb_int_long
            // https://github.com/scipy/scipy/blob/v0.19.1/scipy/special/_comb.pyx
            
            // Compute binom(N, k) for integers.
            //
            // we're disregarding overflow as that practically should not
            // happen unless the number of samples processed is in excess
            // of 4 billion 
            uint64_t val, j, M, nterms;
            uint64_t k = 2;

            M = N + 1;
            nterms = k < (N - k) ? k : N - k;

            val = 1;

            for(j = 1; j < nterms + 1; j++) {
                val *= M - j;
                val /= j;
            }
            return val;
}

// minor variant of internal function in api... repeated here for testing
void _testv_stripes_to_condensed_form(double* stripes[], uint32_t n, uint32_t m, double* cf) {
    uint64_t comb_N = _testv_comb_2(n);
    for(unsigned int stripe = 0; stripe < m; stripe++) {
        // compute the (i, j) position of each element in each stripe
        uint64_t i = 0;
        uint64_t j = stripe + 1;
        for(uint64_t k = 0; k < n; k++, i++, j++) {
            if(j == n) {
                i = 0;
                j = n - (stripe + 1);
            }
            // determine the position in the condensed form vector for a given (i, j)
            // based off of
            // https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
            uint64_t comb_N_minus_i = _testv_comb_2(n - i);
            cf[comb_N - comb_N_minus_i + (j - i - 1)] = stripes[stripe][k];
        }
    }
}


#ifndef API_ONLY
void test_bptree_simple_result(const su::BPTree &tree) {
    unsigned int exp_nparens = 8;

    std::vector<bool> exp_structure;
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(false);
    exp_structure.push_back(false);

    std::vector<uint32_t> exp_openclose;
    exp_openclose.push_back(7);
    exp_openclose.push_back(6);
    exp_openclose.push_back(3);
    exp_openclose.push_back(2);
    exp_openclose.push_back(5);
    exp_openclose.push_back(4);
    exp_openclose.push_back(1);
    exp_openclose.push_back(0);

    std::vector<std::string> exp_names;
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("c"));
    exp_names.push_back(std::string("123:foo; bar"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("b"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());

    std::vector<double> exp_lengths;
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(1.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(2.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(tree.get_openclose() == exp_openclose);
    ASSERT(tree.lengths == exp_lengths);
    ASSERT(tree.names == exp_names);
}

void test_bptree_constructor_simple() {
    SUITE_START("bptree constructor simple");
                                //01234567
                                //11101000
    su::BPTree tree("(('123:foo; bar':1,b:2)c);");

    test_bptree_simple_result(tree);

    SUITE_END();
}

void test_bptree_constructor_simple_cpp() {
    SUITE_START("bptree constructor simple cpp");
                                //01234567
                                //11101000
    su::BPTree tree(std::string("(('123:foo; bar':1,b:2)c);"));

    test_bptree_simple_result(tree);

    SUITE_END();
}

void test_bptree_constructor_from_existing() {
    SUITE_START("bptree constructor from_existing");
                                //01234567
                                //11101000
    su::BPTree existing("(('123:foo; bar':1,b:2)c);");
    su::BPTree tree(existing.get_structure(), existing.lengths, existing.names);
 
    test_bptree_simple_result(tree);

    SUITE_END();
}

void test_bptree_constructor_inmem() {
    SUITE_START("bptree constructor inmem");
                                //01234567
                                //11101000
    su::BPTree tree_org("(('123:foo; bar':1,b:2)c);");
    test_bptree_simple_result(tree_org);


    bool t_structure[] = { true, true, true, false, true, false, false, false };
    double t_lengths[] = { 0, 0, 1, 0, 2, 0, 0, 0 };
    const char*  t_names[] = { "", "c", "123:foo; bar", "", "b", "", "", "" };
    const support_bptree_t support_tree = {t_structure, t_lengths, (char**) t_names, 8};

    su::BPTree tree(support_tree.structure,
                    support_tree.lengths,
                    support_tree.names,
                    support_tree.n_parens);
    test_bptree_simple_result(tree);

    SUITE_END();
}

void test_bptree_mask() {
    SUITE_START("bptree mask");
                                //01234567
                                //11101000
                                //111000
    std::vector<bool> mask = {true, true, true, true, false, false, true, true};
    su::BPTree base("(('123:foo; bar':1,b:2)c);");
    su::BPTree tree = base.mask(mask, base.lengths);
    unsigned int exp_nparens = 6;

    std::vector<bool> exp_structure;
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(false);
    exp_structure.push_back(false);

    std::vector<uint32_t> exp_openclose;
    exp_openclose.push_back(5);
    exp_openclose.push_back(4);
    exp_openclose.push_back(3);
    exp_openclose.push_back(2);
    exp_openclose.push_back(1);
    exp_openclose.push_back(0);

    std::vector<std::string> exp_names;
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("c"));
    exp_names.push_back(std::string("123:foo; bar"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());

    std::vector<double> exp_lengths;
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(1.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(tree.get_openclose() == exp_openclose);
    ASSERT(tree.lengths == exp_lengths);
    ASSERT(tree.names == exp_names);

    SUITE_END();
}

void test_bptree_constructor_single_descendent() {
    SUITE_START("bptree constructor single descendent");

    su::BPTree tree("(((a)b)c,((d)e)f,g)r;");

    unsigned int exp_nparens = 16;

    bool structure_arr[] = {1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0};
    std::vector<bool> exp_structure = _bool_array_to_vector(structure_arr, exp_nparens);

    double length_arr[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> exp_lengths = _double_array_to_vector(length_arr, exp_nparens);

    std::string names_arr[] = {"r", "c", "b", "a", "", "", "", "f", "e", "d", "", "", "", "g", "", ""};
    std::vector<std::string> exp_names = _string_array_to_vector(names_arr, exp_nparens);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(vec_almost_equal(tree.lengths, exp_lengths));
    ASSERT(tree.names == exp_names);

    SUITE_END();
}

void test_bptree_constructor_complex() {
    SUITE_START("bp tree constructor complex");
    su::BPTree tree("(((a:1,b:2.5)c:6,d:8,(e),(f,g,(h:1,i:2)j:1)k:1.2)l,m:2)r;");

    unsigned int exp_nparens = 30;

    bool structure_arr[] = {1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0};
    std::vector<bool> exp_structure = _bool_array_to_vector(structure_arr, exp_nparens);

    double length_arr[] = {0, 0, 6, 1, 0, 2.5, 0, 0, 8, 0, 0, 0, 0, 0, 1.2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 2, 0, 0};
    std::vector<double> exp_lengths = _double_array_to_vector(length_arr, exp_nparens);

    std::string names_arr[] = {"r", "l", "c", "a", "", "b", "", "", "d", "", "", "e", "", "", "k", "f", "", "g", "", "j", "h", "", "i", "", "", "", "", "m", "", ""};
    std::vector<std::string> exp_names = _string_array_to_vector(names_arr, exp_nparens);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(vec_almost_equal(tree.lengths, exp_lengths));
    ASSERT(tree.names == exp_names);
    SUITE_END();
}

void test_bptree_constructor_semicolon() {
    SUITE_START("bp tree constructor semicolon");
    su::BPTree tree("((a,(b,c):5)'d','e; foo':10,((f))g)r;");

    unsigned int exp_nparens = 20;

    bool structure_arr[] = {1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0};
    std::vector<bool> exp_structure = _bool_array_to_vector(structure_arr, exp_nparens);

    double length_arr[] = {0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> exp_lengths = _double_array_to_vector(length_arr, exp_nparens);

    std::string names_arr[] = {"r", "d", "a", "", "", "b", "", "c", "", "", "", "e; foo", "", "g", "", "f", "", "", "", ""};
    std::vector<std::string> exp_names = _string_array_to_vector(names_arr, exp_nparens);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(vec_almost_equal(tree.lengths, exp_lengths));
    ASSERT(tree.names == exp_names);
    SUITE_END();
}

void test_bptree_constructor_edgecases() {
    SUITE_START("bp tree constructor edgecases");

    su::BPTree tree1("((a,b));");
    bool structure_arr1[] = {1, 1, 1, 0, 1, 0, 0, 0};
    std::vector<bool> exp_structure1 = _bool_array_to_vector(structure_arr1, 8);

    su::BPTree tree2("(a);");
    bool structure_arr2[] = {1, 1, 0, 0};
    std::vector<bool> exp_structure2 = _bool_array_to_vector(structure_arr2, 4);

    su::BPTree tree3("();");
    bool structure_arr3[] = {1, 1, 0, 0};
    std::vector<bool> exp_structure3 = _bool_array_to_vector(structure_arr3, 4);

    su::BPTree tree4("((a,b),c);");
    bool structure_arr4[] = {1, 1, 1, 0, 1, 0, 0, 1, 0, 0};
    std::vector<bool> exp_structure4 = _bool_array_to_vector(structure_arr4, 10);

    su::BPTree tree5("(a,(b,c));");
    bool structure_arr5[] = {1, 1, 0, 1, 1, 0, 1, 0, 0, 0};
    std::vector<bool> exp_structure5 = _bool_array_to_vector(structure_arr5, 10);

    ASSERT(tree1.get_structure() == exp_structure1);
    ASSERT(tree2.get_structure() == exp_structure2);
    ASSERT(tree3.get_structure() == exp_structure3);
    ASSERT(tree4.get_structure() == exp_structure4);
    ASSERT(tree5.get_structure() == exp_structure5);

    SUITE_END();
}

void test_bptree_constructor_quoted_comma() {
    SUITE_START("quoted comma bug");
    su::BPTree tree("((3,'foo,bar')x,c)r;");
    std::vector<std::string> exp_names = {"r", "x", "3", "", "foo,bar", "", "", "c", "", ""};
    ASSERT(exp_names.size() == tree.names.size());

    for(unsigned int i = 0; i < tree.names.size(); i++) {
        ASSERT(exp_names[i] == tree.names[i]);
    }
    SUITE_END();
}

void test_bptree_constructor_quoted_parens() {
    SUITE_START("quoted parens");
    su::BPTree tree("((3,'foo(b)ar')x,c)r;");
    std::vector<std::string> exp_names = {"r", "x", "3", "", "foo(b)ar", "", "", "c", "", ""};
    ASSERT(exp_names.size() == tree.names.size());

    for(unsigned int i = 0; i < tree.names.size(); i++) {
        ASSERT(exp_names[i] == tree.names[i]);
    }
    SUITE_END();
}
void test_bptree_postorder() {
    SUITE_START("postorderselect");

    // fig1 from https://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    su::BPTree tree("((3,4,(6)5)2,7,((10,100)9)8)1;");
    uint32_t exp[] = {2, 4, 7, 6, 1, 11, 15, 17, 14, 13, 0};
    uint32_t *obs =new uint32_t[tree.nparens / 2];

    for(unsigned int i = 0; i < (tree.nparens / 2); i++)
        obs[i] = tree.postorderselect(i);

    std::vector<uint32_t> exp_v = _uint32_array_to_vector(exp, tree.nparens / 2);
    std::vector<uint32_t> obs_v = _uint32_array_to_vector(obs, tree.nparens / 2);

    ASSERT(obs_v == exp_v);
    delete[] obs;
    SUITE_END();
}

void test_bptree_preorder() {
    SUITE_START("preorderselect");

    // fig1 from https://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    su::BPTree tree("((3,4,(6)5)2,7,((10,100)9)8)1;");
    uint32_t exp[] = {0, 1, 2, 4, 6, 7, 11, 13, 14, 15, 17};
    uint32_t *obs = new uint32_t[tree.nparens / 2];

    for(unsigned int i = 0; i < (tree.nparens / 2); i++)
        obs[i] = tree.preorderselect(i);

    std::vector<uint32_t> exp_v = _uint32_array_to_vector(exp, tree.nparens / 2);
    std::vector<uint32_t> obs_v = _uint32_array_to_vector(obs, tree.nparens / 2);

    ASSERT(obs_v == exp_v);
    delete[] obs;
    SUITE_END();
}

void test_bptree_parent() {
    SUITE_START("parent");

    // fig1 from https://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    su::BPTree tree("((3,4,(6)5)2,7,((10,100)9)8)1;");
    uint32_t exp[] = {0, 1, 1, 1, 1, 1, 6, 6, 1, 0, 0, 0, 0, 13, 14, 14, 14, 14, 13, 0};

    // all the -2 and +1 garbage is to avoid testing the root.
    uint32_t *obs = new uint32_t[tree.nparens - 2];

    for(int i = 0; i < (int(tree.nparens) - 2); i++)
        obs[i] = tree.parent(i+1);

    std::vector<uint32_t> exp_v = _uint32_array_to_vector(exp, tree.nparens - 2);
    std::vector<uint32_t> obs_v = _uint32_array_to_vector(obs, tree.nparens - 2);

    ASSERT(obs_v == exp_v);
    delete[] obs;
    SUITE_END();
}

// need a child class to test protected constructor
class test_biom_inmem : public su::biom_inmem {
        public:
            test_biom_inmem(const su::biom_inmem& other, bool _clean_on_destruction)
               : su::biom_inmem(other, _clean_on_destruction) {}
};

void test_biom_constructor() {
    SUITE_START("biom constructor");

    su::biom table("test.biom");
    uint32_t exp_n_samples = 6;
    uint32_t exp_n_obs = 5;

    std::string sids[] = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};
    std::vector<std::string> exp_sids = _string_array_to_vector(sids, exp_n_samples);

    std::string oids[] = {"GG_OTU_1", "GG_OTU_2","GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
    std::vector<std::string> exp_oids = _string_array_to_vector(oids, exp_n_obs);

    uint32_t s_indptr[] = {0, 2, 5, 9, 11, 12, 15};
    std::vector<uint32_t> exp_s_indptr = _uint32_array_to_vector(s_indptr, exp_n_samples + 1);

    uint32_t o_indptr[] = {0, 1, 6, 9, 13, 15};
    std::vector<uint32_t> exp_o_indptr = _uint32_array_to_vector(o_indptr, exp_n_obs + 1);

    uint32_t exp_nnz = 15;

    ASSERT(table.n_samples == exp_n_samples);
    ASSERT(table.n_obs == exp_n_obs);
    ASSERT(table.nnz == exp_nnz);
    ASSERT(table.get_sample_ids() == exp_sids);
    ASSERT(table.get_obs_ids() == exp_oids);
    ASSERT(table.is_sample_indptr(exp_s_indptr));
    ASSERT(table.is_obs_indptr(exp_o_indptr));

    test_biom_inmem table_copy(table,true);
    ASSERT(table_copy.n_samples == exp_n_samples);
    ASSERT(table_copy.n_obs == exp_n_obs);
    ASSERT(table_copy.get_sample_ids() == exp_sids);
    ASSERT(table_copy.get_obs_ids() == exp_oids);

    uint32_t t_n_samples = su::biom::load_n_samples("test.biom");
    ASSERT(t_n_samples == exp_n_samples);

    SUITE_END();
}

void _exercise_get_obs_data(su::biom_interface &table) {
    double exp0[] = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    std::vector<double> exp0_vec = _double_array_to_vector(exp0, 6);
    double exp1[] = {5.0, 1.0, 0.0, 2.0, 3.0, 1.0};
    std::vector<double> exp1_vec = _double_array_to_vector(exp1, 6);
    double exp2[] = {0.0, 0.0, 1.0, 4.0, 0.0, 2.0};
    std::vector<double> exp2_vec = _double_array_to_vector(exp2, 6);
    double exp3[] = {2.0, 1.0, 1.0, 0.0, 0.0, 1.0};
    std::vector<double> exp3_vec = _double_array_to_vector(exp3, 6);
    double exp3a[] = {1.0, 1.0, 0.0};
    std::vector<double> exp3a_vec = _double_array_to_vector(exp3a, 3);
    double exp3b[] = {1.0/3, 1.0/4, 0.0};
    std::vector<double> exp3b_vec = _double_array_to_vector(exp3b, 3);
    double exp4[] = {0.0, 1.0, 1.0, 0.0, 0.0, 0.0};
    std::vector<double> exp4_vec = _double_array_to_vector(exp4, 6);

    double *out = (double*)malloc(sizeof(double) * 6);
    double *out2 = (double*)malloc(sizeof(double) * 3);
    std::vector<double> obs_vec;

    table.get_obs_data(std::string("GG_OTU_1").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp0_vec));

    table.get_obs_data(std::string("GG_OTU_2").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp1_vec));

    table.get_obs_data(std::string("GG_OTU_3").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp2_vec));

    table.get_obs_data(std::string("GG_OTU_4").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp3_vec));

    table.get_obs_data_range(std::string("GG_OTU_4").c_str(), 1,1+3,false, out2);
    obs_vec = _double_array_to_vector(out2, 3);
    ASSERT(vec_almost_equal(obs_vec, exp3a_vec));

    table.get_obs_data_range(std::string("GG_OTU_4").c_str(), 1,1+3,true, out2);
    obs_vec = _double_array_to_vector(out2, 3);
    ASSERT(vec_almost_equal(obs_vec, exp3b_vec));

    table.get_obs_data(std::string("GG_OTU_5").c_str(), out);
    obs_vec = _double_array_to_vector(out, 6);
    ASSERT(vec_almost_equal(obs_vec, exp4_vec));

    free(out2);
    free(out);
}

void test_biom_filter() {
    SUITE_START("biom filter");

    su::biom table_full("test.biom");
    su::biom_inmem table(table_full,5.0);

    uint32_t exp_n_samples = 2;
    uint32_t exp_n_obs = 3;

    std::string sids[] = {"Sample1", "Sample4"};
    std::vector<std::string> exp_sids = _string_array_to_vector(sids, exp_n_samples);

    std::string oids[] = {"GG_OTU_2", "GG_OTU_3", "GG_OTU_4"};
    std::vector<std::string> exp_oids = _string_array_to_vector(oids, exp_n_obs);

    ASSERT(table.n_samples == exp_n_samples);
    ASSERT(table.n_obs == exp_n_obs);
    ASSERT(table.get_sample_ids() == exp_sids);
    ASSERT(table.get_obs_ids() == exp_oids);

    const double *counts = table.get_sample_counts();
    for (uint32_t i=0; i<table.n_samples; i++) ASSERT(counts[i]>=5.0);

    double exp1[] = {5.0, 2.0};
    std::vector<double> exp1_vec = _double_array_to_vector(exp1, 2);
    double exp2[] = {0.0, 4.0};
    std::vector<double> exp2_vec = _double_array_to_vector(exp2, 2);
    double exp3[] = {2.0, 0.0};
    std::vector<double> exp3_vec = _double_array_to_vector(exp3, 2);

    double *out = (double*)malloc(sizeof(double) * 2);
    std::vector<double> obs_vec;

    table.get_obs_data(std::string("GG_OTU_2").c_str(), out);
    obs_vec = _double_array_to_vector(out, 2);
    ASSERT(vec_almost_equal(obs_vec, exp1_vec));

    table.get_obs_data(std::string("GG_OTU_3").c_str(), out);
    obs_vec = _double_array_to_vector(out, 2);
    ASSERT(vec_almost_equal(obs_vec, exp2_vec));

    table.get_obs_data(std::string("GG_OTU_4").c_str(), out);
    obs_vec = _double_array_to_vector(out, 2);
    ASSERT(vec_almost_equal(obs_vec, exp3_vec));

    SUITE_END();
}

void test_biom_constructor_from_sparse() {
    SUITE_START("biom from sparse constructor");
    uint32_t index[] = {2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2};
    uint32_t indptr[] = {0,  1,  6,  9, 13, 15};
    double data[] = {1., 5., 1., 2., 3., 1., 1., 4., 2., 2., 1., 1., 1., 1., 1.};
    const char* obs_ids[] = {"GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
    const char* samp_ids[] = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};

    su::biom_inmem table(obs_ids, samp_ids, index, indptr, data, 5, 6);
    _exercise_get_obs_data(table);
   
    SUITE_END();
}

void test_biom_constructor_from_dense() {
    SUITE_START("biom from dense constructor");
    double sample1[] = {0, 5, 0, 2, 0};
    double sample2[] = {0, 1, 0, 1, 1};
    double sample3[] = {1, 0, 1, 1, 1};
    double sample4[] = {0, 2, 4, 0, 0};
    double sample5[] = {0, 3, 0, 0, 0};
    double sample6[] = {0, 1, 2, 1, 0};
    const char* obs_ids[] = {"GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
    const char* samp_ids[] = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};

    double* data_ptrs[] = {sample1,sample2,sample3,sample4,sample5, sample6};
    su::biom_inmem table(obs_ids, samp_ids, data_ptrs, 5, 6);
    _exercise_get_obs_data(table);
   
    SUITE_END();
}

void test_biom_nullary() {
    SUITE_START("biom nullary");
    su::biom table;
    SUITE_END();
}

void test_bptree_nullary() {
    SUITE_START("bptree nullary");
    su::BPTree tree;
    SUITE_END();
}

void test_biom_get_obs_data() {
    SUITE_START("biom get obs data");

    su::biom table("test.biom");
    _exercise_get_obs_data(table);
    SUITE_END();
}


void test_bptree_leftchild() {
    SUITE_START("test bptree left child");
    su::BPTree tree("((3,4,(6)5)2,7,((10,100)9)8)1;");

    uint32_t exp[] = {1, 2, 0, 0, 7, 0, 0, 14, 15, 0, 0};
    std::vector<bool> structure = tree.get_structure();

    uint32_t exp_pos = 0;
    for(unsigned int i = 0; i < tree.nparens; i++) {
        if(structure[i])
            ASSERT(tree.leftchild(i) == exp[exp_pos++]);
    }
    SUITE_END();
}

void test_bptree_rightchild() {
    SUITE_START("test bptree right child");
    su::BPTree tree("((3,4,(6)5)2,7,((10,100)9)8)1;");

    uint32_t exp[] = {13, 6, 0, 0, 7, 0, 0, 14, 17, 0, 0};
    std::vector<bool> structure = tree.get_structure();

    uint32_t exp_pos = 0;
    for(unsigned int i = 0; i < tree.nparens; i++) {
        if(structure[i])
            ASSERT(tree.rightchild(i) == exp[exp_pos++]);
    }
    SUITE_END();
}

void test_bptree_rightsibling() {
    SUITE_START("test bptree rightsibling");
    su::BPTree tree("((3,4,(6)5)2,7,((10,100)9)8)1;");

    uint32_t exp[] = {0, 11, 4, 6, 0, 0, 13, 0, 0, 17, 0};
    std::vector<bool> structure = tree.get_structure();

    uint32_t exp_pos = 0;
    for(unsigned int i = 0; i < tree.nparens; i++) {
        if(structure[i])
            ASSERT(tree.rightsibling(i) == exp[exp_pos++]);
    }
    SUITE_END();
}

void test_propstack_constructor() {
    SUITE_START("test propstack constructor");
    su::PropStack<double> ps(10);
    // nothing to test directly...
    SUITE_END();
}

void test_propstack_push_and_pop() {
    SUITE_START("test propstack push and pop");
    su::PropStack<double> ps(10);

    double *vec1 = ps.pop(1);
    double *vec2 = ps.pop(2);
    double *vec3 = ps.pop(3);
    double *vec1_obs;
    double *vec2_obs;
    double *vec3_obs;

    ps.push(1);
    ps.push(2);
    ps.push(3);

    vec3_obs = ps.pop(4);
    vec2_obs = ps.pop(5);
    vec1_obs = ps.pop(6);

    ASSERT(vec1 == vec1_obs);
    ASSERT(vec2 == vec2_obs);
    ASSERT(vec3 == vec3_obs);
    SUITE_END();
}

void test_propstack_get() {
    SUITE_START("test propstack get");
    su::PropStack<double> ps(10);

    double *vec1 = ps.pop(1);
    double *vec2 = ps.pop(2);
    double *vec3 = ps.pop(3);

    double *vec1_obs = ps.get(1);
    double *vec2_obs = ps.get(2);
    double *vec3_obs = ps.get(3);

    ASSERT(vec1 == vec1_obs);
    ASSERT(vec2 == vec2_obs);
    ASSERT(vec3 == vec3_obs);
    SUITE_END();
}

void test_unifrac_set_proportions() {
    SUITE_START("test unifrac set proportions");
    //                           0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
    //                           ( ( ) ( ( ) ( ) ) ( ( ) ( ) ) )
    su::BPTree tree("(GG_OTU_1,(GG_OTU_2,GG_OTU_3),(GG_OTU_5,GG_OTU_4));");
    su::biom table("test.biom");
    su::PropStack<double> ps(table.n_samples);

    double *obs = ps.pop(4); // GG_OTU_2
    double exp4[] = {0.714285714286, 0.333333333333, 0.0, 0.333333333333, 1.0, 0.25};
    set_proportions(obs, tree, 4, table, ps);
    for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp4[i]) < 0.000001);

    obs = ps.pop(6); // GG_OTU_3
    double exp6[] = {0.0, 0.0, 0.25, 0.666666666667, 0.0, 0.5};
    set_proportions(obs, tree, 6, table, ps);
    for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp6[i]) < 0.000001);

    obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
    double exp3[] = {0.71428571, 0.33333333, 0.25, 1.0, 1.0, 0.75};
    set_proportions(obs, tree, 3, table, ps);
    for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp3[i]) < 0.000001);
    SUITE_END();
}

void test_unifrac_set_proportions_range() {
    SUITE_START("test unifrac set proportions range");
    //                           0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
    //                           ( ( ) ( ( ) ( ) ) ( ( ) ( ) ) )
    su::BPTree tree("(GG_OTU_1,(GG_OTU_2,GG_OTU_3),(GG_OTU_5,GG_OTU_4));");
    su::biom table("test.biom");

    const double exp4[] = {0.714285714286, 0.333333333333, 0.0, 0.333333333333, 1.0, 0.25};
    const double exp6[] = {0.0, 0.0, 0.25, 0.666666666667, 0.0, 0.5};
    const double exp3[] = {0.71428571, 0.33333333, 0.25, 1.0, 1.0, 0.75};


    // first the whole table
    {
      su::PropStack<double> ps(table.n_samples);

      double *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, 0, table.n_samples, ps);
      for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, 0, table.n_samples, ps);
      for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, 0, table.n_samples, ps);
      for(unsigned int i = 0; i < table.n_samples; i++)
        ASSERT(fabs(obs[i] - exp3[i]) < 0.000001);
    }

    // beginning
    {
      su::PropStack<double> ps(3);

      double *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, 0, 3, ps);
      for(unsigned int i = 0; i < 3; i++)
        ASSERT(fabs(obs[i] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, 0, 3, ps);
      for(unsigned int i = 0; i < 3; i++)
        ASSERT(fabs(obs[i] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, 0, 3, ps);
      for(unsigned int i = 0; i < 3; i++)
        ASSERT(fabs(obs[i] - exp3[i]) < 0.000001);
    }

    
    // end
    {
      su::PropStack<double> ps(4);

      double *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, 2, table.n_samples, ps);
      for(unsigned int i = 2; i < table.n_samples; i++)
        ASSERT(fabs(obs[i-2] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, 2, table.n_samples, ps);
      for(unsigned int i = 2; i < table.n_samples; i++)
        ASSERT(fabs(obs[i-2] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, 2, table.n_samples, ps);
      for(unsigned int i = 2; i < table.n_samples; i++)
        ASSERT(fabs(obs[i-2] - exp3[i]) < 0.000001);
    }

    
    // middle
    {
      const unsigned int start = 1;
      const unsigned int end = 4;
      su::PropStack<double> ps(end-start);

      double *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, start, end, ps);
      for(unsigned int i =start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, start, end, ps);
      for(unsigned int i = start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, start, end, ps);
      for(unsigned int i = start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp3[i]) < 0.000001);
    }

    SUITE_END();
}

void test_unifrac_set_proportions_range_float() {
    SUITE_START("test unifrac set proportions range float");
    //                           0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
    //                           ( ( ) ( ( ) ( ) ) ( ( ) ( ) ) )
    su::BPTree tree("(GG_OTU_1,(GG_OTU_2,GG_OTU_3),(GG_OTU_5,GG_OTU_4));");
    su::biom table("test.biom");

    const float exp4[] = {0.714285714286, 0.333333333333, 0.0, 0.333333333333, 1.0, 0.25};
    const float exp6[] = {0.0, 0.0, 0.25, 0.666666666667, 0.0, 0.5};
    const float exp3[] = {0.71428571, 0.33333333, 0.25, 1.0, 1.0, 0.75};

    // just midle
    {
      const unsigned int start = 1;
      const unsigned int end = 4;
      su::PropStack<float> ps(end-start);

      float *obs = ps.pop(4); // GG_OTU_2
      set_proportions_range(obs, tree, 4, table, start, end, ps);
      for(unsigned int i =start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp4[i]) < 0.000001);

      obs = ps.pop(6); // GG_OTU_3
      set_proportions_range(obs, tree, 6, table, start, end, ps);
      for(unsigned int i = start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp6[i]) < 0.000001);

      obs = ps.pop(3); // node containing GG_OTU_2 and GG_OTU_3
      set_proportions_range(obs, tree, 3, table, start, end, ps);
      for(unsigned int i = start; i < end; i++)
        ASSERT(fabs(obs[i-start] - exp3[i]) < 0.000001);
    }

    SUITE_END();
}



void test_unifrac_deconvolute_stripes() {
    SUITE_START("test deconvolute stripes");
    std::vector<double*> stripes;
    double s1[] = {1, 1, 1, 1, 1, 1};
    double s2[] = {2, 2, 2, 2, 2, 2};
    double s3[] = {3, 3, 3, 3, 3, 3};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);

    double exp[6][6] = { {0, 1, 2, 3, 2, 1},
                         {1, 0, 1, 2, 3, 2},
                         {2, 1, 0, 1, 2, 3},
                         {3, 2, 1, 0, 1, 2},
                         {2, 3, 2, 1, 0, 1},
                         {1, 2, 3, 2, 1, 0} };
    double **obs = su::deconvolute_stripes(stripes, 6);
    for(unsigned int i = 0; i < 6; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(exp[i][j] == obs[i][j]);
        }
    }
    free(obs);
    SUITE_END();
}
#endif

void test_unifrac_stripes_to_condensed_form_even() {
    SUITE_START("test stripes_to_condensed_form even samples");
    std::vector<double*> stripes;
    double s1[] = {0,  9, 17, 24, 30, 35, 39, 42, 44,  8};
    double s2[] = {1, 10, 18, 25, 31, 36, 40, 43,  7, 16};
    double s3[] = {2, 11, 19, 26, 32, 37, 41,  6, 15, 23};
    double s4[] = {3, 12, 20, 27, 33, 38,  5, 14, 22, 29};
    double s5[] = {4, 13, 21, 28, 34,  4, 13, 21, 28, 34};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[45] = {/* 0, */  0,  1,  2,  3,  4,  5,  6,  7,  8,
                      /* *,  0, */  9, 10, 11, 12, 13, 14, 15, 16,
                      /* *,  *,  0, */ 17, 18, 19, 20, 21, 22, 23,
                      /* *,  *,  *,  0, */ 24, 25, 26, 27, 28, 29,
                      /* *,  *,  *,  *,  0, */ 30, 31, 32, 33, 34,
                      /* *,  *,  *,  *,  *,  0, */ 35, 36, 37, 38,
                      /* *,  *,  *,  *,  *,  *,  0, */ 39, 40, 41,
                      /* *,  *,  *,  *,  *,  *,  *,  0, */ 42, 43,
                      /* *,  *,  *,  *,  *,  *,  *,  *,  0, */ 44};
                      /* *,  *,  *,  *,  *,  *,  *,  *,  *, *,  0 */

    double *obs = NULL;
#ifndef API_ONLY
    obs = (double*)malloc(sizeof(double) * 45);
    su::stripes_to_condensed_form(stripes, 10, obs, 0, 5);
    for(unsigned int i = 0; i < 45; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
#endif
    // test internal version, too
    double* cstripes[] = {s1, s2, s3, s4 ,s5};
    obs = (double*)malloc(sizeof(double) * 45);
    _testv_stripes_to_condensed_form(cstripes, 10, 5, obs);
    for(unsigned int i = 0; i < 45; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
    SUITE_END();
}

void test_unifrac_stripes_to_condensed_form_odd() {
    SUITE_START("test stripes_to_condensed_form odd samples");
    std::vector<double*> stripes;
    double s1[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 0};
    double s2[] = {20, 19, 18, 17, 16, 15, 14 ,13, 12, 11, 1};
    double s3[] = {21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 2};
    double s4[] = {40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 3};
    double s5[] = {41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 4};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[55] = {/* 0, */ 1, 20, 21, 40, 41, 47, 33, 29, 11,  0,
                      /* 1,  0, */ 2, 19, 22, 39, 42, 48, 32, 30,  1,
                      /*20,  2,  0, */ 3, 18, 23, 38, 43, 49, 31,  2,
                      /*21, 19,  3,  0, */ 4, 17, 24, 37, 44, 50,  3,
                      /*40, 22, 18,  4,  0, */ 5, 16, 25 ,36, 45 , 4,
                      /*41, 39, 23, 17,  5,  0, */ 6, 15, 26, 35, 46,
                      /*47, 42, 38, 24, 16,  6,  0, */ 7, 14, 27, 34,
                      /*33, 48, 43, 37, 25, 15,  7,  0,*/  8, 13, 28,
                      /*29, 32, 49, 44, 36, 26, 14,  8,  0, */ 9, 12,
                      /*11, 30, 31, 50, 45, 35, 27, 13,  9,  0,*/ 10};
                      /* 0,  1,  2,  3,  4, 46, 34, 28, 12, 10,  0}; */
    double *obs = NULL;
#ifndef API_ONLY
    obs = (double*)malloc(sizeof(double) * 55);
    su::stripes_to_condensed_form(stripes, 11, obs, 0, 5);
    for(unsigned int i = 0; i < 55; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
#endif
    // test internal version, too
    double* cstripes[] = {s1, s2, s3, s4 ,s5};
    obs = (double*)malloc(sizeof(double) * 55);
    _testv_stripes_to_condensed_form(cstripes, 11, 5, obs);
    for(unsigned int i = 0; i < 55; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
    SUITE_END();
}

void test_unifrac_stripes_to_condensed_form_odd2() {
    SUITE_START("test stripes_to_condensed_form odd(2) samples");
    std::vector<double*> stripes;
    double s1[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9};
    double s2[] = {18, 17, 16, 15, 14, 13, 12 ,11, 10};
    double s3[] = {19, 20, 21, 22, 23, 24, 25, 26, 27};
    double s4[] = {36, 35, 34, 33, 32, 31, 30, 29, 28};
    double s5[] = {31, 30, 29, 28, 36, 35, 34, 33, 32};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[36] = {/* 0, */ 1, 18, 19, 36, 31, 25, 11,  9,
                      /* 1,  0, */ 2, 17, 20, 35, 30, 26, 10,
                      /*20,  2,  0, */ 3, 16, 21, 34, 29, 27,
                      /*21, 19,  3,  0, */ 4, 15, 22, 33, 28,
                      /*40, 22, 18,  4,  0, */ 5, 14, 23 ,32,
                      /*41, 39, 23, 17,  5,  0, */ 6, 13, 24,
                      /*47, 42, 38, 24, 16,  6,  0, */ 7, 12,
                      /*47, 42, 38, 24, 16,  6,  7, 0, */  8};
                      /* 0,  1,  2,  3,  4, 46, 34, 8,  8, 0}; */
    double *obs = NULL;
#ifndef API_ONLY
    obs = (double*)malloc(sizeof(double) * 36);
    su::stripes_to_condensed_form(stripes, 9, obs, 0, 5);
    for(unsigned int i = 0; i < 36; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
#endif
    // test internal version, too
    double* cstripes[] = {s1, s2, s3, s4 ,s5};
    obs = (double*)malloc(sizeof(double) * 36);
    _testv_stripes_to_condensed_form(cstripes, 9, 5, obs);
    for(unsigned int i = 0; i < 36; i++) {
        ASSERT(exp[i] == obs[i]);
    }
    free(obs);
    SUITE_END();
}

#ifndef API_ONLY
class ValidatedMemoryStripes : public su::MemoryStripes {
        private:
           const uint32_t n_stripes;
           mutable std::vector<uint8_t> stripe_status; // 0 new, 1 allocated, 2 deallocated, 3 reallocated, 6 deallocate after rellocation
        public:
           ValidatedMemoryStripes(uint32_t _n_stripes, std::vector<double*> &_stripes) 
           : su::MemoryStripes(_stripes) 
           , n_stripes(_n_stripes)
           , stripe_status(n_stripes)
           {
             for (uint32_t i=0; i<n_stripes; i++) stripe_status[i] = 0;
           }

           virtual const double *get_stripe(uint32_t stripe) const {
              stripe_status[stripe]|=1; 
              return su::MemoryStripes::get_stripe(stripe);
           }
           virtual void release_stripe(uint32_t stripe) const { 
              if (stripe_status[stripe]<2) {
                 stripe_status[stripe]=2;
              } else {
                 stripe_status[stripe]=6;
              }
           }

           bool allInitialized() const {
             bool out = true;
             for (uint32_t i=0; i<n_stripes; i++) out &= (stripe_status[i] != 0);
             return out;
           }

           bool allDealocated() const {
             bool out = true;
             for (uint32_t i=0; i<n_stripes; i++) out &= ((stripe_status[i]&2) == 2);
             return out;
           }

           bool anyRealocated() const {
             bool out = false;
             for (uint32_t i=0; i<n_stripes; i++) out |= (stripe_status[i] >2);
             return out;
           }


};


void test_unifrac_stripes_to_matrix_even() {
    SUITE_START("test stripes_to_matrix even samples");
    std::vector<double*> stripes;
    double s1[] = {0,  9, 17, 24, 30, 35, 39, 42, 44,  8};
    double s2[] = {1, 10, 18, 25, 31, 36, 40, 43,  7, 16};
    double s3[] = {2, 11, 19, 26, 32, 37, 41,  6, 15, 23};
    double s4[] = {3, 12, 20, 27, 33, 38,  5, 14, 22, 29};
    double s5[] = {4, 13, 21, 28, 34,  4, 13, 21, 28, 34};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    // test also double to float conversion
    float exp[100] = {0,  0,  1,  2,  3,  4,  5,  6,  7,  8, 
                      0,  0,  9, 10, 11, 12, 13, 14, 15, 16,  
                      1,  9,  0, 17, 18, 19, 20, 21, 22, 23,
                      2, 10, 17,  0, 24, 25, 26, 27, 28, 29, 
                      3, 11, 18, 24,  0, 30, 31, 32, 33, 34,
                      4, 12, 19, 25, 30,  0, 35, 36, 37, 38,
                      5, 13, 20, 26, 31, 35,  0, 39, 40, 41,
                      6, 14, 21, 27, 32, 36, 39,  0, 42, 43,
                      7, 15, 22, 28, 33, 37, 40, 42,  0, 44,
                      8, 16, 23, 29, 34, 38, 41, 43, 44,  0};
    {
      float *obs = (float*)malloc(sizeof(float) * 100);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix_fp32(vs, 10, 5, obs);
      for(unsigned int i = 0; i < 100; i++) {
        ASSERT(exp[i] == obs[i]);
      }

      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true); 
      ASSERT(vs.anyRealocated() == false);

      free(obs);
    }

    { // small tiles
      float *obs = (float*)malloc(sizeof(float) * 100);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix_fp32(vs, 10, 5, obs, 4);
      for(unsigned int i = 0; i < 100; i++) {
        ASSERT(exp[i] == obs[i]);
      }

      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
    
      free(obs);
    }

    { // large tiles
      float *obs = (float*)malloc(sizeof(float) * 100);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix_fp32(vs, 10, 5, obs, 128);
      for(unsigned int i = 0; i < 100; i++) {
        ASSERT(exp[i] == obs[i]);
      }

      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
    
      free(obs);
    }


    // test also intermediate, 2-step procedure
    double *obsC = (double*)malloc(sizeof(double) * 45);
    su::stripes_to_condensed_form(stripes, 10, obsC, 0, 5);

    float *obs2 = (float*)malloc(sizeof(float) * 100);
    su::condensed_form_to_matrix_fp32(obsC, 10, obs2);

    for(unsigned int i = 0; i < 100; i++) {
        ASSERT(exp[i] == obs2[i]);
    }

    free(obs2);
    free(obsC);
    SUITE_END();
}

void test_unifrac_stripes_to_matrix_odd() {
    SUITE_START("test stripes_to_matrix odd samples");
    std::vector<double*> stripes;
    double s1[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 0};
    double s2[] = {20, 19, 18, 17, 16, 15, 14 ,13, 12, 11, 1};
    double s3[] = {21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 2};
    double s4[] = {40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 3};
    double s5[] = {41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 4};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[121] = { 0,  1, 20, 21, 40, 41, 47, 33, 29, 11,  0, 
                        1,  0,  2, 19, 22, 39, 42, 48, 32, 30,  1,
                       20,  2,  0,  3, 18, 23, 38, 43, 49, 31,  2,
                       21, 19,  3,  0,  4, 17, 24, 37, 44, 50,  3,
                       40, 22, 18,  4,  0,  5, 16, 25 ,36, 45 , 4, 
                       41, 39, 23, 17,  5,  0,  6, 15, 26, 35, 46,
                       47, 42, 38, 24, 16,  6,  0,  7, 14, 27, 34,
                       33, 48, 43, 37, 25, 15,  7,  0,  8, 13, 28, 
                       29, 32, 49, 44, 36, 26, 14,  8,  0,  9, 12,
                       11, 30, 31, 50, 45, 35, 27, 13,  9,  0, 10,
                        0,  1,  2,  3,  4, 46, 34, 28, 12, 10,  0};

    {
      double *obs = (double*)malloc(sizeof(double) * 121);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 11, 5, obs);
      for(unsigned int i = 0; i < 121; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    { // small tiling
      double *obs = (double*)malloc(sizeof(double) * 121);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 11, 5, obs,4);
      for(unsigned int i = 0; i < 121; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    { // large tiling
      double *obs = (double*)malloc(sizeof(double) * 121);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 11, 5, obs,128);
      for(unsigned int i = 0; i < 121; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }


    // test also intermediate, 2-step procedure
    double *obsC = (double*)malloc(sizeof(double) * 55);
    su::stripes_to_condensed_form(stripes, 11, obsC, 0, 5);

    double *obs2 = (double*)malloc(sizeof(double) * 121);
    su::condensed_form_to_matrix(obsC, 11, obs2);
    
    for(unsigned int i = 0; i < 121; i++) {
        ASSERT(exp[i] == obs2[i]);
    }

    free(obs2);
    free(obsC);
    SUITE_END();
}

void test_unifrac_stripes_to_matrix_odd2() {
    SUITE_START("test stripes_to_matrix odd(2) samples");
    std::vector<double*> stripes;
    double s1[] = { 1,  2,  3,  4,  5,  6,  7,  8,  9};
    double s2[] = {18, 17, 16, 15, 14, 13, 12 ,11, 10};
    double s3[] = {19, 20, 21, 22, 23, 24, 25, 26, 27};
    double s4[] = {36, 35, 34, 33, 32, 31, 30, 29, 28};
    double s5[] = {31, 30, 29, 28, 36, 35, 34, 33, 32};
    stripes.push_back(s1);
    stripes.push_back(s2);
    stripes.push_back(s3);
    stripes.push_back(s4);
    stripes.push_back(s5);

    double exp[81] = { 0,  1, 18, 19, 36, 31, 25, 11,  9,
                       1,  0,  2, 17, 20, 35, 30, 26, 10,
                      18,  2,  0,  3, 16, 21, 34, 29, 27,
                      19, 17,  3,  0,  4, 15, 22, 33, 28,
                      36, 20, 16,  4,  0,  5, 14, 23 ,32,
                      31, 35, 21, 15,  5,  0,  6, 13, 24,
                      25, 30, 34, 22, 14,  6,  0,  7, 12,
                      11, 26, 29, 33, 23, 13,  7,  0,  8,
                       9, 10, 27, 28, 32, 24, 12,  8,  0};

    {
      double *obs = (double*)malloc(sizeof(double) * 81);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 9, 5, obs);
      for(unsigned int i = 0; i < 81; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    { // small tile
      double *obs = (double*)malloc(sizeof(double) * 81);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 9, 5, obs,4);
      for(unsigned int i = 0; i < 81; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    { // large tile
      double *obs = (double*)malloc(sizeof(double) * 81);
      ValidatedMemoryStripes vs(5,stripes);
      su::stripes_to_matrix(vs, 9, 5, obs,128);
      for(unsigned int i = 0; i < 81; i++) {
        ASSERT(exp[i] == obs[i]);
      }
      ASSERT(vs.allInitialized() == true);
      ASSERT(vs.allDealocated() == true);
      ASSERT(vs.anyRealocated() == false);
      free(obs);
    }

    // test also intermediate, 2-step procedure
    double *obsC = (double*)malloc(sizeof(double) * 36);
    su::stripes_to_condensed_form(stripes, 9, obsC, 0, 5);

    double *obs2 = (double*)malloc(sizeof(double) * 81);
    su::condensed_form_to_matrix(obsC, 9, obs2);
    
    for(unsigned int i = 0; i < 81; i++) {
        ASSERT(exp[i] == obs2[i]);
    }

    free(obs2);
    free(obsC);
    SUITE_END();
}
#endif


void test_unnormalized_weighted_unifrac() {
    SUITE_START("test unnormalized weighted unifrac");

#ifndef API_ONLY
    su::BPTree tree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table("test.biom");
#endif

    std::vector<double*> exp;
    double stride1[] = {1.52380952, 1.25, 2.75, 1.33333333, 2., 1.07142857};
    double stride2[] = {2.17857143, 2.66666667, 3.25, 1.0, 1.14285714, 1.83333333};
    double stride3[] = {1.9047619, 2.66666667, 1.75, 1.9047619, 2.66666667, 1.75};
#ifndef API_ONLY
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = false; task_p.normalize_sample_counts = true;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::weighted_unnormalized,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
    }
#endif

    // repeat using the API
    double* cstripes[] = {stride1, stride2, stride3};
    mat_t *res = NULL;
    ComputeStatus rc = one_off("test.biom", "test.tre",
		               "weighted_unnormalized", false, 0.0,
			       false, 1, &res);
    ASSERT(rc == okay);
    double *expS = (double*)malloc(sizeof(double) * 15);
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    free(expS);
    destroy_mat(&res);

    SUITE_END();
}

void test_generalized_unifrac() {
    SUITE_START("test generalized unifrac");

#ifndef API_ONLY
    su::BPTree tree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table("test.biom");
#endif

    // weighted normalized unifrac as computed above
    std::vector<double*> w_exp;
    double w_stride1[] = {0.38095238, 0.33333333, 0.73333333, 0.33333333, 0.5, 0.26785714};
    double w_stride2[] = {0.58095238, 0.66666667, 0.86666667, 0.25, 0.28571429, 0.45833333};
    double w_stride3[] = {0.47619048, 0.66666667, 0.46666667, 0.47619048, 0.66666667, 0.46666667};
#ifndef API_ONLY
    w_exp.push_back(w_stride1);
    w_exp.push_back(w_stride2);
    w_exp.push_back(w_stride3);
    std::vector<double*> w_strides = su::make_strides(6);
    std::vector<double*> w_strides_total = su::make_strides(6);
    su::task_parameters w_task_p;
    w_task_p.start = 0; w_task_p.stop = 3; w_task_p.tid = 0; w_task_p.n_samples = 6; w_task_p.bypass_tips = false; w_task_p.normalize_sample_counts = true;
    w_task_p.g_unifrac_alpha = 1.0;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(w_task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::generalized,
                        false,
                        std::ref(w_strides),
                        std::ref(w_strides_total),
                        std::ref(tasks));
#endif

    // as computed by GUniFrac v1.0
    //          Sample1   Sample2   Sample3   Sample4   Sample5   Sample6
    //Sample1 0.0000000 0.4408392 0.6886965 0.7060606 0.5833333 0.3278410
    //Sample2 0.4408392 0.0000000 0.5102041 0.7500000 0.8000000 0.5208125
    //Sample3 0.6886965 0.5102041 0.0000000 0.8649351 0.9428571 0.5952381
    //Sample4 0.7060606 0.7500000 0.8649351 0.0000000 0.5000000 0.4857143
    //Sample5 0.5833333 0.8000000 0.9428571 0.5000000 0.0000000 0.7485714
    //Sample6 0.3278410 0.5208125 0.5952381 0.4857143 0.7485714 0.0000000
    std::vector<double*> d0_exp;
    double d0_stride1[] = {0.4408392, 0.5102041, 0.8649351, 0.5000000, 0.7485714, 0.3278410};
    double d0_stride2[] = {0.6886965, 0.7500000, 0.9428571, 0.4857143, 0.5833333, 0.5208125};
    double d0_stride3[] = {0.7060606, 0.8000000, 0.5952381, 0.7060606, 0.8000000, 0.5952381};
#ifndef API_ONLY
    d0_exp.push_back(d0_stride1);
    d0_exp.push_back(d0_stride2);
    d0_exp.push_back(d0_stride3);
    std::vector<double*> d0_strides = su::make_strides(6);
    std::vector<double*> d0_strides_total = su::make_strides(6);
    su::task_parameters d0_task_p;
    d0_task_p.start = 0; d0_task_p.stop = 3; d0_task_p.tid = 0; d0_task_p.n_samples = 6; d0_task_p.bypass_tips = false; d0_task_p.normalize_sample_counts = true;
    d0_task_p.g_unifrac_alpha = 0.0;

    tasks.clear();
    tasks.push_back(d0_task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::generalized,
                        false,
                        std::ref(d0_strides),
                        std::ref(d0_strides_total),
                        std::ref(tasks));
#endif

    // as computed by GUniFrac v1.0
    //          Sample1   Sample2   Sample3   Sample4   Sample5   Sample6
    //Sample1 0.0000000 0.4040518 0.6285560 0.5869439 0.4082483 0.2995673
    //Sample2 0.4040518 0.0000000 0.4160597 0.7071068 0.7302479 0.4860856
    //Sample3 0.6285560 0.4160597 0.0000000 0.8005220 0.9073159 0.5218198
    //Sample4 0.5869439 0.7071068 0.8005220 0.0000000 0.4117216 0.3485667
    //Sample5 0.4082483 0.7302479 0.9073159 0.4117216 0.0000000 0.6188282
    //Sample6 0.2995673 0.4860856 0.5218198 0.3485667 0.6188282 0.0000000
    std::vector<double*> d05_exp;
    double d05_stride1[] = {0.4040518, 0.4160597, 0.8005220, 0.4117216, 0.6188282, 0.2995673};
    double d05_stride2[] = {0.6285560, 0.7071068, 0.9073159, 0.3485667, 0.4082483, 0.4860856};
    double d05_stride3[] = {0.5869439, 0.7302479, 0.5218198, 0.5869439, 0.7302479, 0.5218198};
#ifndef API_ONLY
    d05_exp.push_back(d05_stride1);
    d05_exp.push_back(d05_stride2);
    d05_exp.push_back(d05_stride3);
    std::vector<double*> d05_strides = su::make_strides(6);
    std::vector<double*> d05_strides_total = su::make_strides(6);
    su::task_parameters d05_task_p;
    d05_task_p.start = 0; d05_task_p.stop = 3; d05_task_p.tid = 0; d05_task_p.n_samples = 6; d05_task_p.bypass_tips = false; d05_task_p.normalize_sample_counts = true;
    d05_task_p.g_unifrac_alpha = 0.5;

    tasks.clear();
    tasks.push_back(d05_task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::generalized,
                        false,
                        std::ref(d05_strides),
                        std::ref(d05_strides_total),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(w_strides[i][j] - w_exp[i][j]) < 0.000001);
            ASSERT(fabs(d0_strides[i][j] - d0_exp[i][j]) < 0.000001);
            ASSERT(fabs(d05_strides[i][j] - d05_exp[i][j]) < 0.000001);
        }
        free(w_strides[i]);
        free(d0_strides[i]);
        free(d05_strides[i]);
    }
#endif

    // repeat using the API
    double* cstripes[3];
    double *expS = (double*)malloc(sizeof(double) * 15);

    mat_t *res = NULL;
    ComputeStatus rc;

    rc =  one_off("test.biom", "test.tre",
	               "generalized", false, 1.0,
		       false, 1, &res);
    ASSERT(rc == okay);
    cstripes[0] = w_stride1; cstripes[1] = w_stride2; cstripes[2] = w_stride3;
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    destroy_mat(&res);

    rc =  one_off("test.biom", "test.tre",
	               "generalized", false, 0.0,
		       false, 1, &res);
    ASSERT(rc == okay);
    cstripes[0] = d0_stride1; cstripes[1] = d0_stride2; cstripes[2] = d0_stride3;
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    destroy_mat(&res);

    rc =  one_off("test.biom", "test.tre",
	               "generalized", false, 0.5,
		       false, 1, &res);
    ASSERT(rc == okay);
    cstripes[0] = d05_stride1; cstripes[1] = d05_stride2; cstripes[2] = d05_stride3;
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    destroy_mat(&res);

    free(expS);

    SUITE_END();
}

void test_vaw_unifrac_weighted_normalized() {
    SUITE_START("test vaw weighted normalized unifrac");

#ifndef API_ONLY
    su::BPTree tree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table("test.biom");
#endif

    // as computed by GUniFrac, the original implementation of VAW-UniFrac
    // could not be found.
    //          Sample1   Sample2   Sample3   Sample4   Sample5   Sample6
    //Sample1 0.0000000 0.4086040 0.6240185 0.4639481 0.2857143 0.2766318
    //Sample2 0.4086040 0.0000000 0.3798594 0.6884992 0.6807616 0.4735781
    //Sample3 0.6240185 0.3798594 0.0000000 0.7713254 0.8812897 0.5047114
    //Sample4 0.4639481 0.6884992 0.7713254 0.0000000 0.6666667 0.2709298
    //Sample5 0.2857143 0.6807616 0.8812897 0.6666667 0.0000000 0.4735991
    //Sample6 0.2766318 0.4735781 0.5047114 0.2709298 0.4735991 0.0000000
    // weighted normalized unifrac as computed above

    std::vector<double*> w_exp;
    double w_stride1[] = {0.4086040, 0.3798594, 0.7713254, 0.6666667, 0.4735991, 0.2766318};
    double w_stride2[] = {0.6240185, 0.6884992, 0.8812897, 0.2709298, 0.2857143, 0.4735781};
    double w_stride3[] = {0.4639481, 0.6807616, 0.5047114, 0.4639481, 0.6807616, 0.5047114};
#ifndef API_ONLY
    w_exp.push_back(w_stride1);
    w_exp.push_back(w_stride2);
    w_exp.push_back(w_stride3);
    std::vector<double*> w_strides = su::make_strides(6);
    std::vector<double*> w_strides_total = su::make_strides(6);
    su::task_parameters w_task_p;
    w_task_p.start = 0; w_task_p.stop = 3; w_task_p.tid = 0; w_task_p.n_samples = 6; w_task_p.bypass_tips = false; w_task_p.normalize_sample_counts = true;
    w_task_p.g_unifrac_alpha = 1.0;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(w_task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::weighted_normalized,
                        true,
                        std::ref(w_strides),
                        std::ref(w_strides_total),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(w_strides[i][j] - w_exp[i][j]) < 0.000001);
        }
        free(w_strides[i]);
    }
#endif

    // repeat using the API
    double* cstripes[] = {w_stride1, w_stride2, w_stride3};
    mat_t *res = NULL;
    ComputeStatus rc = one_off("test.biom", "test.tre",
		               "weighted_normalized", true, 1.0,
			       false, 1, &res);
    ASSERT(rc == okay);
    double *expS = (double*)malloc(sizeof(double) * 15);
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    free(expS);
    destroy_mat(&res);

    SUITE_END();
}


#ifndef API_ONLY
void test_make_strides() {
    SUITE_START("test make stripes");
    std::vector<double*> exp;
    double stride[] = {0., 0., 0.};
    exp.push_back(stride);
    exp.push_back(stride);
    exp.push_back(stride);

    std::vector<double*> obs = su::make_strides(3);
    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(obs[i][j] - exp[i][j]) < 0.000001);
        }
        free(obs[i]);
    }
}
#endif

void test_faith_pd() {
    SUITE_START("test faith PD");

#ifndef API_ONLY
    // Note this tree is binary (opposed to example below)
    su::BPTree tree("((GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1):2,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table("test.biom");
#endif

    // make vector of expectations from faith PD
    double exp[6] = {6., 7., 8., 5., 4., 7.};

#ifndef API_ONLY
    // run faith PD to get obs
    double obs[6] = {0, 0, 0, 0, 0, 0};

    su::faith_pd(table, tree, obs);

    // ASSERT that results = expectation
    for (unsigned int i = 0; i < 6; i++){
        ASSERT(fabs(exp[i]-obs[i]) < 0.000001)
    }
#endif

    // repeat using the API
    r_vec *res = NULL;
    ComputeStatus rc = faith_pd_one_off("test.biom", "test_faith.tre", &res);
    ASSERT(rc == okay);
    for(unsigned int i = 0; i < 6; i++) {
       ASSERT(fabs(res->values[i] - exp[i]) < 0.000001);
    }
    destroy_results_vec(&res);

    SUITE_END();
}

// Synthetic 5-OTU / 6-sample CSR fixture shared by the in-memory API tests
// below. Identical to test/capi_inmem_test.c (table data) and the multi-
// furcating tree from src/tests/wasm/fixtures.hpp (unit branch lengths).
namespace inmem_fixture {
    static const unsigned int N_OBS  = 5;
    static const unsigned int N_SAMP = 6;
    static const char* const  OBS_IDS[]  = {"GG_OTU_1","GG_OTU_2","GG_OTU_3","GG_OTU_4","GG_OTU_5"};
    static const char* const  SAMP_IDS[] = {"Sample1","Sample2","Sample3","Sample4","Sample5","Sample6"};
    static       uint32_t     INDICES[]  = {2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2};
    static       uint32_t     INDPTR[]   = {0, 1, 6, 9, 13, 15};
    static       double       DATA[]     = {1., 5., 1., 2., 3., 1., 1., 4., 2., 2., 1., 1., 1., 1., 1.};
    static const unsigned int NPARENS = 16;
    static       bool         STRUCTURE[] = { true, true, false, true,
                                              true, false, true, false,
                                              false, true, true, false,
                                              true, false, false, false };
    static       double       LENGTHS[]   = { 0., 1., 0., 1.,
                                              1., 0., 1., 0.,
                                              0., 1., 1., 0.,
                                              1., 0., 0., 0. };
    static const char* const  NAMES[]     = {"", "GG_OTU_1", "", "",
                                             "GG_OTU_2", "", "GG_OTU_3", "",
                                             "", "", "GG_OTU_5", "",
                                             "GG_OTU_4", "", "", ""};
    static const uint32_t     GROUPING[6] = {0, 0, 1, 1, 1, 0};

    static support_biom_t make_table() {
        support_biom_t t = {(char**) OBS_IDS, (char**) SAMP_IDS,
                            INDICES, INDPTR, DATA, (int) N_OBS, (int) N_SAMP, 0};
        return t;
    }

    static support_bptree_t make_tree() {
        support_bptree_t t = {STRUCTURE, LENGTHS, (char**) NAMES, (int) NPARENS};
        return t;
    }

    /* Subsampling depth for the seeded cases. Every sample in the fixture has a
     * total count >= this, so none are dropped and the matrix keeps its size.
     */
    static const unsigned int SUBSAMPLE_DEPTH = 3;
    static const int          SUBSAMPLE_SEED  = 42;

    /* One unweighted compute, flattened. seed < 0 takes the unsubsampled v3
     * path; seed >= 0 the subsampled v4 path.
     */
    static ComputeStatus run_matrix(std::vector<float> &out,
                                    unsigned int n_substeps, int seed) {
        const support_biom_t   table = make_table();
        const support_bptree_t tree  = make_tree();

        mat_full_fp32_t* mat = NULL;
        ComputeStatus rc = (seed < 0)
            ? one_off_matrix_inmem_fp32_v3(&table, &tree, "unweighted_fp32",
                                           false, 1.0, false, true, n_substeps,
                                           0, false, NULL, &mat)
            : one_off_matrix_inmem_fp32_v4(&table, &tree, "unweighted_fp32",
                                           false, 1.0, false, true, n_substeps,
                                           SUBSAMPLE_DEPTH, false, seed, NULL, &mat);
        if (rc != okay) return rc;

        const size_t n_els = size_t(mat->n_samples) * size_t(mat->n_samples);
        out.assign(mat->matrix, mat->matrix + n_els);
        destroy_mat_full_fp32(&mat);
        return okay;
    }
}

void test_faith_pd_inmem() {
    SUITE_START("test faith_pd_inmem");

    using namespace inmem_fixture;
    const support_biom_t   table = {(char**) OBS_IDS, (char**) SAMP_IDS,
                                    INDICES, INDPTR, DATA, (int) N_OBS, (int) N_SAMP, 0};
    const support_bptree_t tree  = {STRUCTURE, LENGTHS, (char**) NAMES, (int) NPARENS};

    // Hand-derived from the multifurcating tree (root → {OTU_1, {OTU_2,OTU_3},
    // {OTU_5,OTU_4}}, all branches = 1.0) and the per-sample OTU sets
    // implied by the CSR table.
    double exp[6] = {4., 5., 6., 3., 2., 5.};

    r_vec* res = NULL;
    ComputeStatus rc = faith_pd_inmem(&table, &tree, &res);
    ASSERT(rc == okay);
    ASSERT(res != NULL);
    ASSERT(res->n_samples == N_SAMP);
    for (unsigned int i = 0; i < N_SAMP; i++) {
        ASSERT(fabs(res->values[i] - exp[i]) < 1e-6);
    }
    destroy_results_vec(&res);

    // Error paths.
    r_vec* tmp = NULL;
    ASSERT(faith_pd_inmem(&table, NULL, &tmp) == tree_missing);
    ASSERT(faith_pd_inmem(NULL, &tree, &tmp) == table_missing);

    SUITE_END();
}

void test_subsample_inmem() {
    SUITE_START("test subsample_table_inmem + accessors");

    using namespace inmem_fixture;
    const support_biom_t table = {(char**) OBS_IDS, (char**) SAMP_IDS,
                                  INDICES, INDPTR, DATA, (int) N_OBS, (int) N_SAMP, 0};
    const unsigned int depth = 3;

    ssu_set_random_seed(42);
    opaque_biom_inmem_t* sub = NULL;
    ASSERT(subsample_table_inmem(&table, depth, false, &sub) == okay);
    ASSERT(sub != NULL);

    // All six per-sample counts ({7,3,4,6,3,3}) ≥ depth=3, so every sample
    // survives. Each retained column sums to depth.
    ASSERT(subsampled_n_samples(sub) == N_SAMP);
    unsigned int n_sub_obs  = subsampled_n_obs(sub);
    unsigned int n_sub_samp = subsampled_n_samples(sub);
    double col_sums[6] = {0., 0., 0., 0., 0., 0.};
    for (unsigned int i = 0; i < n_sub_obs; i++) {
        const char* oid = subsampled_get_obs_id(sub, i);
        ASSERT(oid != NULL);
        double row[6] = {0., 0., 0., 0., 0., 0.};
        ASSERT(subsampled_get_obs_data(sub, oid, row));
        for (unsigned int j = 0; j < n_sub_samp; j++) col_sums[j] += row[j];
    }
    for (unsigned int j = 0; j < n_sub_samp; j++) {
        ASSERT(fabs(col_sums[j] - (double) depth) < 1e-9);
    }

    // Out-of-range / unknown-id rejection.
    ASSERT(subsampled_get_sample_id(sub, n_sub_samp) == NULL);
    double junk[6] = {0., 0., 0., 0., 0., 0.};
    ASSERT(!subsampled_get_obs_data(sub, "NOT_A_REAL_OTU", junk));

    // Determinism: identical seed → identical cells.
    ssu_set_random_seed(42);
    opaque_biom_inmem_t* sub2 = NULL;
    ASSERT(subsample_table_inmem(&table, depth, false, &sub2) == okay);
    ASSERT(subsampled_n_obs(sub2) == n_sub_obs);
    for (unsigned int i = 0; i < n_sub_obs; i++) {
        const char* oid = subsampled_get_obs_id(sub, i);
        double r1[6] = {0.}, r2[6] = {0.};
        subsampled_get_obs_data(sub,  oid, r1);
        subsampled_get_obs_data(sub2, oid, r2);
        for (unsigned int j = 0; j < n_sub_samp; j++) ASSERT(r1[j] == r2[j]);
    }

    destroy_subsampled_inmem(&sub);
    ASSERT(sub == NULL);
    destroy_subsampled_inmem(&sub2);

    SUITE_END();
}

void test_permanova_inmem() {
    SUITE_START("test compute_permanova_inmem_fp64/fp32");

    using namespace inmem_fixture;
    const support_biom_t   table = {(char**) OBS_IDS, (char**) SAMP_IDS,
                                    INDICES, INDPTR, DATA, (int) N_OBS, (int) N_SAMP, 0};
    // Unit branch lengths so unifrac produces a non-degenerate distance
    // matrix; an all-zero-lengths tree drives PERMANOVA into a degenerate
    // fstat regime.
    const support_bptree_t tree  = {STRUCTURE, LENGTHS, (char**) NAMES, (int) NPARENS};

    mat_full_fp64_t* dm = NULL;
    ComputeStatus rc = one_off_matrix_inmem_v2(&table, &tree, "unweighted_fp64",
                                               false, 1.0, false, 1,
                                               0, true, NULL, &dm);
    ASSERT(rc == okay);
    ASSERT(dm != NULL);

    ssu_set_random_seed(42);
    double fstat = 0.0, pvalue = 0.0;
    rc = compute_permanova_inmem_fp64(dm->matrix, dm->n_samples, GROUPING,
                                      999, &fstat, &pvalue);
    ASSERT(rc == okay);
    ASSERT(fstat > 0.0);
    ASSERT(pvalue > 0.0 && pvalue <= 1.0);

    // Re-seed and rerun. Under multi-threaded OMP the unpermuted F is
    // computed via parallel reductions and can drift by ULPs across runs;
    // the pvalue can shift correspondingly when the observed F sits near
    // a permutation-tail boundary. Allow a tight numeric tolerance rather
    // than requiring bit-exactness.
    ssu_set_random_seed(42);
    double fstat2 = 0.0, pvalue2 = 0.0;
    compute_permanova_inmem_fp64(dm->matrix, dm->n_samples, GROUPING,
                                 999, &fstat2, &pvalue2);
    ASSERT(fabs(fstat2 - fstat) < 1e-6);
    ASSERT(fabs(pvalue2 - pvalue) < 1e-2);

    // Error paths.
    ASSERT(compute_permanova_inmem_fp64(NULL, dm->n_samples, GROUPING, 9, &fstat, &pvalue) != okay);
    ASSERT(compute_permanova_inmem_fp64(dm->matrix, dm->n_samples, NULL, 9, &fstat, &pvalue) != okay);

    // fp32 variant: same matrix in fp32 layout.
    mat_full_fp32_t* dm32 = NULL;
    rc = one_off_matrix_inmem_fp32_v2(&table, &tree, "unweighted_fp32",
                                      false, 1.0, false, 1,
                                      0, true, NULL, &dm32);
    ASSERT(rc == okay);
    ssu_set_random_seed(42);
    float fstat32 = 0.0f, pvalue32 = 0.0f;
    rc = compute_permanova_inmem_fp32(dm32->matrix, dm32->n_samples, GROUPING,
                                      999, &fstat32, &pvalue32);
    ASSERT(rc == okay);
    ASSERT(fstat32 > 0.0f);
    ASSERT(pvalue32 > 0.0f && pvalue32 <= 1.0f);

    destroy_mat_full_fp64(&dm);
    destroy_mat_full_fp32(&dm32);

    SUITE_END();
}

/* n_substeps arrives straight from the caller and only says how to split the
 * stripe range across tasks, so every value must produce the same matrix rather
 * than a crash. The 6-sample fixture has (6 + 1) / 2 = 3 stripes, so the cases
 * below cover fewer substeps than stripes, exactly as many, more, and zero.
 */
void test_matrix_inmem_substeps() {
    SUITE_START("test one_off_matrix_inmem_fp32 substep bounds");

    using namespace inmem_fixture;

    // n_substeps == 1 comes first and establishes the expected matrix
    const unsigned int cases[] = {1, 2, 3, 4, 8, 64, 0};
    const unsigned int n_cases = sizeof(cases) / sizeof(cases[0]);

    std::vector<float> reference;
    for (unsigned int c = 0; c < n_cases; c++) {
        std::vector<float> got;
        ASSERT(run_matrix(got, cases[c], /*seed*/ -1) == okay);
        ASSERT(got.size() == size_t(N_SAMP) * size_t(N_SAMP));
        if (reference.empty()) {
            reference = got;
        } else {
            // splitting the same stripes over more tasks changes nothing about
            // the arithmetic within a stripe, so this is exact
            ASSERT(got == reference);
        }
    }

    SUITE_END();
}

/* Same bounds question for partial_v3, which differs from the one_off entries in
 * a way that matters: it sizes dm_stripes to the *total* stripe count while
 * asking set_tasks to divide only the caller's sub-range. So the number of
 * stripes the tasks actually divide is stripe_stop - stripe_start, and a
 * sub-range that does not start at stripe 0 is the case where clamping against
 * the total instead of the sub-range still handed a trailing task a start index
 * one past the end of dm_stripes.
 *
 * test.biom has 6 samples, so (6 + 1) / 2 = 3 stripes total and a sub-range of
 * 2 here. Note the overflow this pins is an out-of-bounds *read* of a pointer
 * that an empty stripe range never dereferences, so it needs a sanitizer to be
 * seen -- this test passing is necessary but not sufficient.
 */
void test_partial_substeps() {
    SUITE_START("test partial_v3 substep bounds");

    const unsigned int stripe_start = 1;
    const unsigned int stripe_stop  = 3;
    const unsigned int cases[] = {1, 2, 3, 8, 0};
    const unsigned int n_cases = sizeof(cases) / sizeof(cases[0]);

    std::vector<std::vector<double> > reference;
    for (unsigned int c = 0; c < n_cases; c++) {
        partial_mat_t* pm = NULL;
        ASSERT(partial_v3("test.biom", "test.tre", "unweighted",
                          false, 1.0, false, true,
                          cases[c],                    // n_substeps
                          stripe_start, stripe_stop, &pm) == okay);
        ASSERT(pm != NULL);
        ASSERT(pm->stripe_start == stripe_start);
        ASSERT(pm->stripe_stop == stripe_stop);

        std::vector<std::vector<double> > got;
        for (unsigned int s = 0; s < stripe_stop - stripe_start; s++)
            got.push_back(std::vector<double>(pm->stripes[s],
                                              pm->stripes[s] + pm->n_samples));

        if (reference.empty()) {
            reference = got;
        } else {
            ASSERT(got == reference);
        }
        destroy_partial_mat(&pm);
    }

    SUITE_END();
}

/* Concurrency: several computes in flight in one process must not corrupt
 * shared library state, and each must return the same answer it would have
 * returned on its own.
 *
 * The ASSERT macros mutate non-atomic harness counters, so worker threads never
 * assert; each records into its own slot and the main thread checks every slot
 * after joining.
 */
namespace concurrency_fixture {
    static const unsigned int N_THREADS = 4;
    static const unsigned int N_ITERS   = 25;

    struct outcome {
        unsigned int n_ok       = 0;  // computes that returned okay
        unsigned int n_status   = 0;  // computes that returned something else
        unsigned int n_mismatch = 0;  // computes that disagreed with the reference
    };

    static void check(const std::vector<outcome> &results) {
        for (unsigned int t = 0; t < results.size(); t++) {
            ASSERT(results[t].n_status == 0);
            ASSERT(results[t].n_mismatch == 0);
            ASSERT(results[t].n_ok == N_ITERS);
        }
    }

    // run one worker per thread, join, then check every slot
    template<typename F>
    static void run_workers(F worker) {
        std::vector<outcome> results(N_THREADS);
        std::vector<std::thread> workers;
        for (unsigned int t = 0; t < N_THREADS; t++)
            workers.emplace_back(worker, &results[t]);
        for (unsigned int t = 0; t < N_THREADS; t++)
            workers[t].join();
        check(results);
    }

    /* The unweighted compute is deterministic, so a concurrent result must be
     * bit-identical to the serial one -- no tolerance needed. seed selects the
     * plain or the subsampled path, as in inmem_fixture::run_matrix.
     */
    static void matrix_worker(int seed, const std::vector<float>* reference, outcome* out) {
        for (unsigned int i = 0; i < N_ITERS; i++) {
            std::vector<float> got;
            if (inmem_fixture::run_matrix(got, 1, seed) != okay) {
                out->n_status++;
                continue;
            }
            if (got != *reference) out->n_mismatch++;
            out->n_ok++;
        }
    }

    static void faith_pd_worker(outcome* out) {
        using namespace inmem_fixture;
        const support_biom_t   table = make_table();
        const support_bptree_t tree  = make_tree();
        // same expectation as test_faith_pd_inmem
        const double expected[6] = {4., 5., 6., 3., 2., 5.};

        for (unsigned int i = 0; i < N_ITERS; i++) {
            r_vec* res = NULL;
            if (faith_pd_inmem(&table, &tree, &res) != okay) {
                out->n_status++;
                continue;
            }
            if (res->n_samples != N_SAMP) {
                out->n_mismatch++;
            } else {
                for (unsigned int j = 0; j < N_SAMP; j++) {
                    if (fabs(res->values[j] - expected[j]) > 1e-6) {
                        out->n_mismatch++;
                        break;
                    }
                }
            }
            destroy_results_vec(&res);
            out->n_ok++;
        }
    }

    /* Serial compute: both the expected answer and the thing that installs the
     * SIGUSR1 handler, since su::process_stripes registers it.
     */
    static void serial_reference(std::vector<float> &reference) {
        ASSERT(inmem_fixture::run_matrix(reference, 1, /*seed*/ -1) == okay);
        ASSERT(reference.size() == size_t(inmem_fixture::N_SAMP) * size_t(inmem_fixture::N_SAMP));
    }
}

void test_concurrent_matrix_inmem() {
    SUITE_START("test concurrent one_off_matrix_inmem_fp32");

    using namespace concurrency_fixture;

    // Serial reference first, so the expected answer is known-good.
    std::vector<float> reference;
    serial_reference(reference);

    run_workers([&reference](outcome* o) { matrix_worker(-1, &reference, o); });

    SUITE_END();
}

void test_concurrent_faith_pd_inmem() {
    SUITE_START("test concurrent faith_pd_inmem");

    using namespace concurrency_fixture;

    run_workers(faith_pd_worker);

    SUITE_END();
}

/* PCoA and PERMANOVA also draw from an RNG -- the randomized SVD needs a random
 * matrix, and PERMANOVA needs its permutations. Without a per-call seed both
 * draw from the dependency's process-global generator, which no caller can hold
 * still while another one runs. The seeded entry points take a generator local
 * to the call instead, which is what these tests pin.
 *
 * Both are checked to a tolerance rather than bit-exactly -- see "Ordination
 * reproduces to a tolerance" in README.md. The bounds are far below the signal:
 * a seed that was not local to the call would move the answer by ~0.5, the
 * scale test_pcoa_seeded pins directly.
 */
namespace concurrency_fixture {
    static const int          ORD_SEED   = 7;
    static const unsigned int PCOA_DIMS  = 3;   // < n_samples (6)
    static const unsigned int PERM_PERMS = 99;
    static const double       PCOA_TOL   = 1e-12;
    static const double       FSTAT_TOL  = 1e-6;
    static const double       PVALUE_TOL = 1e-2;

    // the distance matrix both ordination tests run on
    static mat_full_fp64_t* ordination_dm() {
        using namespace inmem_fixture;
        const support_biom_t   table = make_table();
        const support_bptree_t tree  = make_tree();

        mat_full_fp64_t* dm = NULL;
        ASSERT(one_off_matrix_inmem_v4(&table, &tree, "unweighted_fp64",
                                       false, 1.0, false, true, 1, 0, false, -1,
                                       NULL, &dm) == okay);
        ASSERT(dm != NULL);
        return dm;
    }

    static bool collect_pcoa(const double *dm, unsigned int n_samples,
                             std::vector<double> &out) {
        double *ev = NULL, *sa = NULL, *pe = NULL;
        pcoa_seeded(dm, n_samples, PCOA_DIMS, ORD_SEED, &ev, &sa, &pe);
        if (ev == NULL || sa == NULL || pe == NULL) return false;

        out.clear();
        out.insert(out.end(), ev, ev + PCOA_DIMS);
        out.insert(out.end(), sa, sa + (size_t(PCOA_DIMS) * n_samples));
        out.insert(out.end(), pe, pe + PCOA_DIMS);
        free(ev);
        free(sa);
        free(pe);
        return true;
    }

    static void pcoa_worker(const double *dm, unsigned int n_samples,
                            const std::vector<double>* reference, outcome* out) {
        for (unsigned int i = 0; i < N_ITERS; i++) {
            std::vector<double> got;
            if (!collect_pcoa(dm, n_samples, got)) {
                out->n_status++;
                continue;
            }
            if (got.size() != reference->size()) {
                out->n_mismatch++;
            } else {
                for (size_t j = 0; j < got.size(); j++) {
                    if (fabs(got[j] - (*reference)[j]) > PCOA_TOL) {
                        out->n_mismatch++;
                        break;
                    }
                }
            }
            out->n_ok++;
        }
    }

    static void permanova_worker(const double *dm, unsigned int n_samples,
                                 double ref_fstat, double ref_pvalue, outcome* out) {
        for (unsigned int i = 0; i < N_ITERS; i++) {
            double fstat = 0.0, pvalue = 0.0;
            if (compute_permanova_inmem_fp64_seeded(dm, n_samples, inmem_fixture::GROUPING,
                                                    PERM_PERMS, ORD_SEED,
                                                    &fstat, &pvalue) != okay) {
                out->n_status++;
                continue;
            }
            if (fabs(fstat - ref_fstat) > FSTAT_TOL || fabs(pvalue - ref_pvalue) > PVALUE_TOL)
                out->n_mismatch++;
            out->n_ok++;
        }
    }
}

void test_concurrent_pcoa() {
    SUITE_START("test concurrent pcoa_seeded");

    using namespace concurrency_fixture;

    mat_full_fp64_t* dm = ordination_dm();

    // serial reference at the same seed
    std::vector<double> reference;
    ASSERT(collect_pcoa(dm->matrix, dm->n_samples, reference));

    run_workers([&](outcome* o) { pcoa_worker(dm->matrix, dm->n_samples, &reference, o); });

    destroy_mat_full_fp64(&dm);

    SUITE_END();
}

void test_concurrent_permanova_inmem() {
    SUITE_START("test concurrent compute_permanova_inmem_fp64_seeded");

    using namespace concurrency_fixture;

    mat_full_fp64_t* dm = ordination_dm();

    // serial reference at the same seed
    double ref_fstat = 0.0, ref_pvalue = 0.0;
    ASSERT(compute_permanova_inmem_fp64_seeded(dm->matrix, dm->n_samples,
                                               inmem_fixture::GROUPING,
                                               PERM_PERMS, ORD_SEED,
                                               &ref_fstat, &ref_pvalue) == okay);
    ASSERT(ref_fstat > 0.0);
    ASSERT(ref_pvalue > 0.0 && ref_pvalue <= 1.0);

    run_workers([&](outcome* o) {
        permanova_worker(dm->matrix, dm->n_samples, ref_fstat, ref_pvalue, o);
    });

    destroy_mat_full_fp64(&dm);

    SUITE_END();
}

/* A subsampled compute draws from an RNG. Passing the seed per call is what
 * makes a reproducible subsampled compute possible without holding a lock
 * across seed-then-compute: the alternative, ssu_set_random_seed() followed by a
 * v3 call, mutates a process-global generator, so concurrent callers interleave
 * their seeding and neither gets the answer it asked for.
 *
 * Bit-exactness here relies on every call seeing the same OpenMP width, since
 * the subsample draw is distributed across the OpenMP team. That holds within
 * one process: nthreads-var is a per-thread ICV and no one changes it, so the
 * plain std::threads below each get a team of the same size. It is NOT a claim
 * that a subsampled result is reproducible across different thread counts.
 */
void test_concurrent_matrix_inmem_seeded() {
    SUITE_START("test concurrent seeded one_off_matrix_inmem_fp32_v4");

    using namespace inmem_fixture;
    using namespace concurrency_fixture;

    // serial reference at the same seed
    std::vector<float> reference;
    ASSERT(run_matrix(reference, 1, SUBSAMPLE_SEED) == okay);
    ASSERT(reference.size() == size_t(N_SAMP) * size_t(N_SAMP));

    run_workers([&reference](outcome* o) { matrix_worker(SUBSAMPLE_SEED, &reference, o); });

    SUITE_END();
}

/* v4 seeding semantics, serially: an explicit seed is reproducible, and a
 * negative seed keeps the legacy behaviour of drawing from the process-global
 * generator that ssu_set_random_seed() sets.
 */
void test_matrix_inmem_seeded() {
    SUITE_START("test one_off_matrix_inmem_fp32_v4 seeding");

    using namespace inmem_fixture;

    std::vector<float> a, b, c, v3, v4_neg;

    // an explicit seed reproduces, with no seeding call in between
    ASSERT(run_matrix(a, 1, 7) == okay);
    ASSERT(run_matrix(b, 1, 7) == okay);
    ASSERT(a == b);

    // ... and does not depend on the global generator's state
    ssu_set_random_seed(999);
    ASSERT(run_matrix(c, 1, 7) == okay);
    ASSERT(a == c);

    /* A negative seed is the legacy path: same global seed, same answer as v3.
     * run_matrix takes the v3 entry point at seed < 0 without subsampling, so
     * call the two explicitly here -- this is the one place the equivalence of
     * the two entry points is what is under test.
     */
    const support_biom_t   table = make_table();
    const support_bptree_t tree  = make_tree();
    for (int pass = 0; pass < 2; pass++) {
        mat_full_fp32_t* mat = NULL;
        ssu_set_random_seed(SUBSAMPLE_SEED);
        ASSERT((pass == 0
                ? one_off_matrix_inmem_fp32_v3(&table, &tree, "unweighted_fp32", false, 1.0,
                                               false, true, 1, SUBSAMPLE_DEPTH, false,
                                               NULL, &mat)
                : one_off_matrix_inmem_fp32_v4(&table, &tree, "unweighted_fp32", false, 1.0,
                                               false, true, 1, SUBSAMPLE_DEPTH, false,
                                               -1, NULL, &mat)) == okay);
        ASSERT(mat != NULL);
        const size_t n_els = size_t(mat->n_samples) * size_t(mat->n_samples);
        (pass == 0 ? v3 : v4_neg).assign(mat->matrix, mat->matrix + n_els);
        destroy_mat_full_fp32(&mat);
    }
    ASSERT(v3 == v4_neg);

    SUITE_END();
}

/* The progress-reporting path itself, under concurrency: SIGUSR1 sets every
 * flag, and each in-flight compute clears and reports its own via sync_printf.
 *
 * The flags are raised *before* the workers start rather than during the run,
 * so every worker is guaranteed to hit one instead of racing the signal against
 * a compute that may already be finished -- deterministic coverage rather than
 * a test that usually exercises nothing.
 *
 * Emits a few "tid:..." progress lines on stdout; that is the feature working.
 * Keep this last in main(): sig_handler sets all CPU_SETSIZE flags and only the
 * ones belonging to tasks that actually run get consumed, so the leftovers
 * would make later suites emit stray progress lines.
 */
void test_concurrent_matrix_inmem_reporting() {
    SUITE_START("test concurrent progress reporting");

    using namespace concurrency_fixture;

    // Also installs the SIGUSR1 handler, so the raise() below cannot hit the
    // default disposition and kill the test process.
    std::vector<float> reference;
    serial_reference(reference);

    raise(SIGUSR1);

    // Reporting must not disturb the results, and must not crash.
    run_workers([&reference](outcome* o) { matrix_worker(-1, &reference, o); });

    SUITE_END();
}

void test_faith_pd_shear(){
    SUITE_START("test faith PD extra OTUs in tree");

#ifndef API_ONLY
    su::BPTree tree("((GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1,GG_OTU_ex:9):1):2,(GG_OTU_5:1,GG_OTU_4:1,GG_OTU_ex2:12):1);");
    su::biom table("test.biom");
#endif

    // make vector of expectations from faith PD
    double exp[6] = {6., 7., 8., 5., 4., 7.};

#ifndef API_ONLY
    // run faith PD to get obs
    double obs[6] = {0, 0, 0, 0, 0, 0};

    std::unordered_set<std::string> to_keep(table.get_obs_ids().begin(),           \
                                            table.get_obs_ids().end());            \
    su::BPTree tree_sheared = tree.shear(to_keep).collapse();
    su::faith_pd(table, tree_sheared, obs);

    // ASSERT that results = expectation
    for (unsigned int i = 0; i < 6; i++){
        ASSERT(fabs(exp[i]-obs[i]) < 0.000001)
    }
#endif

    // repeat using the API
    r_vec *res = NULL;
    ComputeStatus rc = faith_pd_one_off("test.biom", "test_faith_shear.tre", &res);
    ASSERT(rc == okay);
    for(unsigned int i = 0; i < 6; i++) {
       ASSERT(fabs(res->values[i] - exp[i]) < 0.000001);
    }
    destroy_results_vec(&res);

    SUITE_END();
}

void test_unweighted_unifrac() {
    SUITE_START("test unweighted unifrac");
#ifndef API_ONLY
    su::BPTree tree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table("test.biom");
#endif

    std::vector<double*> exp;
    double stride1[] = {0.2, 0.42857143, 0.71428571, 0.33333333, 0.6, 0.2};
    double stride2[] = {0.57142857, 0.66666667, 0.85714286, 0.4, 0.5, 0.33333333};
    double stride3[] = {0.6, 0.6, 0.42857143, 0.6, 0.6, 0.42857143};
#ifndef API_ONLY
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = false; task_p.normalize_sample_counts = true;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::unweighted,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
    }
#endif

    // repeat using the API
    double* cstripes[] = {stride1, stride2, stride3};
    mat_t *res = NULL;
    ComputeStatus rc = one_off("test.biom", "test.tre",
		               "unweighted", false, 0.0,
			       false, 1, &res);
    ASSERT(rc == okay);
    double *expS = (double*)malloc(sizeof(double) * 15);
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    free(expS);
    destroy_mat(&res);

    SUITE_END();
}

void test_unweighted_unifrac_fast() {
    SUITE_START("test unweighted unifrac no tips");
#ifndef API_ONLY
    su::BPTree tree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table("test.biom");
#endif

    std::vector<double*> exp;
    double stride1[] = {0., 0., 0.5, 0., 0.5, 0.};
    double stride2[] = {0., 0.5, 0.5, 0.5, 0.5, 0.};
    double stride3[] = {0.5, 0.5, 0., 0.5, 0.5, 0.};
#ifndef API_ONLY
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = true; task_p.normalize_sample_counts = true;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::unweighted,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
    }
#endif

    // repeat using the API
    double* cstripes[] = {stride1, stride2, stride3};
    mat_t *res = NULL;
    ComputeStatus rc = one_off("test.biom", "test.tre",
		               "unweighted", false, 0.0,
			       true, 1, &res);
    ASSERT(rc == okay);
    double *expS = (double*)malloc(sizeof(double) * 15);
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    free(expS);
    destroy_mat(&res);

    SUITE_END();
}

void test_unnormalized_unweighted_unifrac() {
    SUITE_START("test unnormalized unweighted unifrac");
#ifndef API_ONLY
    su::BPTree tree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table("test.biom");
#endif

    std::vector<double*> exp;
    double stride1[] = {1,3,5,1,3,1};
    double stride2[] = {4,4,6,2,2,2};
    double stride3[] = {3,3,3,3,3,3};
#ifndef API_ONLY
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = false; task_p.normalize_sample_counts = true;

    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::unweighted_unnormalized,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
    }
#endif

    // repeat using the API
    double* cstripes[] = {stride1, stride2, stride3};
    mat_t *res = NULL;
    ComputeStatus rc = one_off("test.biom", "test.tre",
		               "unweighted_unnormalized", false, 0.0,
			       false, 1, &res);
    ASSERT(rc == okay);
    double *expS = (double*)malloc(sizeof(double) * 15);
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    free(expS);
    destroy_mat(&res);

    SUITE_END();
}

void test_normalized_weighted_unifrac() {
    SUITE_START("test normalized weighted unifrac");
#ifndef API_ONLY
    su::BPTree tree("(GG_OTU_1:1,(GG_OTU_2:1,GG_OTU_3:1):1,(GG_OTU_5:1,GG_OTU_4:1):1);");
    su::biom table("test.biom");
#endif

    std::vector<double*> exp;
    double stride1[] = {0.38095238, 0.33333333, 0.73333333, 0.33333333, 0.5, 0.26785714};
    double stride2[] = {0.58095238, 0.66666667, 0.86666667, 0.25, 0.28571429, 0.45833333};
    double stride3[] = {0.47619048, 0.66666667, 0.46666667, 0.47619048, 0.66666667, 0.46666667};
#ifndef API_ONLY
    exp.push_back(stride1);
    exp.push_back(stride2);
    exp.push_back(stride3);
    std::vector<double*> strides = su::make_strides(6);
    std::vector<double*> strides_total = su::make_strides(6);

    su::task_parameters task_p;
    task_p.start = 0; task_p.stop = 3; task_p.tid = 0; task_p.n_samples = 6; task_p.bypass_tips = false; task_p.normalize_sample_counts = true;


    std::vector<su::task_parameters> tasks;
    tasks.push_back(task_p);
    su::process_stripes(std::ref(table), 
                        std::ref(tree),
                        su::weighted_normalized,
                        false,
                        std::ref(strides),
                        std::ref(strides_total),
                        std::ref(tasks));

    for(unsigned int i = 0; i < 3; i++) {
        for(unsigned int j = 0; j < 6; j++) {
            ASSERT(fabs(strides[i][j] - exp[i][j]) < 0.000001);
        }
        free(strides[i]);
    }
#endif

    // repeat using the API
    double* cstripes[] = {stride1, stride2, stride3};
    mat_t *res = NULL;
    ComputeStatus rc = one_off("test.biom", "test.tre",
		               "weighted_normalized", false, 0.0,
			       false, 1, &res);
    ASSERT(rc == okay);
    double *expS = (double*)malloc(sizeof(double) * 15);
    _testv_stripes_to_condensed_form(cstripes, 6, 3, expS);
    for(unsigned int i = 0; i < 15; i++) {
       ASSERT(fabs(res->condensed_form[i] - expS[i]) < 0.000001);
    }
    free(expS);
    destroy_mat(&res);

    SUITE_END();
}

#ifndef API_ONLY
void test_bptree_shear_simple() {
    SUITE_START("test bptree shear simple");
    su::BPTree tree("((3:2,4:3,(6:5)5:4)2:1,7:6,((10:9,11:10)9:8)8:7)r");

    // simple
    std::unordered_set<std::string> to_keep = {"4", "6", "7", "10", "11"};

    uint32_t exp_nparens = 20;
    std::vector<bool> exp_structure = {true, true, true, false, true, true, false, false, false, true,
                                       false, true, true, true, false, true, false, false, false, false};
    std::vector<std::string> exp_names = {"r", "2", "4", "", "5", "6", "", "", "", "7", "", "8", "9", "10", "",
                                          "11", "", "", "", ""};
    std::vector<double> exp_lengths = {0, 1, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0, 10, 0, 0, 0, 0};

    su::BPTree obs = tree.shear(to_keep);
    ASSERT(obs.get_structure() == exp_structure);
    ASSERT(exp_nparens == obs.nparens);
    ASSERT(vec_almost_equal(exp_lengths, obs.lengths));
    ASSERT(obs.names == exp_names);
    SUITE_END();
}

void test_bptree_shear_deep() {
    SUITE_START("test bptree shear deep");
    su::BPTree tree("((3:2,4:3,(6:5)5:4)2:1,7:6,((10:9,11:10)9:8)8:7)r");

    // deep
    std::unordered_set<std::string> to_keep = {"10", "11"};

    uint32_t exp_nparens = 10;
    std::vector<bool> exp_structure = {true, true, true, true, false, true, false, false, false, false};
    std::vector<std::string> exp_names = {"r", "8", "9", "10", "", "11", "", "", "", ""};
    std::vector<double> exp_lengths = {0, 7, 8, 9, 0, 10, 0, 0, 0, 0};

    su::BPTree obs = tree.shear(to_keep);
    ASSERT(exp_nparens == obs.nparens);
    ASSERT(obs.get_structure() == exp_structure);
    ASSERT(vec_almost_equal(exp_lengths, obs.lengths));
    ASSERT(obs.names == exp_names);
    SUITE_END();
}

void test_test_table_ids_are_subset_of_tree() {
    SUITE_START("test test_table_ids_are_subset_of_tree");

    su::BPTree tree("(a:1,b:2)r;");
    su::biom table("test.biom");
    std::string expected = "GG_OTU_1";
    std::string observed = su::test_table_ids_are_subset_of_tree(table, tree);
    ASSERT(observed == expected);

    su::BPTree tree2("(GG_OTU_1,GG_OTU_5,GG_OTU_6,GG_OTU_2,GG_OTU_3,GG_OTU_4);");
    su::biom table2("test.biom");
    expected = "";
    observed = su::test_table_ids_are_subset_of_tree(table2, tree2);
    ASSERT(observed == expected);
    SUITE_END();
}


void test_bptree_get_tip_names() {
    SUITE_START("test bptree get_tip_names");
    su::BPTree tree("((a:2,b:3,(c:5)d:4)e:1,f:6,((g:9,h:10)i:8)j:7)r");

    std::unordered_set<std::string> expected = {"a", "b", "c", "f", "g", "h"};
    std::unordered_set<std::string> observed = tree.get_tip_names();
    ASSERT(observed == expected);
    SUITE_END();
}

void test_bptree_collapse_simple() {
    SUITE_START("test bptree collapse simple");
    su::BPTree tree("((3:2,4:3,(6:5)5:4)2:1,7:6,((10:9,11:10)9:8)8:7)r");

    uint32_t exp_nparens = 18;
    std::vector<bool> exp_structure = {true, true, true, false, true, false, true, false, false,
                                       true, false, true, true, false, true, false, false, false};
    std::vector<std::string> exp_names = {"r", "2", "3", "", "4", "", "6", "", "", "7", "", "9", "10", "", "11", "", "", ""};
    std::vector<double> exp_lengths = {0, 1, 2, 0, 3, 0, 9, 0, 0, 6, 0, 15, 9, 0, 10, 0, 0, 0};

    su::BPTree obs = tree.collapse();

    ASSERT(obs.get_structure() == exp_structure);
    ASSERT(exp_nparens == obs.nparens);
    ASSERT(vec_almost_equal(exp_lengths, obs.lengths));
    ASSERT(obs.names == exp_names);
    SUITE_END();
}

void test_bptree_collapse_edge() {
    SUITE_START("test bptree collapse edge case against root");

    su::BPTree tree("((a),b)r;");
    su::BPTree exp("(a,b)r;");
    su::BPTree obs = tree.collapse();
    ASSERT(obs.get_structure() == exp.get_structure());
    ASSERT(obs.names == exp.names);
    ASSERT(vec_almost_equal(obs.lengths, exp.lengths));

    SUITE_END();
}

void test_unifrac_sample_counts() {
    SUITE_START("test unifrac sample counts");
    su::biom table("test.biom");
    const double* obs = table.get_sample_counts();
    const double exp[] = {7, 3, 4, 6, 3, 4};
    for(unsigned int i = 0; i < 6; i++)
        ASSERT(obs[i] == exp[i]);
    SUITE_END();
}

void test_set_tasks() {
    SUITE_START("test set tasks");
    std::vector<su::task_parameters> obs(1);
    std::vector<su::task_parameters> exp(1);

    exp[0].g_unifrac_alpha = 1.0;
    exp[0].n_samples = 100;
    exp[0].bypass_tips = false;
    exp[0].start = 0;
    exp[0].stop = 100;
    exp[0].tid = 0;

    set_tasks(obs, 1.0, 100, 0, 100, false, true, 1);
    ASSERT(obs[0].g_unifrac_alpha == exp[0].g_unifrac_alpha);
    ASSERT(obs[0].n_samples == exp[0].n_samples);
    ASSERT(obs[0].start == exp[0].start);
    ASSERT(obs[0].stop == exp[0].stop);
    ASSERT(obs[0].tid == exp[0].tid);

    std::vector<su::task_parameters> obs2(2);
    std::vector<su::task_parameters> exp2(2);

    exp2[0].g_unifrac_alpha = 1.0;
    exp2[0].n_samples = 100;
    exp2[0].bypass_tips = false;
    exp2[0].start = 0;
    exp2[0].stop = 50;
    exp2[0].tid = 0;
    exp2[1].g_unifrac_alpha = 1.0;
    exp2[1].n_samples = 100;
    exp2[1].bypass_tips = false;
    exp2[1].start = 50;
    exp2[1].stop = 100;
    exp2[1].tid = 1;

    set_tasks(obs2, 1.0, 100, 0, 100, false, true, 2);
    for(unsigned int i=0; i < 2; i++) {
        ASSERT(obs2[i].g_unifrac_alpha == exp2[i].g_unifrac_alpha);
        ASSERT(obs2[i].n_samples == exp2[i].n_samples);
        ASSERT(obs2[i].start == exp2[i].start);
        ASSERT(obs2[i].stop == exp2[i].stop);
        ASSERT(obs2[i].tid == exp2[i].tid);
    }

    std::vector<su::task_parameters> obs3(3);
    std::vector<su::task_parameters> exp3(3);

    exp3[0].g_unifrac_alpha = 1.0;
    exp3[0].n_samples = 100;
    exp3[0].bypass_tips = false;
    exp3[0].start = 25;
    exp3[0].stop = 50;
    exp3[0].tid = 0;
    exp3[1].g_unifrac_alpha = 1.0;
    exp3[1].n_samples = 100;
    exp3[1].bypass_tips = false;
    exp3[1].start = 50;
    exp3[1].stop = 75;
    exp3[1].tid = 1;
    exp3[2].g_unifrac_alpha = 1.0;
    exp3[2].n_samples = 100;
    exp3[2].bypass_tips = false;
    exp3[2].start = 75;
    exp3[2].stop = 100;
    exp3[2].tid = 2;

    set_tasks(obs3, 1.0, 100, 25, 100, false, true, 3);
    for(unsigned int i=0; i < 3; i++) {
        ASSERT(obs3[i].g_unifrac_alpha == exp3[i].g_unifrac_alpha);
        ASSERT(obs3[i].n_samples == exp3[i].n_samples);
        ASSERT(obs3[i].start == exp3[i].start);
        ASSERT(obs3[i].stop == exp3[i].stop);
        ASSERT(obs3[i].tid == exp3[i].tid);
    }

    std::vector<su::task_parameters> obs4(3);
    std::vector<su::task_parameters> exp4(3);

    exp4[0].g_unifrac_alpha = 1.0;
    exp4[0].n_samples = 100;
    exp4[0].bypass_tips = false;
    exp4[0].start = 26;
    exp4[0].stop = 51;
    exp4[0].tid = 0;
    exp4[1].g_unifrac_alpha = 1.0;
    exp4[1].n_samples = 100;
    exp4[1].bypass_tips = false;
    exp4[1].start = 51;
    exp4[1].stop = 76;
    exp4[1].tid = 1;
    exp4[2].g_unifrac_alpha = 1.0;
    exp4[2].n_samples = 100;
    exp4[2].bypass_tips = false;
    exp4[2].start = 76;
    exp4[2].stop = 100;
    exp4[2].tid = 2;

    set_tasks(obs4, 1.0, 100, 26, 100, false, true, 3);
    for(unsigned int i=0; i < 3; i++) {
        ASSERT(obs4[i].g_unifrac_alpha == exp4[i].g_unifrac_alpha);
        ASSERT(obs4[i].n_samples == exp4[i].n_samples);
        ASSERT(obs4[i].start == exp4[i].start);
        ASSERT(obs4[i].stop == exp4[i].stop);
        ASSERT(obs4[i].tid == exp4[i].tid);
    }

    // set_tasks boundary bug
    std::vector<su::task_parameters> obs16(16);
    std::vector<su::task_parameters> exp16(16);
    set_tasks(obs16, 1.0, 9511, 0, 0, false, true, 16);
    exp16[15].start = 4459;
    exp16[15].stop = 4756;
    ASSERT(obs16[15].start == exp16[15].start);
    ASSERT(obs16[15].stop == exp16[15].stop);
    SUITE_END();
}

void test_bptree_cstyle_constructor() {
    SUITE_START("test bptree constructor from c-style data");
                                //01234567
                                //11101000
    bool structure[] = {true, true, true, false, true, false, false, false};
    double lengths[] = {0, 0, 1, 0, 2, 0, 0, 0};
    const char* names[] = {"", "c", "123:foo; bar", "", "b", "", "", ""};
    su::BPTree tree(structure, lengths, names, 8);

    unsigned int exp_nparens = 8;

    std::vector<bool> exp_structure;
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(true);
    exp_structure.push_back(false);
    exp_structure.push_back(false);
    exp_structure.push_back(false);

    std::vector<uint32_t> exp_openclose;
    exp_openclose.push_back(7);
    exp_openclose.push_back(6);
    exp_openclose.push_back(3);
    exp_openclose.push_back(2);
    exp_openclose.push_back(5);
    exp_openclose.push_back(4);
    exp_openclose.push_back(1);
    exp_openclose.push_back(0);

    std::vector<std::string> exp_names;
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("c"));
    exp_names.push_back(std::string("123:foo; bar"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string("b"));
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());
    exp_names.push_back(std::string());

    std::vector<double> exp_lengths;
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(1.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(2.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);
    exp_lengths.push_back(0.0);

    ASSERT(tree.nparens == exp_nparens);
    ASSERT(tree.get_structure() == exp_structure);
    ASSERT(tree.get_openclose() == exp_openclose);
    ASSERT(tree.lengths == exp_lengths);
    ASSERT(tree.names == exp_names);

    SUITE_END();
}

void test_bptree_constructor_newline_bug() {
    SUITE_START("test bptree constructor newline bug");
    su::BPTree tree("((362be41f31fd26be95ae43a8769b91c0:0.116350803,(a16679d5a10caa9753f171977552d920:0.105836235,((a7acc2abb505c3ee177a12e514d3b994:0.008268754,(4e22aa3508b98813f52e1a12ffdb74ad:0.03144211,8139c4ac825dae48454fb4800fb87896:0.043622957)0.923:0.046588301)0.997:0.120902074,((2d3df7387323e2edcbbfcb6e56a02710:0.031543994,3f6752aabcc291b67a063fb6492fd107:0.091571442)0.759:0.016335166,((d599ebe277afb0dfd4ad3c2176afc50e:5e-09,84d0affc7243c7d6261f3a7d680b873f:0.010245188)0.883:0.048993011,51121722488d0c3da1388d1b117cd239:0.119447926)0.763:0.035660204)0.921:0.058191474)0.776:0.02854575)0.657:0.052060833)0.658:0.032547569,(99647b51f775c8ddde8ed36a7d60dbcd:0.173334268,(f18a9c8112372e2916a66a9778f3741b:0.194813398,(5833416522de0cca717a1abf720079ac:5e-09,(2bf1067d2cd4f09671e3ebe5500205ca:0.031692682,(b32621bcd86cb99e846d8f6fee7c9ab8:0.031330707,1016319c25196d73bdb3096d86a9df2f:5e-09)0.058:0.01028612)0.849:0.010284866)0.791:0.041353384)0.922:0.109470534):0.022169824000000005)root;\n\n");
    SUITE_END();
}
#endif

int main(int argc, char** argv) {
#ifndef API_ONLY
    test_bptree_constructor_simple();
    test_bptree_constructor_simple_cpp();
    test_bptree_constructor_newline_bug();
    test_bptree_constructor_from_existing();
    test_bptree_constructor_single_descendent();
    test_bptree_constructor_complex();
    test_bptree_constructor_semicolon();
    test_bptree_constructor_edgecases();
    test_bptree_constructor_quoted_comma();
    test_bptree_constructor_quoted_parens();
    test_bptree_cstyle_constructor();
    test_bptree_nullary();
    test_bptree_postorder();
    test_bptree_preorder();
    test_bptree_parent();
    test_bptree_leftchild();
    test_bptree_rightchild();
    test_bptree_rightsibling();
    test_bptree_get_tip_names();
    test_bptree_mask();
    test_bptree_shear_simple();
    test_bptree_shear_deep();
    test_bptree_collapse_simple();
    test_bptree_collapse_edge();
    test_bptree_constructor_inmem();

    test_biom_constructor();
    test_biom_constructor_from_sparse();
    test_biom_constructor_from_dense();
    test_biom_nullary();
    test_biom_get_obs_data();
    test_biom_filter();

    test_propstack_constructor();
    test_propstack_push_and_pop();
    test_propstack_get();

    test_unifrac_set_proportions();
    test_unifrac_set_proportions_range();
    test_unifrac_set_proportions_range_float();
    test_unifrac_deconvolute_stripes();
#endif

    test_unifrac_stripes_to_condensed_form_even();
    test_unifrac_stripes_to_condensed_form_odd();
    test_unifrac_stripes_to_condensed_form_odd2();

#ifndef API_ONLY
    test_unifrac_stripes_to_matrix_even();
    test_unifrac_stripes_to_matrix_odd();
    test_unifrac_stripes_to_matrix_odd2();
#endif

    test_unweighted_unifrac();
    test_unweighted_unifrac_fast();
    test_unnormalized_unweighted_unifrac();
    test_unnormalized_weighted_unifrac();
    test_normalized_weighted_unifrac();
    test_generalized_unifrac();
    test_vaw_unifrac_weighted_normalized();

#ifndef API_ONLY
    test_unifrac_sample_counts();
    test_set_tasks();
    test_test_table_ids_are_subset_of_tree();
#endif

    test_faith_pd();
    test_faith_pd_shear();
    test_faith_pd_inmem();
    test_subsample_inmem();
    test_permanova_inmem();
    test_matrix_inmem_substeps();
    test_partial_substeps();
    test_matrix_inmem_seeded();
    test_concurrent_matrix_inmem();
    test_concurrent_matrix_inmem_seeded();
    test_concurrent_faith_pd_inmem();
    test_concurrent_pcoa();
    test_concurrent_permanova_inmem();
    // must stay last; see the comment on the function
    test_concurrent_matrix_inmem_reporting();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);

    printf("\n THE END.\n");

    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
