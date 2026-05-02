/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 2 WASM smoke test.
 *
 * Goal: prove the unifrac-binaries WASM static archive (libssu_wasm.a)
 * links cleanly against scikit-bio-binaries' WASM static archive
 * (libskbb_wasm.a) under emcc, and that pure C entry points round-trip
 * correctly under node.js.
 *
 * Exercises:
 *   1. load_bptree_opaque — pure C++ tree parsing from a Newick string.
 *   2. get_bptree_opaque_els — confirms the parsed tree has the expected
 *      tip count.
 *   3. destroy_bptree_opaque — exercises the destructor path.
 *   4. ssu_set_random_seed — minimally exercises the skbb_set_random_seed
 *      hook through the unifrac wrapper. This validates that
 *      libskbb_wasm.a is actually linked in.
 *
 * Intentionally avoids any UniFrac math, PCoA, PERMANOVA, or Faith's PD
 * computation. Stage 3+ will cover those with native-generated ground
 * truth.
 */

#include "tests/wasm/check_macros.hpp"
#include "api.hpp"

int main(void) {
    // (1) Newick parse — three tips at distance 1 from a binary inner node.
    opaque_bptree_t* tree = nullptr;
    load_bptree_opaque("((a:1,b:1):1,c:1):0;", &tree);
    CHECK(tree != nullptr);

    // (2) Tip count matches the input. get_bptree_opaque_els returns the
    // number of tips, not n_parens.
    int n_tips = get_bptree_opaque_els(tree);
    CHECK_EQ(n_tips, 3);

    // (3) Destructor sets the pointer to NULL.
    destroy_bptree_opaque(&tree);
    CHECK(tree == nullptr);

    // (4) Random seed plumbing. No observable effect by itself; we just
    // need this to link without unresolved skbb_* symbols.
    ssu_set_random_seed(42u);
    ssu_set_random_seed(0u);

    std::fprintf(stdout, "OK test_smoke_wasm\n");
    return 0;
}
