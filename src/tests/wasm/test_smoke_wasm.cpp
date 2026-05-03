/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * WASM smoke test: prove libssu_wasm.a links cleanly against
 * libskbb_wasm.a under emcc and that pure C entry points round-trip
 * correctly under node.js.
 *
 * Exercises load_bptree_opaque + get_bptree_opaque_els +
 * destroy_bptree_opaque (pure C++ tree parsing, no skbb call), then
 * ssu_set_random_seed (skbb_set_random_seed via su::set_random_seed —
 * confirms libskbb_wasm.a is actually linked in).
 */

#include "tests/wasm/check_macros.hpp"
#include "api.hpp"

int main(void) {
    opaque_bptree_t* tree = nullptr;
    load_bptree_opaque("((a:1,b:1):1,c:1):0;", &tree);
    CHECK(tree != nullptr);

    int n_tips = get_bptree_opaque_els(tree);
    CHECK_EQ(n_tips, 3);

    destroy_bptree_opaque(&tree);
    CHECK(tree == nullptr);

    ssu_set_random_seed(42u);

    std::fprintf(stdout, "OK test_smoke_wasm\n");
    return 0;
}
