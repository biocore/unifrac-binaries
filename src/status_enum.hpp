#ifndef _UNIFRAC_STATUS_H
#define _UNIFRAC_STATUS_H

/* Append only -- the numeric values are part of the ABI. */
typedef enum compute_status {okay=0, tree_missing, table_missing, table_empty, unknown_method, table_and_tree_do_not_overlap, output_error, invalid_method, grouping_missing, unsupported_device} ComputeStatus;
typedef enum io_status {read_okay=0, write_okay, open_error, read_error, magic_incompatible, bad_header, unexpected_end, write_error} IOStatus;
typedef enum merge_status {merge_okay=0, incomplete_stripe_set, sample_id_consistency, square_mismatch, partials_mismatch, stripes_overlap} MergeStatus;

#endif
