#include "paf.h"
#include <ctype.h>
#include "bioioC.h"

/*
 * Library functions for manipulating paf files.
 */

// Fast int64 to string conversion, returns pointer to end of written string
static inline char *int64_to_str(char *buf, int64_t val) {
    if (val == 0) {
        *buf++ = '0';
        return buf;
    }
    char *start = buf;
    bool neg = val < 0;
    if (neg) {
        *buf++ = '-';
        val = -val;
        start = buf;
    }
    while (val > 0) {
        *buf++ = '0' + (val % 10);
        val /= 10;
    }
    // Reverse the digits
    char *end = buf - 1;
    while (start < end) {
        char tmp = *start;
        *start++ = *end;
        *end-- = tmp;
    }
    return buf;
}

// Fast string to int64 conversion - assumes valid input, no error checking
static inline int64_t str_to_int64(const char *s) {
    int64_t val = 0;
    bool neg = false;
    if (*s == '-') {
        neg = true;
        s++;
    }
    while (*s >= '0' && *s <= '9') {
        val = val * 10 + (*s++ - '0');
    }
    return neg ? -val : val;
}

void cigar_destruct(Cigar *c) {
    if (c != NULL) {
        free(c->recs);
        free(c);
    }
}

void paf_destruct(Paf *paf) {
    if(paf->cigar_string) { // Cleanup the string representing the cigar ops if we are storing it
        free(paf->cigar_string);
    }
    if(paf->cigar) { // cleanup the cigar as a linked list, if stored
        cigar_destruct(paf->cigar);
    }
    // Cleanup names
    free(paf->query_name);
    free(paf->target_name);
    free(paf);
}

Cigar *cigar_parse(char *cigar_string) {
    if(cigar_string[0] == '\0') { // If is the empty string
        return NULL;
    }
    // First pass: count operations
    int64_t count = 0;
    for (char *s = cigar_string; *s != '\0'; s++) {
        if (*s == 'M' || *s == 'I' || *s == 'D' || *s == '=' || *s == 'X') {
            count++;
        }
    }
    // Allocate container and array
    Cigar *cigar = st_malloc(sizeof(Cigar));
    cigar->recs = st_malloc(count * sizeof(CigarRecord));
    cigar->length = count;
    cigar->start = 0;
    cigar->capacity = count;
    // Second pass: fill records
    int64_t idx = 0;
    char *s = cigar_string;
    while (*s != '\0') {
        int64_t len = 0;
        while (*s >= '0' && *s <= '9') {
            len = len * 10 + (*s++ - '0');
        }
        CigarOp op;
        switch (*s) {
            case 'M': op = match; break;
            case '=': op = sequence_match; break;
            case 'X': op = sequence_mismatch; break;
            case 'I': op = query_insert; break;
            case 'D': op = query_delete; break;
            default: st_errAbort("Got an unexpected character paf cigar string: %c\n", *s); op = match; break;
        }
        cigar->recs[idx].length = len;
        cigar->recs[idx].op = op;
        idx++;
        s++;
    }
    assert(idx == count);
    return cigar;
}

static Cigar *cigar_construct_single(int64_t length, CigarOp op) {
    Cigar *c = st_malloc(sizeof(Cigar));
    c->recs = st_malloc(sizeof(CigarRecord));
    c->length = 1;
    c->start = 0;
    c->capacity = 1;
    c->recs[0].length = length;
    c->recs[0].op = op;
    return c;
}

static void cigar_reverse(Cigar *c) {
    if (c == NULL || c->length <= 1) return;
    int64_t lo = c->start;
    int64_t hi = c->start + c->length - 1;
    while (lo < hi) {
        CigarRecord tmp = c->recs[lo];
        c->recs[lo] = c->recs[hi];
        c->recs[hi] = tmp;
        lo++;
        hi--;
    }
}

Paf *paf_parse(char *paf_string, bool parse_cigar_string) {
    Paf *paf = st_calloc(1, sizeof(Paf));

    char *saveptr;
    char *token;

    // Field 0: query_name (must copy — stored in Paf)
    token = strtok_r(paf_string, "\t", &saveptr);
    paf->query_name = stString_copy(token);

    // Fields 1-3: query_length, query_start, query_end
    paf->query_length = str_to_int64(strtok_r(NULL, "\t", &saveptr));
    paf->query_start = str_to_int64(strtok_r(NULL, "\t", &saveptr));
    paf->query_end = str_to_int64(strtok_r(NULL, "\t", &saveptr));

    // Field 4: strand
    token = strtok_r(NULL, "\t", &saveptr);
    char strand = token[0];
    if(strand != '+' && strand != '-') {
        st_errAbort("Got an unexpected strand character (%c) in a paf string\n", strand);
    }
    paf->same_strand = strand == '+';

    // Field 5: target_name (must copy — stored in Paf)
    token = strtok_r(NULL, "\t", &saveptr);
    paf->target_name = stString_copy(token);

    // Fields 6-8: target_length, target_start, target_end
    paf->target_length = str_to_int64(strtok_r(NULL, "\t", &saveptr));
    paf->target_start = str_to_int64(strtok_r(NULL, "\t", &saveptr));
    paf->target_end = str_to_int64(strtok_r(NULL, "\t", &saveptr));

    // Fields 9-11: num_matches, num_bases, mapping_quality
    paf->num_matches = str_to_int64(strtok_r(NULL, "\t", &saveptr));
    paf->num_bases = str_to_int64(strtok_r(NULL, "\t", &saveptr));
    paf->mapping_quality = str_to_int64(strtok_r(NULL, "\t", &saveptr));

    // Set the following to default values to distinguish them from when they are initialized and 0
    paf->tile_level = -1;
    paf->chain_id = -1;
    paf->chain_score = -1;

    // Parse optional tags — format is always XX:T:value
    // Direct character indexing avoids stString_splitByString overhead
    while((token = strtok_r(NULL, "\t", &saveptr)) != NULL) {
        if(token[2] != ':' || token[4] != ':') {
            continue; // Skip malformed tags
        }
        char tag0 = token[0], tag1 = token[1];
        char *value = token + 5;

        if(tag0 == 't' && tag1 == 'p') {
            paf->type = value[0];
            assert(paf->type == 'P' || paf->type == 'S' || paf->type == 'I');
        } else if(tag0 == 'A' && tag1 == 'S') {
            paf->score = str_to_int64(value);
        } else if(tag0 == 'c' && tag1 == 'g') {
            if(parse_cigar_string) {
                paf->cigar = cigar_parse(value);
            } else {
                paf->cigar_string = stString_copy(value);
            }
        } else if(tag0 == 't' && tag1 == 'l') {
            paf->tile_level = str_to_int64(value);
        } else if(tag0 == 'c' && tag1 == 'n') {
            paf->chain_id = str_to_int64(value);
        } else if(tag0 == 's' && tag1 == '1') {
            paf->chain_score = str_to_int64(value);
        }
    }

    return paf;
}

Paf *paf_read_with_buffer(FILE *fh, bool parse_cigar_string, char **paf_buffer, int64_t *paf_length_buffer) {
    int64_t i = stFile_getLineFromFileWithBufferUnlocked(paf_buffer, paf_length_buffer, fh);
    if (i == -1 && strlen(*paf_buffer) == 0) {
        return NULL;
    }
    Paf *paf = paf_parse(*paf_buffer, parse_cigar_string);
    return paf;
}

Paf *paf_read(FILE *fh, bool parse_cigar_string) {
    int64_t buf_size = 4096;  // If too small, will get realloced
    char *buf = st_malloc(buf_size);
    Paf *paf = paf_read_with_buffer(fh, parse_cigar_string, &buf, &buf_size);
    free(buf);
    return paf;
}

Paf *paf_read2(FILE *fh) {
    return paf_read(fh, 1);
}

int64_t cigar_number_of_records(Paf *paf) {
    return cigar_count(paf->cigar);
}

void paf_stats_calc(Paf *paf, int64_t *matches, int64_t *mismatches, int64_t *query_inserts, int64_t *query_deletes,
                    int64_t *query_insert_bases, int64_t *query_delete_bases, bool zero_counts) {
    if(zero_counts) { // If requested to set these counts to zero
        (*matches) = 0; (*mismatches) = 0; (*query_inserts) = 0; (*query_deletes) = 0;
        (*query_insert_bases) = 0; (*query_delete_bases) = 0;
    }
    for (int64_t idx = 0; idx < cigar_count(paf->cigar); idx++) {
        CigarRecord *c = cigar_get(paf->cigar, idx);
        if(c->op == sequence_match || c->op == match) {
            *matches += c->length;
        }
        else if(c->op == sequence_mismatch) {
            *mismatches += c->length;
        }
        else if(c->op == query_insert) {
            (*query_inserts)++;
            (*query_insert_bases)+=c->length;
        }
        else {
            assert(c->op == query_delete);
            (*query_deletes)++;
            (*query_delete_bases)+=c->length;
        }
    }
}

static void paf_pretty_print2(char *seq, int64_t i, int64_t j, FILE *fh) {
    char c = seq[j];
    seq[j] = '\0';
    fprintf(fh, "%s\n", &(seq[i]));
    seq[j] = c;
}

void paf_pretty_print(Paf *paf, char *query_seq, char *target_seq, FILE *fh, bool include_alignment) {
    int64_t matches, mismatches, query_inserts, query_deletes, query_insert_bases, query_delete_bases;
    paf_stats_calc(paf, &matches, &mismatches, &query_inserts, &query_deletes, &query_insert_bases, &query_delete_bases, 1);
    fprintf(fh, "Query:%s\tQ-start:%" PRIi64 "\tQ-length:%" PRIi64 "\tTarget:%s\tT-start:%" PRIi64 "\tT-length:%"
            PRIi64 "\tSame-strand:%i\tScore:%" PRIi64 "\tIdentity:%f\tIdentity-with-gaps%f\tAligned-bases:%"
            PRIi64 "\tQuery-inserts:%" PRIi64 "\tQuery-deletes:%" PRIi64 "\n",
            paf->query_name, paf->query_start, paf->query_end-paf->query_start, paf->target_name, paf->target_start,
            paf->target_end-paf->target_start,
            (int)paf->same_strand, paf->score, (float)matches/(matches+mismatches),
            (float)matches/(matches+mismatches+query_insert_bases+query_delete_bases), matches+mismatches,
            query_inserts, query_deletes);

    // Now print a base level alignment
    if(include_alignment) {
        int64_t max_align_length = paf->query_end - paf->query_start + paf->target_end - paf->target_start;
        char *query_align = st_malloc(sizeof(char) * (max_align_length + 1));
        char *target_align = st_malloc(sizeof(char) * (max_align_length + 1));
        char *star_align = st_malloc(sizeof(char) * (max_align_length + 1));
        int64_t i = 0, j = paf->target_start, k = 0;
        for (int64_t ci = 0; ci < cigar_count(paf->cigar); ci++) {
            CigarRecord *c = cigar_get(paf->cigar, ci);
            for (int64_t l = 0; l < c->length; l++) {
                char m = '-', n = '-';
                if (c->op != query_insert) {
                    m = target_seq[j++];
                }
                if (c->op != query_delete) {
                    n = paf->same_strand ? query_seq[paf->query_start + i++] :
                        stString_reverseComplementChar(query_seq[paf->query_end - (++i)]);
                }
                target_align[k] = m;
                query_align[k] = n;
                star_align[k++] = toupper(m) == toupper(n) ? '*' : ' ';
            }
        }
        assert(k <= max_align_length);
        int64_t window = 150;
        for (int64_t l = 0; l < k; l += window) {
            paf_pretty_print2(target_align, l, l + window < k ? l + window : k, fh);
            paf_pretty_print2(query_align, l, l + window < k ? l + window : k, fh);
            paf_pretty_print2(star_align, l, l + window < k ? l + window : k, fh);
        }
        free(target_align);
        free(query_align);
        free(star_align);
    }
}

static int64_t paf_write_to_buffer(Paf *paf, char *paf_buffer) {
    char *p = paf_buffer;

    // Query name and coords
    int64_t len = strlen(paf->query_name);
    memcpy(p, paf->query_name, len); p += len;
    *p++ = '\t';
    p = int64_to_str(p, paf->query_length); *p++ = '\t';
    p = int64_to_str(p, paf->query_start); *p++ = '\t';
    p = int64_to_str(p, paf->query_end); *p++ = '\t';
    *p++ = paf->same_strand ? '+' : '-'; *p++ = '\t';

    // Target name and coords
    len = strlen(paf->target_name);
    memcpy(p, paf->target_name, len); p += len;
    *p++ = '\t';
    p = int64_to_str(p, paf->target_length); *p++ = '\t';
    p = int64_to_str(p, paf->target_start); *p++ = '\t';
    p = int64_to_str(p, paf->target_end); *p++ = '\t';

    // Core metrics
    p = int64_to_str(p, paf->num_matches); *p++ = '\t';
    p = int64_to_str(p, paf->num_bases); *p++ = '\t';
    p = int64_to_str(p, paf->mapping_quality);

    // Optional tags
    if(paf->type != '\0' || paf->tile_level != -1) {
        char t = paf->type;
        if (t == '\0') t = paf->tile_level > 1 ? 'S' : 'P';
        memcpy(p, "\ttp:A:", 6); p += 6;
        *p++ = t;
    }
    if(paf->score != INT_MAX) {
        memcpy(p, "\tAS:i:", 6); p += 6;
        p = int64_to_str(p, paf->score);
    }
    if(paf->tile_level != -1) {
        memcpy(p, "\ttl:i:", 6); p += 6;
        p = int64_to_str(p, paf->tile_level);
    }
    if(paf->chain_id != -1) {
        memcpy(p, "\tcn:i:", 6); p += 6;
        p = int64_to_str(p, paf->chain_id);
    }
    if(paf->chain_score != -1) {
        memcpy(p, "\ts1:i:", 6); p += 6;
        p = int64_to_str(p, paf->chain_score);
    }

    // CIGAR
    if(paf->cigar) {
        memcpy(p, "\tcg:Z:", 6); p += 6;
        for (int64_t ci = 0; ci < cigar_count(paf->cigar); ci++) {
            CigarRecord *c = cigar_get(paf->cigar, ci);
            p = int64_to_str(p, c->length);
            switch(c->op) {
                case match: *p++ = 'M'; break;
                case query_insert: *p++ = 'I'; break;
                case query_delete: *p++ = 'D'; break;
                case sequence_match: *p++ = '='; break;
                case sequence_mismatch: *p++ = 'X'; break;
                default: *p++ = 'N'; break;
            }
        }
    } else if(paf->cigar_string) {
        memcpy(p, "\tcg:Z:", 6); p += 6;
        len = strlen(paf->cigar_string);
        memcpy(p, paf->cigar_string, len); p += len;
    }

    *p++ = '\n';
    return p - paf_buffer;
}

static int64_t paf_estimate_buffer_size(Paf *paf) {
    // estimate size of buffer needed
    int64_t cigar_size = paf->cigar_string != NULL ? strlen(paf->cigar_string) :
                         12 * cigar_number_of_records(paf);
    return cigar_size + 300 + strlen(paf->query_name) + strlen(paf->target_name);
}

void paf_write_with_buffer(Paf *paf, FILE *fh, char **paf_buffer, int64_t *paf_length_buffer) {
    int64_t buf_size = paf_estimate_buffer_size(paf);
    if(buf_size > *paf_length_buffer) {
        *paf_length_buffer = buf_size * 2 + 1; // Double the buffer size to avoid constant reallocs
        *paf_buffer = realloc(*paf_buffer, *paf_length_buffer);
    }
    int64_t i = paf_write_to_buffer(paf, *paf_buffer);
    fwrite(*paf_buffer, 1, i, fh);
}

void paf_write(Paf *paf, FILE *fh) {
    int64_t buf_size = paf_estimate_buffer_size(paf);
    char stack_buf[4096];
    char *buf = buf_size <= 4096 ? stack_buf : st_malloc(buf_size); // Use stack buffer for small records, heap for large cigars
    int64_t i = paf_write_to_buffer(paf, buf);
    fwrite(buf, 1, i, fh);
    if(buf_size > 4096) {
        free(buf);
    }
}

char *paf_print(Paf *paf) {
    int64_t buf_size = paf_estimate_buffer_size(paf);
    char *buf = st_malloc(buf_size + 1);
    int64_t i = paf_write_to_buffer(paf, buf);
    buf[i-1] = '\0'; // Replace newline with a termination character
    return buf;
}

void paf_check(Paf *paf) {
    if(paf->query_start < 0 || paf->query_start >= paf->query_length) {
        st_errAbort("Paf query start coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->query_start > paf->query_end || paf->query_end > paf->query_length) {
        st_errAbort("Paf query end coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->target_start < 0 || paf->target_start >= paf->target_length) {
        st_errAbort("Paf target start coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->target_start > paf->target_end || paf->target_end > paf->target_length) {
        st_errAbort("Paf target end coordinates are invalid, %s", paf_print(paf));
    }
    if(paf->cigar != NULL) {
        // Check that cigar alignment, if present, matches the alignment:
        int64_t i=0, j=0;
        for (int64_t ci = 0; ci < cigar_count(paf->cigar); ci++) {
            CigarRecord *r = cigar_get(paf->cigar, ci);
            if(r->op != query_delete) {
                i += r->length;
            }
            if(r->op != query_insert) {
                j += r->length;
            }
        }
        if(i != paf->query_end - paf->query_start) {
            st_errAbort("Paf cigar alignment does not match query length: %" PRIi64 " vs. %" PRIi64 " %s", i,
                        paf->query_end - paf->query_start, paf_print(paf));
        }
        if(j != paf->target_end - paf->target_start) {
            st_errAbort("Paf cigar alignment does not match target length: %" PRIi64 " vs. %" PRIi64 " %s", j,
                        paf->target_end - paf->target_start, paf_print(paf));
        }
    }
}

static void swap(void **a, void **b) {
    void *c = *a;
    *a = *b;
    *b = c;
}

void paf_invert(Paf *paf) {
    // Swap the query and target coordinates
    swap((void **)&paf->query_start, (void **)&paf->target_start);
    swap((void **)&paf->query_end, (void **)&paf->target_end);
    swap((void **)&paf->query_length, (void **)&paf->target_length);
    swap((void **)&paf->query_name, (void **)&paf->target_name);

    // Switch the query and target in the cigar
    for (int64_t ci = 0; ci < cigar_count(paf->cigar); ci++) {
        CigarRecord *c = cigar_get(paf->cigar, ci);
        if(c->op == query_insert) {
            c->op = query_delete;
        }
        else if(c->op == query_delete) {
            c->op = query_insert;
        }
    }
    // Now reverse the order if the ordering is swapped
    if(!paf->same_strand) {
        cigar_reverse(paf->cigar);
    }
}

stList *read_pafs(FILE *fh, bool parse_cigar_string) {
    stList *pafs = stList_construct3(0, (void (*)(void *))paf_destruct);
    Paf *paf;
    while((paf = paf_read(fh, parse_cigar_string)) != NULL) {
        stList_append(pafs, paf);
    }
    return pafs;
}

void write_pafs(FILE *fh, stList *pafs) {
    for(int64_t i=0; i<stList_length(pafs); i++) {
        paf_write(stList_get(pafs, i), fh);
    }
}

int64_t paf_get_number_of_aligned_bases(Paf *paf) {
    int64_t aligned_bases = 0;
    for (int64_t ci = 0; ci < cigar_count(paf->cigar); ci++) {
        CigarRecord *r = cigar_get(paf->cigar, ci);
        if(r->op == match || r->op == sequence_match || r->op == sequence_mismatch) {
            aligned_bases += r->length;
        }
    }
    return aligned_bases;
}

static void cigar_trim(int64_t *query_c, int64_t *target_c, Cigar *c, int64_t end_bases_to_trim, int q_sign, int t_sign) {
    int64_t bases_trimmed = 0;
    while(c->length > 0 && ((cigar_get(c, 0)->op != match && cigar_get(c, 0)->op != sequence_match && cigar_get(c, 0)->op != sequence_mismatch) || bases_trimmed < end_bases_to_trim)) {
        CigarRecord *r = cigar_get(c, 0);
        if(r->op == match || r->op == sequence_match || r->op == sequence_mismatch) { // can trim this alignment
            if(bases_trimmed + r->length > end_bases_to_trim) {
                int64_t i = end_bases_to_trim - bases_trimmed;
                r->length -= i;
                (*query_c) += q_sign*i;
                (*target_c) += t_sign*i;
                assert(r->length > 0);
                break;
            }
            bases_trimmed += r->length;
            (*query_c) += q_sign*r->length;
            (*target_c) += t_sign*r->length;
        }
        else if(r->op == query_insert) {
            (*query_c) += q_sign*r->length;
        }
        else {
            assert(r->op == query_delete);
            (*target_c) += t_sign*r->length;
        }
        c->start++;
        c->length--;
    }
}

void paf_trim_ends(Paf *paf, int64_t end_bases_to_trim) {
    if(paf->same_strand) {
        // Trim front end
        cigar_trim(&(paf->query_start), &(paf->target_start), paf->cigar, end_bases_to_trim, 1, 1);
        // Trim back end
        cigar_reverse(paf->cigar);
        cigar_trim(&(paf->query_end), &(paf->target_end), paf->cigar, end_bases_to_trim, -1, -1);
        cigar_reverse(paf->cigar);
    }
    else {
        // Trim front end
        cigar_trim(&(paf->query_end), &(paf->target_start), paf->cigar, end_bases_to_trim, -1, 1);
        // Trim back end
        cigar_reverse(paf->cigar);
        cigar_trim(&(paf->query_start), &(paf->target_end), paf->cigar, end_bases_to_trim, 1, -1);
        cigar_reverse(paf->cigar);
    }
}

void paf_trim_end_fraction(Paf *p, float percentage) {
    // trim
    assert(percentage >= 0 && percentage <= 1.0);
    int64_t aligned_bases = paf_get_number_of_aligned_bases(p);
    int64_t end_bases_to_trim = (aligned_bases * percentage)/2.0;
    st_logDebug("For alignment of %" PRIi64 " query bases, %"
    PRIi64 " target bases and %" PRIi64 " aligned bases trimming %" PRIi64 " bases from each paf end\n",
            p->query_end - p->query_start, p->target_end - p->target_start, aligned_bases, end_bases_to_trim);
    paf_trim_ends(p, end_bases_to_trim);
}

Paf *paf_shatter2(Paf *paf, int64_t query_start, int64_t target_start, int64_t length) {
    Paf *s_paf = st_calloc(1, sizeof(Paf));

    s_paf->query_name = stString_copy(paf->query_name);
    s_paf->query_length = paf->query_length;
    s_paf->query_start = query_start;
    s_paf->query_end = query_start + length;

    s_paf->target_name = stString_copy(paf->target_name);
    s_paf->target_length = paf->target_length;
    s_paf->target_start = target_start;
    s_paf->target_end = target_start + length;

    s_paf->same_strand = paf->same_strand;
    s_paf->cigar = cigar_construct_single(length, match);

    s_paf->score = paf->score;
    s_paf->mapping_quality = paf->mapping_quality;
    s_paf->num_matches = length;
    s_paf->num_bases = length;
    s_paf->tile_level = paf->tile_level;
    s_paf->type = paf->type;
    s_paf->chain_id = paf->chain_id;

    paf_check(s_paf);

    return s_paf;
}

stList *paf_shatter(Paf *paf) {
    int64_t query_coordinate = paf->same_strand ? paf->query_start : paf->query_end;
    int64_t target_coordinate = paf->target_start;
    stList *matches = stList_construct3(0, (void (*)(void *))paf_destruct);
    for (int64_t ci = 0; ci < cigar_count(paf->cigar); ci++) {
        CigarRecord *p = cigar_get(paf->cigar, ci);
        assert(p->length >= 1);
        if (p->op == match) {
            if (paf->same_strand) {
                stList_append(matches, paf_shatter2(paf, query_coordinate, target_coordinate, p->length));
                query_coordinate += p->length;
            } else {
                query_coordinate -= p->length;
                stList_append(matches, paf_shatter2(paf, query_coordinate, target_coordinate, p->length));
            }
            target_coordinate += p->length;
        }
        else if (p->op == query_insert) {
            query_coordinate += paf->same_strand ? p->length : -p->length;
        }
        else {
            assert(p->op == query_delete);
            target_coordinate += p->length;
        }
    }
    assert(target_coordinate == paf->target_end);
    if (paf->same_strand) {
        assert(query_coordinate == paf->query_end);
    }
    else {
        assert(query_coordinate == paf->query_start);
    }

    return matches;
}

/*
 * Functions used by paf_tile and paf_to_bed
 */

void sequenceCountArray_destruct(SequenceCountArray *seq_count_array) {
    free(seq_count_array->name);
    free(seq_count_array->counts);
    free(seq_count_array);
}

SequenceCountArray *get_alignment_count_array(stHash *seq_names_to_alignment_count_arrays, Paf *paf) {
    SequenceCountArray *seq_count_array = stHash_search(seq_names_to_alignment_count_arrays, paf->query_name);
    if(seq_count_array == NULL) { // If the counts have not been initialized yet
        seq_count_array = st_calloc(1, sizeof(SequenceCountArray));
        seq_count_array->name = stString_copy(paf->query_name);
        seq_count_array->length = paf->query_length;
        seq_count_array->counts = st_calloc(paf->query_length, sizeof(uint16_t)); // sets all the counts to zero
        stHash_insert(seq_names_to_alignment_count_arrays, seq_count_array->name, seq_count_array); // adds to the hash
    }
    else {
        assert(seq_count_array->length == paf->query_length); // Check the name is unique
    }
    return seq_count_array;
}

void increase_alignment_level_counts(SequenceCountArray *seq_count_array, Paf *paf) {
    int64_t i = paf->query_start;
    for (int64_t ci = 0; ci < cigar_count(paf->cigar); ci++) {
        CigarRecord *c = cigar_get(paf->cigar, ci);
        if(c->op != query_delete) {
            if(c->op != query_insert) { // Is some kind of match
                assert(c->op == match || c->op == sequence_match || c->op == sequence_mismatch);
                for(int64_t j=0; j<c->length; j++) {
                    assert(i + j < paf->query_end && i + j >= 0 && i + j < paf->query_length);
                    assert(i + j < seq_count_array->length);
                    if(seq_count_array->counts[i + j] < INT16_MAX - 1) { // prevent overflow
                        seq_count_array->counts[i + j]++;
                    }
                }
            }
            i += c->length;
        }
    }
    assert(i == paf->query_end);
}

void interval_destruct(Interval *interval) {
    free(interval->name);
    free(interval);
}

Interval *decode_fasta_header(char *fasta_header) {
    Interval *i = st_calloc(1, sizeof(Interval));
    stList *attributes = fastaDecodeHeader(fasta_header);
    //Decode attributes
    int64_t j = sscanf((const char *) stList_peek(attributes), "%" PRIi64 "", &i->start);
    (void) j;
    free(stList_pop(attributes));
    assert(j == 1);
    j = sscanf((const char *) stList_peek(attributes), "%" PRIi64 "", &i->length);
    free(stList_pop(attributes));
    assert(j == 1);
    //Now relabel attributes
    i->name = fastaEncodeHeader(attributes);
    stList_destruct(attributes);
    return i;
}

int cmp_intervals(const void *i, const void *j) {
    Interval *x = (Interval *)i, *y = (Interval *)j;
    int k = strcmp(x->name, y->name);
    return k == 0 ? (x->start < y->start ? -1 : (x->start > y->start ? 1 : 0)) : k;
}

static int64_t count_mismatch_records(int64_t target_offset, char *target_seq,
                                       int64_t query_offset, char *query_seq,
                                       int64_t length, bool same_strand) {
    int64_t count = 0;
    bool prev_match = false, first = true;
    for(int64_t i=0; i<length; i++) {
        bool is_match = toupper(target_seq[target_offset + i]) == toupper(same_strand ?
                            query_seq[query_offset + i] :
                            stString_reverseComplementChar(query_seq[query_offset - i]));
        if(first || is_match != prev_match) {
            count++;
            first = false;
        }
        prev_match = is_match;
    }
    return count;
}

static int64_t fill_mismatch_records(int64_t target_offset, char *target_seq,
                                      int64_t query_offset, char *query_seq,
                                      int64_t length, bool same_strand,
                                      CigarRecord *dest) {
    int64_t idx = -1;
    bool prev_match = false, first = true;
    for(int64_t i=0; i<length; i++) {
        bool is_match = toupper(target_seq[target_offset + i]) == toupper(same_strand ?
                            query_seq[query_offset + i] :
                            stString_reverseComplementChar(query_seq[query_offset - i]));
        if(first || is_match != prev_match) {
            idx++;
            dest[idx].op = is_match ? sequence_match : sequence_mismatch;
            dest[idx].length = 1;
            first = false;
        } else {
            dest[idx].length++;
        }
        prev_match = is_match;
    }
    return idx + 1;
}

void paf_encode_mismatches(Paf *paf, char *query_seq, char *target_seq) {
    Cigar *cigar = paf->cigar;
    if(cigar == NULL) return;

    // First pass: count output records
    int64_t total = 0;
    int64_t qi = 0, tj = paf->target_start;
    for(int64_t idx = 0; idx < cigar->length; idx++) {
        CigarRecord *r = cigar_get(cigar, idx);
        if(r->op == match) {
            total += count_mismatch_records(tj, target_seq,
                        paf->same_strand ? paf->query_start+qi : paf->query_end-(qi+1),
                        query_seq, r->length, paf->same_strand);
            qi += r->length;
            tj += r->length;
        } else {
            if(r->op == query_insert) {
                qi += r->length;
            } else if(r->op == query_delete) {
                tj += r->length;
            } else {
                assert(r->op == sequence_match || r->op == sequence_mismatch);
                qi += r->length;
                tj += r->length;
            }
            total++;
        }
    }

    // Second pass: allocate and fill
    CigarRecord *new_recs = st_malloc(total * sizeof(CigarRecord));
    int64_t out = 0;
    qi = 0; tj = paf->target_start;
    for(int64_t idx = 0; idx < cigar->length; idx++) {
        CigarRecord *r = cigar_get(cigar, idx);
        if(r->op == match) {
            out += fill_mismatch_records(tj, target_seq,
                        paf->same_strand ? paf->query_start+qi : paf->query_end-(qi+1),
                        query_seq, r->length, paf->same_strand, &new_recs[out]);
            qi += r->length;
            tj += r->length;
        } else {
            new_recs[out++] = *r;
            if(r->op == query_insert) {
                qi += r->length;
            } else if(r->op == query_delete) {
                tj += r->length;
            } else {
                qi += r->length;
                tj += r->length;
            }
        }
    }
    assert(out == total);

    // Swap in new array
    free(cigar->recs);
    cigar->recs = new_recs;
    cigar->length = total;
    cigar->start = 0;
    cigar->capacity = total;
}

void paf_remove_mismatches(Paf *paf) {
    Cigar *cigar = paf->cigar;
    if(cigar == NULL) return;
    int64_t write = 0;
    for(int64_t read = 0; read < cigar->length; read++) {
        CigarRecord *r = cigar_get(cigar, read);
        if(r->op == sequence_match || r->op == sequence_mismatch || r->op == match) {
            if(write > 0 && cigar_get(cigar, write - 1)->op == match) {
                cigar_get(cigar, write - 1)->length += r->length;
            } else {
                CigarRecord *w = cigar_get(cigar, write);
                w->length = r->length;
                w->op = match;
                write++;
            }
        } else {
            if(write != read) {
                *cigar_get(cigar, write) = *r;
            }
            write++;
        }
    }
    cigar->length = write;
}

static int64_t paf_trim_unreliable_ends2(Cigar *cigar, int64_t *matches, int64_t *mismatches, double identity_threshold,
                                 bool less_than, int64_t max_trim) {
    /*
     * Get the longest prefix with an identity less than the given identity threshold.
     * Also calculate the number of matches / mismatches in the alignment.
     * Returns the index of the last element in that prefix, or -1 if no such prefix exists.
     */
    (*matches)=0; (*mismatches)=0;
    int64_t trim_idx = -1;
    for(int64_t idx = 0; idx < cigar_count(cigar); idx++) {
        CigarRecord *c = cigar_get(cigar, idx);
        // First add the op to the number of matches / mismatches
        if(c->op == sequence_match || c->op == match) {
            (*matches) += c->length;
        }
        else { // Note indels are treated as mismatches
            (*mismatches) += c->length;
        }
        if(max_trim >= 0 && *matches + *mismatches > max_trim) { // Don't trim if longer a threshold number of bases
            break;
        }
        double prefix_identity = ((float)(*matches))/((*matches) + (*mismatches));
        if((less_than && prefix_identity < identity_threshold) ||
           (!less_than && prefix_identity >= identity_threshold)) { // Including this op, does the prefix have identity
            // < identity threshold (or >= if less_than is false)
            trim_idx = idx;
        }
    }
    return trim_idx; // Return index of longest prefix with identity < identity threshold
}

static void paf_trim_upto(Paf *paf, int64_t trim_count) {
    /*
     * Remove the first trim_count cigar records from the front, adjusting coordinates
     */
    for(int64_t i = 0; i < trim_count; i++) {
        CigarRecord *r = cigar_get(paf->cigar, i);
        if (r->op != query_insert) {
            paf->target_start += r->length;
        }
        if (r->op != query_delete) {
            if (paf->same_strand) {
                paf->query_start += r->length;
            } else {
                paf->query_end -= r->length;
            }
        }
    }
    paf->cigar->start += trim_count;
    paf->cigar->length -= trim_count;
}

void paf_trim_unreliable_prefix(Paf *paf, float identity_threshold, float identity, int64_t max_trim) {
    /*
     * Trim a prefix of the paf with identity < identity_threshold. Will not trim more than max_trim of columns of the
     * alignment.
     */

    // Calculate largest prefix with identity less than the identity threshold
    int64_t matches, mismatches;
    int64_t trim_idx = paf_trim_unreliable_ends2(paf->cigar, &matches, &mismatches, identity_threshold, 1, max_trim);

    if(trim_idx >= 0) { // If there is a prefix with low identity

        // Find the longest suffix of the prefix [0..trim_idx] with identity >= alignment identity
        // by iterating backward from trim_idx
        int64_t suffix_matches = 0, suffix_mismatches = 0;
        int64_t best_suffix_start = -1;
        for(int64_t i = trim_idx; i >= 0; i--) {
            CigarRecord *r = cigar_get(paf->cigar, i);
            if(r->op == sequence_match || r->op == match) {
                suffix_matches += r->length;
            } else {
                suffix_mismatches += r->length;
            }
            double suffix_identity = ((float)suffix_matches) / (suffix_matches + suffix_mismatches);
            if(suffix_identity >= identity) {
                best_suffix_start = i;
            }
        }

        int64_t trim_count;
        if(best_suffix_start >= 0) {
            trim_count = best_suffix_start; // trim up to but not including best_suffix_start
        } else {
            trim_count = trim_idx + 1; // trim everything including trim_idx
        }

        // Now trim the prefix
        if(trim_count > 0) {
            paf_trim_upto(paf, trim_count);
        }
    }
}

void paf_trim_unreliable_tails(Paf *paf, float score_fraction, float max_fraction_to_trim) {
    /*
     * Algorithm, adapted from proposal by Bob Harris:
     *
     * First compute the average identity, i, over the whole alignment (m/(m+mm+i+d),
     * (b) identify the longest prefix for which the identity is significantly less than i
     * (i - i * score_fraction), then (c) shorten the prefix by re-including the longest suffix of the prefix that has
     * an identity equal or higher to i, (d) trim this prefix, (e) repeat for the suffix, still using the original i.
     * This is not quite symmetric, but it works okay in practice.
     *
     * The tails trimmed will not exceed more than max_fraction_to_trim of the original alignment.
     */

    // Get global match/mismatch stats
    int64_t matches, mismatches;
    paf_trim_unreliable_ends2(paf->cigar, &matches, &mismatches, 0, 1, -1);

    double identity = ((float)matches)/(matches + mismatches); // Calculate the identity
    double identity_threshold = identity - (identity * score_fraction); // Calculate the identity threshold below which to
    // trim the tail
    int64_t max_trim = (matches + mismatches) * max_fraction_to_trim;

    paf_trim_unreliable_prefix(paf, identity_threshold, identity, max_trim); // Trim the prefix
    paf_invert(paf); // Invert the paf
    paf_trim_unreliable_prefix(paf, identity_threshold, identity, max_trim); // So trimming the suffix
    paf_invert(paf); // And invert it back

    // Debug output - check the final identity of the trimmed alignment is not less than the starting alignment
    int64_t trimmed_matches, trimmed_mismatches;
    paf_trim_unreliable_ends2(paf->cigar, &trimmed_matches, &trimmed_mismatches, 0, 1, -1);
    double final_identity =
            ((float) trimmed_matches) / (trimmed_matches + trimmed_mismatches); // Calculate the identity
    if((trimmed_matches != matches || trimmed_mismatches != mismatches) && final_identity > identity + 0.1) {
        st_logDebug("Trimming unreliable prefix, got: %"
        PRIi64
        " matches and %"
        PRIi64
        " mismatches, an alignment identity of %f and trim threshold of %f, after trimming got identity of %f with %"
        PRIi64
        " matches and %"
        PRIi64
        " mismatches, using a max trim of %"
        PRIi64
        " bases\n",
                matches, mismatches, identity, identity_threshold, final_identity, trimmed_matches, trimmed_mismatches, max_trim);
    }
    assert(final_identity >= identity);
}
