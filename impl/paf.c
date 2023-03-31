#include "paf.h"
#include <ctype.h>
#include "../inc/paf.h"
#include "bioioC.h"

/*
 * Library functions for manipulating paf files.
 */

void paf_destruct(Paf *paf) {
    Cigar *c = paf->cigar;
    while(c != NULL) { // Cleanup the individual cigar records
        Cigar *c2 = c->next;
        free(c);
        c = c2;
    }
    free(paf);
}

static Cigar *parse_cigar_record(char **c) {
    Cigar *cigar = st_calloc(1, sizeof(Cigar)); // Allocate a cigar record
    // Calculate the number of characters representing the length of the record
    int64_t i=0;
    while(isdigit((*c)[i])) {
        i++;
    }
    char t = (*c)[i]; // The type of the cigar operation
    switch(t) {
    case 'M': // match
        cigar->op = match;
        break;
    case '=': // exact match
        cigar->op = sequence_match;
        break;
    case 'X': // snp match
        cigar->op = sequence_mismatch;
        break;
    case 'I':
        cigar->op = query_insert;
        break;
    case 'D':
        cigar->op = query_delete;
        break;
    default:
        st_errAbort("Got an unexpected character paf cigar string: %c\n", t);
        break;
    }
    (*c)[i] = ' ';
    cigar->length = atoll(*c);
    *c = &((*c)[i+1]);
    return cigar;
}

static Cigar *parse_cigar(char *cigar_string) {
    if(cigar_string[0] == '\0') { // If is the empty string
        return NULL;
    }
    Cigar *cigar = parse_cigar_record(&cigar_string);
    Cigar *pCigar = cigar;
    while(cigar_string[0] != '\0') {
        Cigar *nCigar = parse_cigar_record(&cigar_string);
        pCigar->next = nCigar;
        pCigar = nCigar;
    }
    return cigar;
}

static Cigar *cigar_reverse(Cigar *c) {
    if(c == NULL) {
        return NULL;
    }
    Cigar *p = NULL;
    // p -> c -> n -> q
    while(c->next != NULL) {
        Cigar *n = c->next;
        c->next = p; // p <- c , n -> q
        p = c;
        c = n;
    }
    c->next = p; // p <- c <- n
    return c;
}

Paf *paf_parse(char *paf_string) {
    Paf *paf = st_calloc(1, sizeof(Paf));

    stList *tokens = stString_split(paf_string); // Tokenize the record

    // Get query coordinates
    paf->query_name = stString_copy(stList_get(tokens, 0));
    paf->query_length = atoll(stList_get(tokens, 1));
    paf->query_start = atoll(stList_get(tokens, 2));
    paf->query_end = atoll(stList_get(tokens, 3));

    // Is the alignment forward or reverse
    char strand = ((char *)stList_get(tokens, 4))[0];
    if(strand != '+' && strand != '-') {
        st_errAbort("Got an unexpected strand character (%c) in a paf string: %s\n", strand, paf_string);
    }
    paf->same_strand = strand == '+';

    // Get target coordinates
    paf->target_name = stString_copy(stList_get(tokens, 5));
    paf->target_length = atoll(stList_get(tokens, 6));
    paf->target_start = atoll(stList_get(tokens, 7));
    paf->target_end = atoll(stList_get(tokens, 8));

    // Get the core alignment metric attributes of the record
    paf->num_matches = atoll(stList_get(tokens, 9));
    paf->num_bases = atoll(stList_get(tokens, 10));
    paf->mapping_quality = atoll(stList_get(tokens, 11));

    paf->tile_level = -1;
    paf->chain_id = -1;

    // Parse the remaining optional tags
    for(int64_t i=12; i<stList_length(tokens); i++) {
        stList *tag = stString_splitByString(stList_get(tokens, i), ":");
        char *type = stList_get(tag, 0);
        if(strcmp(type, "tp") == 0) {
            paf->type = ((char *)stList_get(tag, 2))[0];
            assert(paf->type == 'P' || paf->type == 'S' || paf->type == 'I');
        } else if (strcmp(type, "AS") == 0) {
            paf->score = atoll(stList_get(tag, 2));
        } else if(strcmp(type, "cg") == 0) {
            paf->cigar = parse_cigar(stList_get(tag, 2));
        }
        else if(strcmp(type, "tl") == 0) {
            paf->tile_level = atoll(stList_get(tag, 2));
        }
        else if(strcmp(type, "cn") == 0) {
            paf->chain_id = atoll(stList_get(tag, 2));
        }
        stList_destruct(tag);
    }

    // Cleanup
    stList_destruct(tokens);

    return paf;
}

Paf *paf_read(FILE *fh) {
    char *c = stFile_getLineFromFile(fh);
    if(c == NULL) {
        return NULL;
    }
    Paf *paf = paf_parse(c);
    free(c);

    paf_check(paf);

    return paf;
}

int64_t cigar_number_of_records(Paf *paf) {
    int64_t i=0;
    Cigar *c = paf->cigar;
    while(c != NULL) {
        i++;
        c = c->next;
    }
    return i;
}

char *paf_print(Paf *paf) {
    // Generous estimate of size needed for each paf record.
    int64_t buf_size = 12 * cigar_number_of_records(paf) + 140 + strlen(paf->query_name) + strlen(paf->target_name);
    char *buffer = st_malloc(sizeof(char) * buf_size); // Giving a generous
    int64_t i = sprintf(buffer, "%s\t%" PRIi64 "\t%" PRIi64"\t%" PRIi64"\t%c\t%s\t%" PRIi64"\t%" PRIi64"\t%" PRIi64
                                "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64,
                        paf->query_name, paf->query_length, paf->query_start, paf->query_end,
                        paf->same_strand ? '+' : '-',
                        paf->target_name, paf->target_length, paf->target_start, paf->target_end,
                        paf->num_matches, paf->num_bases, paf->mapping_quality);
    if(paf->type != '\0' || paf->tile_level != -1) {
        // if paf type not specified, use tile_level
        char t = paf->type;
        if (t == '\0') {
            t = paf->tile_level > 1 ? 'S' : 'P';
        }
        // sanity check (assumption: secondary alignment iff tile != 1)
        assert(paf->type != 'S' || paf->tile_level == -1 || paf->tile_level != 1);
        i += sprintf(buffer+i, "\ttp:A:%c", t);
    }
    if(paf->score != INT_MAX) {
        i += sprintf(buffer+i, "\tAS:i:%" PRIi64, paf->score);
    }
    if(paf->tile_level != -1) {
        i += sprintf(buffer+i, "\ttl:i:%" PRIi64, paf->tile_level);
    }
    if(paf->chain_id != -1) {
        i += sprintf(buffer+i, "\tcn:i:%" PRIi64, paf->chain_id);
    }
    if(i > buf_size) {
        st_errAbort("Size of paf record exceeded buffer size\n");
    }
    if(paf->cigar != NULL) {
        i += sprintf(buffer+i, "\tcg:Z:");
        Cigar *c = paf->cigar;
        while(c != NULL) {
            char op_char;
            switch(c->op) {
                case match:
                    op_char = 'M';
                    break;
                case query_insert:
                    op_char = 'I';
                    break;
                case query_delete:
                    op_char = 'D';
                    break;
                case sequence_match:
                    op_char = '=';
                    break;
                case sequence_mismatch:
                    op_char = 'X';
                    break;
            }
            i += sprintf(buffer+i, "%" PRIi64 "%c", c->length, op_char);
            c = c->next;
            if(i > buf_size) {
                st_errAbort("Size of paf record exceeded buffer size\n");
            }
        }
    }
    if(i > buf_size) {
        st_errAbort("Size of paf record exceeded buffer size\n");
    }
    return buffer;
}

void paf_stats_calc(Paf *paf, char *query_seq, char *target_seq,
                    int64_t *matches, int64_t *mismatches, int64_t *query_inserts, int64_t *query_deletes) {
    paf_encode_mismatches(paf, query_seq, target_seq);
    (*matches)=0; (*mismatches)=0; (*query_inserts)=0; (*query_deletes)=0;
    Cigar *c = paf->cigar;
    while(c != NULL) {
        if(c->op == sequence_match) {
            *matches += c->length;
        }
        else if(c->op == sequence_mismatch) {
            *mismatches += c->length;
        }
        else if(c->op == query_insert) {
            (*query_inserts)++;
        }
        else {
            assert(c->op == query_delete);
            (*query_deletes)++;
        }
        c = c->next;
    }
}

static void paf_pretty_print2(char *seq, int64_t i, int64_t j, FILE *fh) {
    char c = seq[j];
    seq[j] = '\0';
    fprintf(fh, "%s\n", &(seq[i]));
    seq[j] = c;
}

void paf_pretty_print(Paf *paf, char *query_seq, char *target_seq, FILE *fh, bool include_alignment) {
    int64_t matches, mismatches, query_inserts, query_deletes;
    paf_stats_calc(paf, query_seq, target_seq, &matches, &mismatches, &query_inserts, &query_deletes);
    fprintf(fh, "Query:%s\tQ-start:%" PRIi64 "\tQ-length:%" PRIi64 "\tTarget:%s\tT-start:%" PRIi64 "\tT-length:%" PRIi64 "\tSame-strand:%i\tScore:%" PRIi64 "\tIdentity:%f\tAligned-bases:%" PRIi64 "\tQuery-inserts:%" PRIi64 "\tQuery-deletes:%" PRIi64 "\n",
            paf->query_name, paf->query_start, paf->query_end-paf->query_start, paf->target_name, paf->target_start, paf->target_end-paf->target_start,
            (int)paf->same_strand, paf->score, (float)matches/(matches+mismatches), matches+mismatches, query_inserts, query_deletes);

    // Now print a base level alignment
    if(include_alignment) {
        Cigar *c = paf->cigar;
        int64_t max_align_length = paf->query_end - paf->query_start + paf->target_end - paf->target_start;
        char *query_align = st_malloc(sizeof(char) * (max_align_length + 1));
        char *target_align = st_malloc(sizeof(char) * (max_align_length + 1));
        char *star_align = st_malloc(sizeof(char) * (max_align_length + 1));
        int64_t i = 0, j = paf->target_start, k = 0;
        while (c != NULL) {
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
            c = c->next;
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

void paf_write(Paf *paf, FILE *fh) {
    char *c = paf_print(paf);
    fprintf(fh, "%s\n", c);
    free(c);
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
        Cigar *cigar = paf->cigar;
        while(cigar != NULL) {
            if(cigar->op != query_delete) {
                i += cigar->length;
            }
            if(cigar->op != query_insert) {
                j += cigar->length;
            }
            cigar = cigar->next;
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
    Cigar *c = paf->cigar;
    while(c != NULL) {
        if(c->op == query_insert) {
            c->op = query_delete;
        }
        else if(c->op == query_delete) {
            c->op = query_insert;
        }
        c = c->next;
    }
    // Now reverse the order if the ordering is swapped
    if(!paf->same_strand) {
        paf->cigar = cigar_reverse(paf->cigar);
    }
}

stList *read_pafs(FILE *fh) {
    stList *pafs = stList_construct3(0, (void (*)(void *))paf_destruct);
    Paf *paf;
    while((paf = paf_read(fh)) != NULL) {
        paf_check(paf);
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
    Cigar *c = paf->cigar;
    while(c != NULL) {
        if(c->op == match) {
            aligned_bases += c->length;
        }
        c = c->next;
    }
    return aligned_bases;
}

static Cigar *cigar_trim(int64_t *query_c, int64_t *target_c, Cigar *c, int64_t end_bases_to_trim, int q_sign, int t_sign) {
    int64_t bases_trimmed = 0;
    while(c != NULL && (c->op != match || bases_trimmed < end_bases_to_trim)) {
        if(c->op == match) { // can trim this alignment
            if(bases_trimmed + c->length > end_bases_to_trim) {
                int64_t i = end_bases_to_trim - bases_trimmed;
                c->length -= i;
                (*query_c) += q_sign*i;
                (*target_c) += t_sign*i;
                assert(c->length > 0);
                break;
            }
            bases_trimmed += c->length;
            (*query_c) += q_sign*c->length;
            (*target_c) += t_sign*c->length;
        }
        else if(c->op == query_insert) {
            (*query_c) += q_sign*c->length;
        }
        else {
            assert(c->op == query_delete);
            (*target_c) += t_sign*c->length;
        }
        Cigar *c2 = c;
        c = c->next;
        free(c2);
    }
    return c;
}

void paf_trim_ends(Paf *paf, int64_t end_bases_to_trim) {
    if(paf->same_strand) {
        // Trim front end
        paf->cigar = cigar_trim(&(paf->query_start), &(paf->target_start), paf->cigar, end_bases_to_trim, 1, 1);
        // Trim back end
        paf->cigar = cigar_reverse(cigar_trim(&(paf->query_end), &(paf->target_end), cigar_reverse(paf->cigar), end_bases_to_trim, -1, -1));
    }
    else {
        // Trim front end
        paf->cigar = cigar_trim(&(paf->query_end), &(paf->target_start), paf->cigar, end_bases_to_trim, -1, 1);
        // Trim back end
        paf->cigar = cigar_reverse(cigar_trim(&(paf->query_start), &(paf->target_end), cigar_reverse(paf->cigar), end_bases_to_trim, 1, -1));
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
    s_paf->cigar = st_calloc(1, sizeof(Cigar)); // Allocate a cigar record
    s_paf->cigar->length = length;
    s_paf->cigar->op = match;

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
    Cigar *p = paf->cigar;
    int64_t query_coordinate = paf->same_strand ? paf->query_start : paf->query_end;
    int64_t target_coordinate = paf->target_start;
    stList *matches = stList_construct3(0, (void (*)(void *))paf_destruct);
    while (p != NULL) {
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
        p = p->next;
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

SequenceCountArray *get_alignment_count_array(stHash *seq_names_to_alignment_count_arrays, Paf *paf) {
    SequenceCountArray *seq_count_array = stHash_search(seq_names_to_alignment_count_arrays, paf->query_name);
    if(seq_count_array == NULL) { // If the counts have not been initialized yet
        seq_count_array = st_calloc(1, sizeof(SequenceCountArray));
        seq_count_array->name = paf->query_name;
        seq_count_array->length = paf->query_length;
        seq_count_array->counts = st_calloc(paf->query_length, sizeof(uint16_t)); // sets all the counts to zero
        stHash_insert(seq_names_to_alignment_count_arrays, paf->query_name, seq_count_array); // adds to the hash
    }
    else {
        assert(seq_count_array->length == paf->query_length); // Check the name is unique
    }
    return seq_count_array;
}

void increase_alignment_level_counts(SequenceCountArray *seq_count_array, Paf *paf) {
    Cigar *c = paf->cigar;
    int64_t i = paf->query_start;
    while(c != NULL) {
        if(c->op != query_delete) {
            if(c->op == match) {
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
        c = c->next;
    }
    assert(i == paf->query_end);
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

static Cigar *paf_make_match_encoding_mismatches(int64_t target_offset, char *target_seq,
                                                 int64_t query_offset, char *query_seq,
                                                 int64_t length, bool same_strand,
                                                 Cigar **pCigar) {
    Cigar *cigar = NULL;
    // Calculate the length of the match
    int64_t match_run_length=0, mismatch_run_length=0;
    for(int64_t i=0; i<length; i++) {
        if(toupper(target_seq[target_offset + i]) == toupper(same_strand ?
                                                        query_seq[query_offset + i] :
                                                        stString_reverseComplementChar(query_seq[query_offset - i]))) {
            if(mismatch_run_length) {
                assert(match_run_length == 0);
                cigar = st_calloc(1, sizeof(Cigar)); // Allocate a cigar record
                cigar->op = sequence_mismatch;
                cigar->length = mismatch_run_length;
                *pCigar = cigar;
                pCigar = &(cigar->next);
                mismatch_run_length = 0;
            }
            match_run_length++;
        }
        else {
            if(match_run_length) {
                assert(mismatch_run_length == 0);
                cigar = st_calloc(1, sizeof(Cigar)); // Allocate a cigar record
                cigar->op = sequence_match;
                cigar->length = match_run_length;
                *pCigar = cigar;
                pCigar = &(cigar->next);
                match_run_length = 0;
            }
            mismatch_run_length++;
        }
    }
    cigar = st_calloc(1, sizeof(Cigar)); // Allocate a cigar record
    *pCigar = cigar; // Connect the previous cigar to it
    if(mismatch_run_length) {
        assert(match_run_length == 0);
        cigar->op = sequence_mismatch;
        cigar->length = mismatch_run_length;
    }
    if(match_run_length) {
        assert(mismatch_run_length == 0);
        cigar->op = sequence_match;
        cigar->length = match_run_length;
    }
    return cigar;
}

void paf_encode_mismatches(Paf *paf, char *query_seq, char *target_seq) {
    Cigar *c = paf->cigar, **pCigar = &(paf->cigar);
    int64_t i=0, j=paf->target_start;
    while(c != NULL) {
        if(c->op == match) {
            Cigar *c2 = paf_make_match_encoding_mismatches(j, target_seq,
                                                           paf->same_strand ? paf->query_start+i : paf->query_end-(i+1),
                                                           query_seq, c->length, paf->same_strand, pCigar);
            // Update the i and j coordinates
            i += c->length;
            j += c->length;
            // Connect up the new matches
            c2->next = c->next;
            // Cleanup the old match
            free(c);
            // Now set c to be the end of the new run of matches
            c = c2;
        }
        else if(c->op == query_insert) {
            i += c->length;
        }
        else if (c->op == query_delete) {
            j += c->length;
        }
        else { // Case the paf already has sequence matches and mismatches encoded
            assert(c->op == sequence_match || c->op == sequence_mismatch);
            i += c->length;
            j += c->length;
        }
        // Move to the next operation
        pCigar = &(c->next);
        c = c->next;
    }
}

void paf_remove_mismatches(Paf *paf) {
    Cigar *c = paf->cigar;
    while(c != NULL) {
        if(c->op == sequence_match || c->op == sequence_mismatch) {
            c->op = match; // relabel it a match
            while(c->next != NULL && (c->next->op == sequence_match || c->next->op == sequence_mismatch)) { // remove any
                // mismatches / matches
                Cigar *c2 = c->next;
                c->length += c2->length;
                c->next = c2->next;
                free(c2);
            }
        }
        c = c->next;
    }
}

Cigar *paf_trim_unreliable_ends2(Cigar *c, int64_t *matches, int64_t *mismatches, double identity_threshold,
                                 bool less_than, int64_t max_trim) {
    /*
     * Get the longest prefix with an identity less than the given identity threshold.
     * Also calculate the number of matches / mismatches in the alignment
     */
    (*matches)=0; (*mismatches)=0;
    Cigar *c2 = NULL;
    while(c != NULL) { // For each cigar op in sequence
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
            c2 = c;
        }
        c = c->next;
    }
    return c2; // Return longest prefix with identity < identity threshold
}

void paf_trim_upto(Paf *paf, Cigar *trim_upto) {
    /*
     * Remove the prefix of cigar operations up to but excluding the given op, trim_upto
     */
    while(paf->cigar != NULL && paf->cigar != trim_upto) {
        if (paf->cigar->op != query_insert) {
            paf->target_start += paf->cigar->length;
        }
        if (paf->cigar->op != query_delete) {
            if (paf->same_strand) {
                paf->query_start += paf->cigar->length;
            } else {
                paf->query_end -= paf->cigar->length;
            }
        }
        paf->cigar = paf->cigar->next;
    }
    assert(paf->cigar == trim_upto);
}

void paf_trim_unreliable_prefix(Paf *paf, float identity_threshold, int64_t max_trim) {
    /*
     * Trim a prefix of the paf with identity < identity_threshold. Will not trim more than max_trim of columns of the
     * alignment.
     */

    // Calculate largest prefix with identity less than the identity threshold
    int64_t matches, mismatches;
    Cigar *c = paf_trim_unreliable_ends2(paf->cigar, &matches, &mismatches, identity_threshold, 1, max_trim);

    if(c != NULL) { // If there is a prefix with low identity

        // Trim back the prefix to avoid chopping off a suffix of the prefix with identity
        // higher than the given threshold identity
        paf->cigar = cigar_reverse(paf->cigar); // Reverse the linked list of cigar ops
        Cigar *c2 = paf_trim_unreliable_ends2(c, &matches, &mismatches, identity_threshold, 0, -1); // Get the longest
                                              // suffix of the prefix with identity >= the identity threshold
        paf->cigar = cigar_reverse(paf->cigar); // Reverse back the linked list of cigar ops
        if(c2 != NULL) { // If c2 (the longest suffix) is not null, we trim up to c2 but not including it
            c = c2;
        }
        else { // otherwise, we want to trim everything including c, so we set c to c->next so
            // that we include c in the trim
            c = c->next;
        }

        // Now trim the prefix
        assert(c != NULL);
        paf_trim_upto(paf, c);
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

    paf_trim_unreliable_prefix(paf, identity_threshold, max_trim); // Trim the prefix
    paf_invert(paf); // Invert the paf
    paf_trim_unreliable_prefix(paf, identity_threshold, max_trim); // So trimming the suffix
    paf_invert(paf); // And invert it back

    // Debug output - check the final identity of the trimmed alignment is not less than the starting alignment
    int64_t trimmed_matches, trimmed_mismatches;
    paf_trim_unreliable_ends2(paf->cigar, &trimmed_matches, &trimmed_mismatches, 0, 1, -1);
    double final_identity = ((float)trimmed_matches)/(trimmed_matches + trimmed_mismatches); // Calculate the identity
    st_logDebug("Trimming unreliable prefix, got: %" PRIi64 " matches and %" PRIi64
    " mismatches, an alignment identity of %f and trim threshold of %f, after trimming got identity of %f with %" PRIi64
    " matches and %" PRIi64 " mismatches, using a max trim of %" PRIi64 " bases\n",
        matches, mismatches, identity, identity_threshold, final_identity, trimmed_matches, trimmed_mismatches, max_trim);
    assert(final_identity >= identity);
}
