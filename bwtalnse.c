#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <pthread.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bwa.h"
#include "kthread.h"
#include "kstring.h"

//From bwtaln.c
bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);
//From bwase.c
void bwase_initialize();
void bwa_cal_pac_pos_core(const bntseq_t *bns,
                          const bwt_t *bwt,
                          bwa_seq_t *seq,
                          const int max_mm,
                          const float fnr);
bwtint_t bwa_sa2pos(const bntseq_t *bns,
                    const bwt_t *bwt,
                    bwtint_t sapos,
                    int ref_len,
                    int *strand);
bwa_cigar_t *bwa_refine_gapped_core(bwtint_t l_pac,
                                    const ubyte_t *pacseq,
                                    int len,
                                    ubyte_t *seq,
                                    int ref_shift,
                                    bwtint_t *_rb,
                                    int *n_cigar);
char *bwa_cal_md1(int n_cigar,
                  bwa_cigar_t *cigar,
                  int len,
                  bwtint_t pos,
                  ubyte_t *seq,
                  bwtint_t l_pac,
                  const ubyte_t *pacseq,
                  kstring_t *str,
                  int *_nm);
void bwa_correct_trimmed(bwa_seq_t *s);
void bwa_print_sam1(const bntseq_t *bns,
                    bwa_seq_t *p,
                    const bwa_seq_t *mate,
                    int mode,
                    int max_top2);

void bwa_aln2seq_alnse(int n_aln,
                       const bwt_aln1_t *aln,
                       bwa_seq_t *s,
                       int set_main,
                       uint32_t n_multi)
{
    int i, cnt;
    if (n_aln == 0) {
        s->type = BWA_TYPE_NO_MATCH;
        s->c1 = s->c2 = 0;
        return;
    }
    /*
     The first best scoring alignment will be reported as main alignment.
     Rest n_multi of alignments will be reported in the XA tag.
    */
    if (set_main) {
        const bwt_aln1_t *maln = aln;
        s->n_mm = maln->n_mm; s->n_gapo = maln->n_gapo; s->n_gape = maln->n_gape;
        s->ref_shift = (int)maln->n_del - (int)maln->n_ins;
        s->score = maln->score;
        //s->sa = maln->k + (bwtint_t)((maln->l - maln->k + 1) * drand48());
        //First occurrence is main alignment
        s->sa = maln->k;
        //Count number of occurrences of all equally best scoring alignments
        for (i = cnt = 0; i < n_aln; ++i) {
            const bwt_aln1_t *p = aln + i;
            if (p->score > maln->score) break;
            cnt += p->l - p->k + 1;
        }
        s->c1 = cnt;
        for (; i < n_aln; ++i) cnt += aln[i].l - aln[i].k + 1;
        s->c2 = cnt - s->c1;
        s->type = s->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
    }
    if (n_multi) {
        int k, rest, n_occ, z = 0;
        for (k = 0, n_occ = 0; k < n_aln; ++k) {
            const bwt_aln1_t *q = aln + k;
            n_occ += q->l - q->k + 1;
        }
        n_occ--; //Subtract main alignment occurrence
        if (s->multi) free(s->multi);
        /*
          Take the next rest alignments and report them in the XA tag
          TODO understand why bwt_aln1_t->l - bwt_aln1_t->k + 1 computes
               the number of occurrences of an alignment
               I think it's the SA interval
        */
        rest = n_occ > n_multi? n_multi : n_occ;
        s->multi = calloc(rest, sizeof(bwt_multi1_t));
        s->n_multi = rest;
        for (k = 0; k < n_aln; ++k) {
            const bwt_aln1_t *q = aln + k;
            bwtint_t l;
            /*
             Report alignments starting from the second occurrence of the
             first alignment.
            */
            for (l = k?q->k:q->k+1; l <= q->l; ++l) {
                s->multi[z].pos = l;
                s->multi[z].gap = q->n_gapo + q->n_gape;
                s->multi[z].ref_shift = (int)q->n_del - (int)q->n_ins;
                s->multi[z++].mm = q->n_mm;
                if (rest == z) goto velociraptor;
            }
        }
        //https://xkcd.com/292/
        velociraptor:
    }
}

void bwa_cal_sa_reg_gap1(bwt_t *const bwt,
                         bwa_seq_t *p,
                         const gap_opt_t *opt,
                         uint32_t n_occ,
                         gap_stack_t *stack,
                         bwt_width_t *seed_w,
                         bwt_width_t *w)
{
    int j;
    gap_opt_t local_opt = *opt;
    if (opt->fnr > 0.0)
        local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
    //maximum value of gap opens can not be greater that the maximum number of mismatches
    if (local_opt.max_diff < local_opt.max_gapo)
        local_opt.max_gapo = local_opt.max_diff;
    p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;

    memset(w, 0, (local_opt.max_len + 1) * sizeof(bwt_width_t));
    //Set w to something TODO underdtand
    bwt_cal_width(bwt, p->len, p->seq, w);
    //Disable seeding if seed len > read len
    local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
    //Do same as with w but just for the first seed len bases
    if (p->len > opt->seed_len) {
        memset(seed_w, 0, (opt->seed_len + 1) * sizeof(bwt_width_t));
        bwt_cal_width(bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), seed_w);
    }
    /*
    *  we need to complement will be complemented back in bwa_refine_gapped1
    */
    for (j = 0; j < p->len; ++j)
        p->seq[j] = p->seq[j] > 3? 4 : 3 - p->seq[j];
    // core function
    p->aln = bwt_match_gap(bwt,
                           p->len,
                           p->seq,
                           w,
                           p->len <= opt->seed_len? 0 : seed_w,
                           &local_opt,
                           &p->n_aln,
                           stack);
    //Set main hit and up to n_occ other hits
    bwa_aln2seq_alnse(p->n_aln, p->aln, p, 1, n_occ);
}

/*
Get main hit position, strand, and mapping quality.
Position and strand are also computed for other hits.
Here is where bwa_seq_t->type = BWA_TYPE_NO_MATCH if read unmapped.
*/
void bwa_cal_pac_pos1(const bntseq_t *bns,
                      bwt_t *const bwt,
                      bwa_seq_t *p, //Main hit
                      int max_mm,   //Max mismatch
                      float fnr)    //Max difference given prob and mm rate
{
    if (p->type == BWA_TYPE_NO_MATCH)
        return;
    int j, strand, n_multi;
    bwa_cal_pac_pos_core(bns, bwt, p, max_mm, fnr);
    for (j = n_multi = 0; j < p->n_multi; ++j) {
        bwt_multi1_t *q = p->multi + j;
        //Get hit position and strand
        q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len + q->ref_shift, &strand);
        q->strand = strand;
        /*q->pos should be different to the main alignment position
          and not -1 (duh!!)
        */
        if (q->pos != p->pos && q->pos != (bwtint_t)-1)
            p->multi[n_multi++] = *q;
    }
    p->n_multi = n_multi;
}

void bwa_refine_gapped1(const bntseq_t *bns, bwa_seq_t *s, const ubyte_t *pacseq)
{
    int j, k, nm;
    if (s->type == BWA_TYPE_NO_MATCH || s->type == BWA_TYPE_MATESW )
        return;
    kstring_t str = {0};
    //We need to complement
    for (j = 0; j < s->len; ++j)
        s->seq[j] = s->seq[j] > 3? 4 : 3 - s->seq[j];
    seq_reverse(s->len, s->seq, 0); // IMPORTANT: s->seq is reversed here!!!
    //Main hit
    if (s->n_gapo) {
        s->cigar = bwa_refine_gapped_core(bns->l_pac,
                                          pacseq,
                                          s->len,
                                          s->strand? s->rseq : s->seq,
                                          s->ref_shift,
                                          &s->pos,
                                          &s->n_cigar);
        // this squeezes out gapped alignments which failed the CIGAR generation
        if (s->cigar == 0) {
            s->type = BWA_TYPE_NO_MATCH;
            return;
        }
    }
    // generate MD tag
    s->md = bwa_cal_md1(s->n_cigar,
                        s->cigar,
                        s->len,
                        s->pos,
                        s->strand? s->rseq : s->seq,
                        bns->l_pac,
                        pacseq,
                        &str,
                        &nm);
    s->nm = nm;
    //Other hits
    for (j = k = 0; j < s->n_multi; ++j) {
        bwt_multi1_t *q = s->multi + j;
        int n_cigar = 0;
        nm = 0;
        if (q->gap) { //q->gap == bwt_aln1_t->n_gapo + bwt_aln1_t->n_gape
            q->cigar = bwa_refine_gapped_core(bns->l_pac,
                                              pacseq,
                                              s->len,
                                              q->strand? s->rseq : s->seq,
                                              q->ref_shift,
                                              &q->pos,
                                              &n_cigar);
            if (!q->cigar) continue;
        }
        q->n_cigar = n_cigar;
        q->md = bwa_cal_md1(q->n_cigar,
                            q->cigar,
                            s->len,
                            q->pos,
                            q->strand ? s->rseq : s->seq,
                            bns->l_pac,
                            pacseq,
                            &str,
                            &nm);
        s->multi[k++] = *q;
    }
    s->n_multi = k;
    free(str.s);
    // correct for trimmed reads
    bwa_correct_trimmed(s);
}

typedef struct {
    bwt_width_t *w;
    bwt_width_t *seed_w;
    gap_stack_t *stack;
} fordata_t;

static fordata_t *fordata_init(const gap_opt_t *opt)
{
    fordata_t *fdata = calloc(opt->n_threads, sizeof(fordata_t));
    for (uint32_t i = 0; i < opt->n_threads; i++) {
        fordata_t *data = &fdata[i];
        data->stack     = gap_init_stack(opt->max_diff,
                                         opt->max_gapo,
                                         opt->max_gape,
                                         opt);
        data->seed_w    = calloc(opt->seed_len+1, sizeof(bwt_width_t));
        data->w         = calloc(opt->max_len+1, sizeof(bwt_width_t));
    }
    return fdata;
}

static void fordata_destroy(fordata_t *fdata, uint32_t n)
{
    for (uint32_t i = 0; i < n; i++) {
        fordata_t data = fdata[i];
        free(data.seed_w); free(data.w);
        gap_destroy_stack(data.stack);
    }
    free(fdata);
}

typedef struct {
    const gap_opt_t *opt;
    int64_t     id;
    bwt_t       *bwt;
    bntseq_t    *bns;
    bwa_seqio_t *ks;
    ubyte_t     *pacseq;
    uint64_t    tseq;
    uint32_t    n_occ; //Report up to n_occ hist
    uint8_t     no_aln;
    void        *forpool;
    //Per for thread data
    fordata_t   *fdata;
} pipeline_t;

typedef struct {
    const pipeline_t *p;
    bwa_seq_t *seqs;
    uint64_t nseqs;
} step_t;

static void worker_for(void *data, long i, int tid)
{
    step_t *t = (step_t *)data;
    fordata_t *fdata = (fordata_t *)t->p->fdata;
    //Compute SA coordinates //TODO Understand
    bwa_cal_sa_reg_gap1(t->p->bwt,
                        &t->seqs[i],
                        t->p->opt,
                        t->p->n_occ,
                        fdata[tid].stack,
                        fdata[tid].seed_w,
                        fdata[tid].w);
    //Get alignment position and strand.
    bwa_cal_pac_pos1(t->p->bns,
                     t->p->bwt,
                     &t->seqs[i],
                     t->p->opt->max_diff,
                     t->p->opt->fnr);
    bwa_refine_gapped1(t->p->bns, &t->seqs[i], t->p->pacseq);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
    pipeline_t *p = (pipeline_t*)shared;
    const gap_opt_t *opt = p->opt;
    if (0 == step) {      //Load sequences
        bwa_seq_t *seqs;
        int nseqs;
        bwa_seqio_t *ks = p->ks;
        seqs = bwa_read_seq(ks, 0x200000, &nseqs, opt->mode, opt->trim_qual);
        if ( 0 < nseqs) {
            step_t *t = calloc(1, sizeof(step_t));
            t->p = p;
            t->seqs = seqs;
            t->nseqs = nseqs;
            return t;
        }
    }
    else if (1 == step) { //Compute alignments
        step_t *t = (step_t *)in;
        kt_forpool(p->forpool, worker_for, in, t->nseqs);
        return t;
    }
    else if (2 == step) { //Write to output
        step_t *t = (step_t *)in;
        p->tseq += t->nseqs;
        fprintf(stderr, "[bwa_alnse_core] %ld sequences processed\n", p->tseq);
        int no_aln = p->no_aln;
        for (uint64_t i = 0; i < t->nseqs; ++i) {
            bwa_seq_t *seq = t->seqs + i;
            //Print only mapped sequences should be here
            if ( no_aln && seq->type == BWA_TYPE_NO_MATCH) continue;
            bwa_print_sam1(p->bns, seq, 0, opt->mode, opt->max_top2);
            fflush(stdout);
        }
        bwa_free_read_seq(t->nseqs, t->seqs);
        free(t);
    }
    return 0;
}

static void bwa_alnse_core(const char *prefix,
                           const char *fn_fa,
                           const gap_opt_t *opt,
                           int n_occ,
                           const char *rg_line,
                           int no_aln)
{
    bwa_seqio_t *ks;
    clock_t t;
    bwt_t     *bwt;
    bntseq_t  *bns;
    ubyte_t   *pacseq;
    void      *forpool;
    fordata_t *fdata;
    // initialization
    ks = bwa_open_reads(opt->mode, fn_fa);
    { // load BWT and BNS
        char *str = (char*)calloc(strlen(prefix) + 10, 1);
        strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
        strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
        bns = bns_restore(prefix);
        srand48(bns->seed);
        free(str);
    }
    { // Load packed sequence
        pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
        err_rewind(bns->fp_pac);
        err_fread_noeof(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
    }
    { // Create ktfor pool and data
        forpool = kt_forpool_init(opt->n_threads);
        //tdata = calloc(opt->n_threads, 4);
        fdata   = fordata_init(opt);
    }
    bwase_initialize();
    // Initiate piepeline object
    pipeline_t p = {0};
    p.opt = opt, p.id = 0, p.ks = ks, p.bwt = bwt;
    p.n_occ = n_occ, p.bns = bns, p.no_aln = no_aln;
    p.pacseq = pacseq, p.forpool = forpool, p.fdata = fdata;

    t = clock();
    bwa_print_sam_hdr(bns, rg_line);
    kt_pipeline(3, worker_pipeline, &p, 3);
    fprintf(stderr, "[%s] Processed %lu sequences in:\n\t", __func__, p.tseq);
    fprintf(stderr, "%.3f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

    bwt_destroy(bwt);
    bwa_seq_close(ks);
    bns_destroy(bns);
    free(pacseq);
    kt_forpool_destroy(forpool);
    fordata_destroy(fdata, opt->n_threads);
}

#define OPTSTR "n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:b012IYB:J:r:u"

int bwa_alnse(int argc, char *argv[])
{
    int c, opte = -1, n_occ = 5, no_aln = 0;
    gap_opt_t *opt;
    char *prefix, *rg_line = 0;

    opt = gap_init_opt();
    while ((c = getopt(argc, argv, OPTSTR)) >= 0) {
        switch (c) {
        case 'n':
            if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
            else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
            break;
        case 'o': opt->max_gapo = atoi(optarg);       break;
        case 'e': opte = atoi(optarg);                break;
        case 'M': opt->s_mm = atoi(optarg);           break;
        case 'O': opt->s_gapo = atoi(optarg);         break;
        case 'E': opt->s_gape = atoi(optarg);         break;
        case 'd': opt->max_del_occ = atoi(optarg);    break;
        case 'i': opt->indel_end_skip = atoi(optarg); break;
        case 'l': opt->seed_len = atoi(optarg);       break;
        case 'k': opt->max_seed_diff = atoi(optarg);  break;
        case 'm': opt->max_entries = atoi(optarg);    break;
        case 't': opt->n_threads = atoi(optarg);      break;
        case 'L': opt->mode |= BWA_MODE_LOGGAP;       break;
        case 'R': opt->max_top2 = atoi(optarg);       break;
        case 'q': opt->trim_qual = atoi(optarg);      break;
        case 'N': opt->mode |= BWA_MODE_NONSTOP;
                  opt->max_top2 = 0x7fffffff;         break;
        case 'f': xreopen(optarg, "wb", stdout);      break;
        case 'b': opt->mode |= BWA_MODE_BAM;          break;
        case '0': opt->mode |= BWA_MODE_BAM_SE;       break;
        case '1': opt->mode |= BWA_MODE_BAM_READ1;    break;
        case '2': opt->mode |= BWA_MODE_BAM_READ2;    break;
        case 'I': opt->mode |= BWA_MODE_IL13;         break;
        case 'Y': opt->mode |= BWA_MODE_CFY;          break;
        case 'B': opt->mode |= atoi(optarg) << 24;    break;
        case 'J': n_occ = atoi(optarg);               break;
        case 'a': opt->max_len = atoi(optarg);        break;
        case 'r':
            if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
            break;
        case 'u': no_aln = 1;                         break;
        default: return 1;
        }
    }
    //TODO check for correct value of opt->max_len
    if (opte > 0) {
        opt->max_gape = opte;
        opt->mode &= ~BWA_MODE_GAPE;
    }
    if (optind + 2 > argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   bwa alnse [options] <prefix> <in.fq>\n\n");
        fprintf(stderr, "Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
                BWA_AVG_ERR, opt->fnr);
        fprintf(stderr, "         -o INT    maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
        fprintf(stderr, "         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
        fprintf(stderr, "         -i INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
        fprintf(stderr, "         -d INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
        fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
        fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);
        fprintf(stderr, "         -m INT    maximum entries in the queue [%d]\n", opt->max_entries);
        fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
        fprintf(stderr, "         -M INT    mismatch penalty [%d]\n", opt->s_mm);
        fprintf(stderr, "         -O INT    gap open penalty [%d]\n", opt->s_gapo);
        fprintf(stderr, "         -E INT    gap extension penalty [%d]\n", opt->s_gape);
        fprintf(stderr, "         -R INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
        fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
        fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");
        fprintf(stderr, "         -B INT    length of barcode\n");
        fprintf(stderr, "         -L        log-scaled gap penalty for long deletions\n");
        fprintf(stderr, "         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
        fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
        fprintf(stderr, "         -b        the input read file is in the BAM format\n");
        fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
        fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
        fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
        fprintf(stderr, "         -Y        filter Casava-filtered sequences\n");
        fprintf(stderr, "         -J INT    Report up to INT additional hits in XA tag [5]\n");
        fprintf(stderr, "         -a INT    Maximum read length [255]\n");
        fprintf(stderr, "         -r STR    RG_line\n");
        fprintf(stderr, "         -u        Do not output unmapped reads\n");
        fprintf(stderr, "\n");
        return 1;
    }
    if (opt->fnr > 0.0) {
        int i, k;
        for (i = 17, k = 0; i <= 250; ++i) {
            int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
            if (l != k) fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l);
            k = l;
        }
    }
    if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
        fprintf(stderr, "[%s] fail to locate the index\n", __func__);
        free(opt);
        return 1;
    }
    bwa_alnse_core(prefix, argv[optind+1], opt, n_occ, rg_line, no_aln);
    free(opt); free(prefix);
    return 0;
}
