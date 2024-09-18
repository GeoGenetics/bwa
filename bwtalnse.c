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

//From bwtaln.c
bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);

void bwa_cal_sa_reg_gap1(bwt_t *const bwt,
                         bwa_seq_t *p,
                         const gap_opt_t *opt)
{
    int j, max_l = 0, max_len = p->len;
    gap_stack_t *stack;
    bwt_width_t *w, *seed_w;
    gap_opt_t local_opt = *opt;
    // initiate priority stack
    if (opt->fnr > 0.0)
        local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
    if (local_opt.max_diff < local_opt.max_gapo)
        local_opt.max_gapo = local_opt.max_diff;
    stack = gap_init_stack(local_opt.max_diff,
                           local_opt.max_gapo,
                           local_opt.max_gape,
                           &local_opt);
    seed_w = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
    w = 0;
    p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
    if (max_l < p->len) {
        max_l = p->len;
        w = (bwt_width_t*)realloc(w, (max_l + 1) * sizeof(bwt_width_t));
        memset(w, 0, (max_l + 1) * sizeof(bwt_width_t));
    }
    bwt_cal_width(bwt, p->len, p->seq, w);
    if (opt->fnr > 0.0)
        local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
    local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
    if (p->len > opt->seed_len)
        bwt_cal_width(bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), seed_w);
    // core function
    for (j = 0; j < p->len; ++j) // we need to complement
        p->seq[j] = p->seq[j] > 3? 4 : 3 - p->seq[j];
    p->aln = bwt_match_gap(bwt,
                           p->len,
                           p->seq,
                           w,
                           p->len <= opt->seed_len? 0 : seed_w,
                           &local_opt,
                           &p->n_aln,
                           stack);
    // clean up the unused data in the record
    free(p->name); free(p->seq); free(p->rseq); free(p->qual);
    p->name = 0; p->seq = p->rseq = p->qual = 0;
    free(seed_w); free(w);
    gap_destroy_stack(stack);
}

typedef struct {
    const gap_opt_t *opt;
    int64_t     id;
    bwt_t       *bwt;
    bwa_seqio_t *ks;
    uint64_t    tseq;
} pipeline_t;

typedef struct {
    const pipeline_t *p;
    bwa_seq_t *seqs;
    uint64_t nseqs;
} step_t;

static void worker_for(void *data, long i, int tid)
{
    step_t *t = (step_t *)data;
    bwa_cal_sa_reg_gap1(t->p->bwt, &t->seqs[i], t->p->opt);
}

static void *worker_pipeline(void *shared, int step, void *in)
{
    pipeline_t *p = (pipeline_t*)shared;
    const gap_opt_t *opt = p->opt;
    step_t *t = (step_t *)in;
    //int32_t i;
    if (0 == step) { //Load sequences
        bwa_seq_t *seqs;
        int nseqs;
        bwa_seqio_t *ks = p->ks;
        seqs = bwa_read_seq(ks, 0x200000, &nseqs, opt->mode, opt->trim_qual);
        fprintf(stderr, "[%s] Loaded %d sequences\n", __func__, nseqs);
        p->tseq += nseqs;
        if ( 0 < nseqs) {
            if (t) free(t);
            t = calloc(1, sizeof(step_t));
            t->p = p;
            t->seqs = seqs;
            t->nseqs = nseqs;
            return t;
        }
    }
    else if (1 == step) {
        fprintf(stderr, "[%s] Calculating SA coordinates\n", __func__);
        kt_for(p->opt->n_threads, worker_for, in, t->nseqs);
    }
    else if (2 == step) {
        
    }
    return 0;
}

static void bwa_alnse_core(const char *prefix,
                           const char *fn_fa,
                           const gap_opt_t *opt,
                           int n_occ,
                           const char *rg_line)
{
    fprintf(stderr, "%s %d %s\n", __func__, n_occ, rg_line);
    bwa_seqio_t *ks;
    clock_t t;
    bwt_t *bwt;
    // initialization
    ks = bwa_open_reads(opt->mode, fn_fa);
    { // load BWT
        char *str = (char*)calloc(strlen(prefix) + 10, 1);
        strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
        free(str);
    }
    // Initiate piepeline object
    pipeline_t p = {0};
    p.opt = opt, p.id = 0, p.ks = ks, p.bwt = bwt;
    t = clock();
    kt_pipeline(2, worker_pipeline, &p, 3);
    fprintf(stderr, "[%s] Processed %lu sequences in:\n\t", __func__, p.tseq);
    fprintf(stderr, "%.3f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    /*
		t = clock();
		fprintf(stderr, "[bwa_aln_core] write to the disk... ");
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			err_fwrite(&p->n_aln, 4, 1, stdout);
			if (p->n_aln) err_fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[bwa_aln_core] %lld sequences have been processed.\n", tot_seqs);
	}

	// destroy
	bwt_destroy(bwt);
	bwa_seq_close(ks);
    */
}

#define OPTSTR "n:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:b012IYB:J:r:"

int bwa_alnse(int argc, char *argv[])
{
    int c, opte = -1, n_occ = 5;
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
        case 'r':
            if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
            break;
        default: return 1;
        }
    }
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
        fprintf(stderr, "         -r STR    RG_line\n");
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
    bwa_alnse_core(prefix, argv[optind+1], opt, n_occ, rg_line);
    free(opt); free(prefix);
    fprintf(stderr, "%s\n", __func__);
    return 0;
}
