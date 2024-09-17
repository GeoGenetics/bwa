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

static void bwa_alnse_core(const char *prefix,
                           const char *fn_fa,
                           const gap_opt_t *opt,
                           int n_occ,
                           const char *rg_line)
{
    fprintf(stderr, "%s %d %s\n", __func__, n_occ, rg_line);
    sleep(100);
    int i, n_seqs;
    long long tot_seqs = 0;
    bwa_seq_t *seqs;
    bwa_seqio_t *ks;
    clock_t t;
    bwt_t *bwt;
    /*
	// initialization
	ks = bwa_open_reads(opt->mode, fn_fa);

	{ // load BWT
		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
		free(str);
	}

	// core loop
	err_fwrite(SAI_MAGIC, 1, 4, stdout);
	err_fwrite(opt, sizeof(gap_opt_t), 1, stdout);
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();

		fprintf(stderr, "[bwa_aln_core] calculate SA coordinate... ");

#ifdef HAVE_PTHREAD
		if (opt->n_threads <= 1) { // no multi-threading at all
			bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
		} else {
			pthread_t *tid;
			pthread_attr_t attr;
			thread_aux_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
			for (j = 0; j < opt->n_threads; ++j) {
				data[j].tid = j; data[j].bwt = bwt;
				data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
				pthread_create(&tid[j], &attr, worker, data + j);
			}
			for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else
		bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
#endif

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

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
