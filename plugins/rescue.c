/*  plugins/rescue.c -- predict the rescue of protein-truncating variants.

    The variants of a phased VCF that were annotated with the consequences of
    Ensembl VEP or bcftools csq are checked for splice changes with SpliceAI-style
    models, per haplotype. A protein-truncating variant (stop gained, frameshift,
    splice donor or acceptor) that is rescued by another variant of the haplotype,
    or that induces an exon skip or an intron retention, is tagged, together with
    whether the transcript is predicted to be degraded by nonsense-mediated decay.

    The predictions are made by librescueranger, which is not part of bcftools:
    this plugin is only compiled when the library is installed, see INSTALL.
    This file takes care of the command line and the files. */

#include <getopt.h>
#include <limits.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <rescueranger/rescueranger.h>
#include "bcftools.h"

typedef struct
{
    char *fname;                // input VCF
    char *output_fname, *scores_fname;
    int output_type, clevel;
    int record_cmd_line;
    RescueRangerOptions opts;
    int argc;
    char **argv;
}
args_t;

const char *about(void)
{
    return "Predict the rescue of protein-truncating variants from splice changes.\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Predict the rescue of protein-truncating variants (PTVs) by other variants of the same\n"
        "       haplotype, using SpliceAI-style models. The VCF must be phased and annotated with the\n"
        "       consequences by Ensembl VEP (INFO/CSQ) or bcftools csq (INFO/BCSQ). Adds the INFO tags\n"
        "       with the rescued and induced PTVs and the flags of nonsense-mediated decay (NMD).\n"
        "\n"
        "Usage: bcftools +rescue [Options] <in.vcf> -m <models> -g <annot.gff> -f <ref.fa>\n"
        "Plugin options:\n"
        "   -f, --fasta-ref FILE            Reference sequence in fasta format\n"
        "   -g, --gff-annot FILE            GFF3 annotation matching the VCF, see the csq command for the format\n"
        "   -m, --models DIR                Directory with the SpliceAI models\n"
        "   -o, --output FILE               Write output to a file [standard output]\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "       --af-tag TAG                Allele frequency tag [INFO/AF]\n"
        "       --af-threshold FLOAT        Allele frequency threshold [0.01]\n"
        "       --annotation vep|bcftools   Program that annotated the consequences of the variants [vep]\n"
        "       --csq-tag TAG               INFO field with the consequences [CSQ for vep, BCSQ for bcftools]\n"
        "       --exclude-unphased          Leave out heterozygous variants with an unknown phase (by default they are on both haplotypes)\n"
        "       --gt-tag TAG                Phased genotype tag [FORMAT/GT]\n"
        "       --no-version                Do not append version and command line to the header\n"
        "       --prediction-radius INT     Prediction window radius in bases; other values than 50 give scores not comparable to SpliceAI [50]\n"
        "       --raw-scores FILE           Write the raw SpliceAI scores to a TSV file\n"
        "       --reference-scores FILE     Precomputed reference scores (cpliceai_reference). Optional: without them, and for genes that\n"
        "                                       are not in them, the reference predictions are made from the sequence while running\n"
        "\n"
        "Example:\n"
        "   bcftools +rescue in.vcf.gz -m models/ -g annot.gff3.gz -f ref.fa -o out.vcf --raw-scores scores.tsv\n"
        "\n";
}

static int annotate(args_t *args)
{
    int ret = 1;
    htsFile *fp_in = NULL, *fp_out = NULL;
    bcf_hdr_t *hdr = NULL;
    FILE *scores_tsv = NULL;
    RescueRanger *rr = NULL;

    if ( !(fp_in = hts_open(args->fname, "r")) )
    {
        fprintf(stderr, "Could not open variants file: %s\n", args->fname);
        goto done;
    }
    if ( !(hdr = bcf_hdr_read(fp_in)) )
    {
        fprintf(stderr, "Could not read VCF header from: %s\n", args->fname);
        goto done;
    }

    // Loads the models, the annotation and the reference, and adds the INFO tags of the output to the header
    if ( !(rr = rr_init(&args->opts, hdr)) ) goto done;

    char wmode[8];
    set_wmode(wmode, args->output_type, args->output_fname, args->clevel);
    if ( !(fp_out = hts_open(args->output_fname, wmode)) )
    {
        fprintf(stderr, "Could not open annotated VCF for writing: %s\n", args->output_fname);
        goto done;
    }
    if ( args->record_cmd_line ) bcf_hdr_append_version(hdr, args->argc, args->argv, "bcftools_plugin");
    if ( bcf_hdr_write(fp_out, hdr)!=0 )
    {
        fprintf(stderr, "Could not write VCF header to: %s\n", args->output_fname);
        goto done;
    }

    if ( args->scores_fname )
    {
        if ( !(scores_tsv = fopen(args->scores_fname, "w")) )
        {
            fprintf(stderr, "Could not open raw scores TSV for writing: %s\n", args->scores_fname);
            goto done;
        }
        if ( rr_set_scores_output(rr, scores_tsv)!=0 ) goto done;
    }

    if ( rr_annotate(rr, fp_in, fp_out)==0 ) ret = 0;

done:
    rr_destroy(rr);
    if ( scores_tsv && fclose(scores_tsv)!=0 ) ret = 1;
    if ( fp_out && hts_close(fp_out)!=0 ) ret = 1;
    if ( hdr ) bcf_hdr_destroy(hdr);
    if ( fp_in ) hts_close(fp_in);
    return ret;
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1, sizeof(args_t));
    args->argc = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type  = FT_VCF;
    args->record_cmd_line = 1;
    args->clevel = -1;
    rr_opts_init(&args->opts);
    static struct option loptions[] =
    {
        {"fasta-ref",required_argument,NULL,'f'},
        {"gff-annot",required_argument,NULL,'g'},
        {"models",required_argument,NULL,'m'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"af-tag",required_argument,NULL,1},
        {"af-threshold",required_argument,NULL,2},
        {"annotation",required_argument,NULL,3},
        {"csq-tag",required_argument,NULL,4},
        {"exclude-unphased",no_argument,NULL,5},
        {"gt-tag",required_argument,NULL,6},
        {"no-version",no_argument,NULL,7},
        {"prediction-radius",required_argument,NULL,8},
        {"raw-scores",required_argument,NULL,9},
        {"reference-scores",required_argument,NULL,10},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    long radius;
    while ((c = getopt_long(argc, argv, "f:g:m:o:O:h?", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'f': args->opts.fasta_fname = optarg; break;
            case 'g': args->opts.gff_fname = optarg; break;
            case 'm': args->opts.model_dir = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default:
                    {
                        args->clevel = strtol(optarg,&tmp,10);
                        if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                    }
                }
                if ( optarg[1] )
                {
                    args->clevel = strtol(optarg+1,&tmp,10);
                    if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                }
                break;
            case  1 : args->opts.af_tag = optarg; break;
            case  2 :
                args->opts.af_threshold = strtod(optarg, &tmp);
                if ( tmp==optarg || *tmp || args->opts.af_threshold<0 || args->opts.af_threshold>1 ) error("Could not parse: --af-threshold %s\n", optarg);
                break;
            case  3 :
                if ( rr_annotation_parse(optarg, &args->opts.annotation)!=0 ) error("--annotation must be vep or bcftools, got %s\n", optarg);
                break;
            case  4 : args->opts.csq_tag = optarg; break;
            case  5 : args->opts.exclude_unphased = 1; break;
            case  6 : args->opts.gt_tag = optarg; break;
            case  7 : args->record_cmd_line = 0; break;
            case  8 :
                radius = strtol(optarg, &tmp, 10);
                if ( tmp==optarg || *tmp || radius<1 || radius>INT_MAX ) error("Could not parse: --prediction-radius %s\n", optarg);
                args->opts.prediction_radius = (int) radius;
                break;
            case  9 : args->scores_fname = optarg; break;
            case 10 : args->opts.reference_scores = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s", usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s", usage_text());
    else args->fname = argv[optind];

    if ( !args->opts.model_dir ) error("Missing the -m option\n");
    if ( !args->opts.gff_fname ) error("Missing the -g option\n");
    if ( !args->opts.fasta_fname ) error("Missing the -f option\n");

    rr_set_log_level(RR_LOG_INFO);
    // TensorFlow logs independently of the library; 1 drops its INFO output
    setenv("TF_CPP_MIN_LOG_LEVEL", "1", 1);

    int ret = annotate(args);
    free(args);
    return ret;
}
