/*  plugins/lof.c -- LOFTEE confidence of the loss-of-function consequences of bcftools csq.

    A reimplementation of the LOFTEE plugin of Ensembl VEP (v1.0.4_GRCh38,
    https://github.com/konradjk/loftee) for the consequences in INFO/BCSQ. Each
    stop_gained, frameshift, splice_donor and splice_acceptor consequence of a
    protein-coding transcript is high confidence (HC) unless one of the LOFTEE
    filters applies, which makes it low confidence (LC). The flags are only
    informative. They go to INFO/LOF, INFO/LOF_FILTERS and INFO/LOF_FLAGS.

    Where this differs from LOFTEE:
        - the transcripts come from the GFF3 given to bcftools csq, Ensembl or
          GENCODE, rather than from the Ensembl API
        - END_TRUNC follows LOFTEE-2 (the relaxed cutoffs of the gnomAD v4
          flagship paper) instead of the GERP-weighted distance to the stop
          codon, and START_TRUNC filters a truncation too close to the start of
          the CDS, as LOFTEE-2 does
        - for a frameshift, START_TRUNC is measured from the new stop codon,
          located from the protein predicted by csq, and STOP_END_TRUNC filters
          a frameshift whose new stop falls where END_TRUNC would apply
        - there are no splice predictions (MaxEntScan, de novo donors), so no
          OS consequences
        - the PhyloCSF, ancestral allele and incomplete CDS checks are always
          on. A CDS is incomplete if the GFF phase trims it or if the reference
          lacks its start or stop codon
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <getopt.h>
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"
#include "regidx.h"
#include "gff.h"

// The relaxed LOFTEE-2 cutoffs: a stop gain or a frameshift must be further than this many bases from
#define MIN_DIST_LAST_EXON      50      // the start of the last coding exon (the 50 bp rule)
#define MIN_DIST_CDS_START      50      // the start of the CDS
#define MIN_DIST_CDS_START_1EX 250      // the start of the CDS, single-exon transcripts
#define MIN_DIST_TSCRIPT_END   200      // the end of the transcript, single-exon transcripts

#define MIN_INTRON_SIZE 15

#define LOF_STOP_GAINED  1
#define LOF_FRAMESHIFT   2
#define LOF_DONOR        4
#define LOF_ACCEPTOR     8

#define FLAG   0
#define FILTER 1

typedef struct
{
    uint32_t beg, end;      // 0-based, inclusive
    int tbeg;               // offset of the exon's first base in the spliced transcript
}
exon_t;

// Exons and CDS of a transcript in the transcript orientation, made on first use and kept in gf_tscript_t.aux
typedef struct
{
    int nexon, tlen, ok;
    exon_t *exon;
    int cds_tbeg, cds_tend;     // the CDS as transcript offsets, including the stop codon
    int cds_ibeg, cds_iend;     // the first and the last coding exon
    uint32_t cds_beg, cds_end;  // the CDS in genomic coordinates
    int incomplete;
}
tscript_t;

typedef struct
{
    int nexon;
    float *orf, *max;       // corresponding_orf_score and max_score by exon, NAN if missing
}
phylocsf_t;

typedef struct
{
    uint32_t beg, end;      // reference bases changed by the variant, both flanking bases of an insertion
    uint32_t pos;           // the first changed base, the position LOFTEE measures distances from
    char ref, alt;          // alleles of an SNV, 0 otherwise
}
variant_t;

// One INFO/BCSQ consequence: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change
typedef struct
{
    int type;               // LOF_* bits
    char *tscript, *aa, *dna;
}
bcsq_t;

typedef struct
{
    int off, dlen;          // CDS offset and length change of an indel
}
indel_t;

typedef struct
{
    bcf_hdr_t *hdr_in, *hdr_out;
    char *csq_tag;
    gff_t *gff;
    regidx_t *idx_tscript, *idx_exon;
    regitr_t *itr_tscript, *itr_exon;
    faidx_t *ref, *anc;
    void *phylo_idx;        // transcript id without the version -> index to phylo
    phylocsf_t *phylo;
    int nphylo, mphylo;
    int rid;                // the current chromosome and its names in the GFF and the fasta files
    const char *chr_gff, *chr_ref, *chr_anc;
    kstring_t chr_gff_str, chr_ref_str, chr_anc_str;
    char *bcsq;
    int mbcsq;
    kstring_t lof, filters, flags;      // INFO tags of the current record
    kstring_t csq_filters, csq_flags;   // of the current consequence
    kstring_t tmp;
    indel_t *indel;
    int nindel, mindel;
    int nmultiallelic, nunknown_tscript, nno_ancestral;
}
args_t;

static args_t args;

const char *about(void)
{
    return "LOFTEE confidence (HC/LC) of the loss-of-function consequences of bcftools csq.\n";
}

const char *usage(void)
{
    return
        "\n"
        "About: Assess the loss-of-function consequences of bcftools csq like the LOFTEE plugin of VEP.\n"
        "       Adds INFO/LOF with HC or LC, INFO/LOF_FILTERS with the LOFTEE filters that make it LC and\n"
        "       INFO/LOF_FLAGS with the informative LOFTEE flags, one value per INFO/BCSQ consequence and\n"
        "       '.' for consequences that are not a loss of function or have no filters or flags.\n"
        "       Works with local and haplotype-aware csq output; a compound consequence is assessed\n"
        "       where csq prints it, and the \"@POS\" references to it are '.'. Multiallelic records are\n"
        "       skipped, split them before running bcftools csq.\n"
        "Usage: bcftools +lof [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -a, --ancestral FILE        Ancestral sequence indexed by samtools faidx, e.g. human_ancestor.fa.gz\n"
        "   -c, --custom-tag STRING     The INFO tag written by bcftools csq [BCSQ]\n"
        "   -f, --fasta-ref FILE        Reference sequence in fasta format\n"
        "       --force                 Skip transcripts with an inconsistent CDS phase, as bcftools csq --force\n"
        "   -g, --gff-annot FILE        The GFF3 annotation given to bcftools csq, Ensembl or GENCODE\n"
        "   -p, --phylocsf FILE         PhyloCSF scores of the exons, tab-delimited with a header naming the\n"
        "                               columns transcript, exon, corresponding_orf_score and max_score\n"
        "\n"
        "Example:\n"
        "   # the PhyloCSF scores from the SQLite database of LOFTEE\n"
        "   sqlite3 -header -separator $'\\t' loftee.sql 'SELECT * FROM phylocsf_summary' | bgzip > phylocsf.tsv.gz\n"
        "\n"
        "   bcftools csq -l -f ref.fa -g annot.gff3.gz in.vcf.gz -Ou | \\\n"
        "       bcftools +lof -Oz -o out.vcf.gz -- -f ref.fa -g annot.gff3.gz -a human_ancestor.fa.gz -p phylocsf.tsv.gz\n"
        "\n";
}

static inline char complement(char c)
{
    switch (c)
    {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
    }
    return 'N';
}

// Fetch the bases beg..end (0-based, inclusive) in upper case, reverse complemented for the reverse
// strand. The output must have room for end-beg+2 bytes. Returns -1 if not all bases are in the fasta
static int fetch_seq(faidx_t *fai, const char *chr, uint32_t beg, uint32_t end, int strand, char *out)
{
    hts_pos_t i, len;
    char *seq = faidx_fetch_seq64(fai, chr, beg, end, &len);
    if ( !seq || len != end - beg + 1 ) { free(seq); return -1; }
    for (i=0; i<len; i++)
    {
        char base = toupper((unsigned char)seq[strand==STRAND_REV ? len-i-1 : i]);
        out[i] = strand==STRAND_REV ? complement(base) : base;
    }
    out[len] = 0;
    free(seq);
    return 0;
}

static inline int is_stop(const char *codon)
{
    return !strcmp(codon,"TAA") || !strcmp(codon,"TAG") || !strcmp(codon,"TGA");
}

// Transcript offset of a genomic position, -1 if not in an exon
static int tscript_offset(gf_tscript_t *tr, tscript_t *tx, uint32_t pos, int *iexon)
{
    int i;
    for (i=0; i<tx->nexon; i++)
    {
        exon_t *exon = &tx->exon[i];
        if ( pos < exon->beg || pos > exon->end ) continue;
        if ( iexon ) *iexon = i;
        return exon->tbeg + (tr->strand==STRAND_FWD ? pos - exon->beg : exon->end - pos);
    }
    return -1;
}

// The exon of a transcript offset, 0 <= toff < tlen
static int tscript_exon(tscript_t *tx, int toff)
{
    int i = tx->nexon - 1;
    while ( i>0 && tx->exon[i].tbeg > toff ) i--;
    return i;
}

// Genomic position of a transcript offset
static uint32_t tscript_pos(gf_tscript_t *tr, tscript_t *tx, int toff, int *iexon)
{
    int i = tscript_exon(tx, toff);
    if ( iexon ) *iexon = i;
    exon_t *exon = &tx->exon[i];
    return tr->strand==STRAND_FWD ? exon->beg + (toff - exon->tbeg) : exon->end - (toff - exon->tbeg);
}

// Bases toff..toff+len-1 of the spliced transcript
static int tscript_seq(gf_tscript_t *tr, tscript_t *tx, int toff, int len, char *out)
{
    if ( toff < 0 || toff + len > tx->tlen ) return -1;
    int n = 0;
    while ( n < len )
    {
        int i;
        uint32_t pos = tscript_pos(tr, tx, toff + n, &i);
        exon_t *exon = &tx->exon[i];
        int nbases = exon->tbeg + (exon->end - exon->beg + 1) - (toff + n);
        if ( nbases > len - n ) nbases = len - n;
        uint32_t beg = tr->strand==STRAND_FWD ? pos : pos - nbases + 1;
        if ( fetch_seq(args.ref, args.chr_ref, beg, beg + nbases - 1, tr->strand, out + n) < 0 ) return -1;
        n += nbases;
    }
    return 0;
}

static tscript_t *tscript_get(gf_tscript_t *tr)
{
    if ( tr->aux ) return (tscript_t*) tr->aux;

    tscript_t *tx = calloc(1, sizeof(*tx));
    tr->aux = tx;
    if ( !tr->ncds || !tr->nexons ) return tx;

    // gff.c numbers the exons in the transcript orientation
    tx->nexon = tr->nexons;
    tx->exon  = calloc(tx->nexon, sizeof(*tx->exon));
    int i, n = 0;
    if ( regidx_overlap(args.idx_exon, args.chr_gff, tr->beg, tr->end, args.itr_exon) )
    {
        while ( regitr_overlap(args.itr_exon) )
        {
            gf_exon_t *exon = regitr_payload(args.itr_exon, gf_exon_t*);
            if ( exon->tr != tr || exon->iexon >= tx->nexon ) continue;
            tx->exon[exon->iexon].beg = exon->beg;
            tx->exon[exon->iexon].end = exon->end;
            n++;
        }
    }
    if ( n != tx->nexon ) return tx;
    for (i=0; i<tx->nexon; i++)
    {
        tx->exon[i].tbeg = tx->tlen;
        tx->tlen += tx->exon[i].end - tx->exon[i].beg + 1;
    }

    uint32_t beg = UINT32_MAX, end = 0;
    for (i=0; i<tr->ncds; i++)
    {
        gf_cds_t *cds = tr->cds[i];
        if ( !cds->len ) continue;
        if ( beg > cds->beg ) beg = cds->beg;
        if ( end < cds->beg + cds->len - 1 ) end = cds->beg + cds->len - 1;
    }
    int fwd = tr->strand==STRAND_FWD;
    tx->cds_tbeg = tscript_offset(tr, tx, fwd ? beg : end, NULL);
    tx->cds_tend = tscript_offset(tr, tx, fwd ? end : beg, NULL);
    if ( tx->cds_tbeg < 0 || tx->cds_tend < tx->cds_tbeg ) return tx;

    // The CDS includes the stop codon in Ensembl GFFs but not in GENCODE
    char codon[4];
    int has_stop = !tscript_seq(tr, tx, tx->cds_tend - 2, 3, codon) && is_stop(codon);
    if ( !has_stop && !tr->trim && !tscript_seq(tr, tx, tx->cds_tend + 1, 3, codon) && is_stop(codon) )
    {
        tx->cds_tend += 3;
        has_stop = 1;
    }
    int has_start = !tscript_seq(tr, tx, tx->cds_tbeg, 3, codon) && !strcmp(codon,"ATG");
    tx->incomplete = tr->trim || !has_start || !has_stop;

    beg = tscript_pos(tr, tx, tx->cds_tbeg, &tx->cds_ibeg);
    end = tscript_pos(tr, tx, tx->cds_tend, &tx->cds_iend);
    tx->cds_beg = fwd ? beg : end;
    tx->cds_end = fwd ? end : beg;
    tx->ok = 1;
    return tx;
}

static void tscript_destroy(void)
{
    regitr_t *itr = regitr_init(args.idx_tscript);
    while ( regitr_loop(itr) )
    {
        gf_tscript_t *tr = regitr_payload(itr, gf_tscript_t*);
        tscript_t *tx = (tscript_t*) tr->aux;
        if ( !tx ) continue;
        free(tx->exon);
        free(tx);
        tr->aux = NULL;
    }
    regitr_destroy(itr);
}

static gf_tscript_t *find_tscript(const char *id, variant_t *var)
{
    if ( !args.chr_gff || !regidx_overlap(args.idx_tscript, args.chr_gff, var->beg, var->end, args.itr_tscript) ) return NULL;
    while ( regitr_overlap(args.itr_tscript) )
    {
        gf_tscript_t *tr = regitr_payload(args.itr_tscript, gf_tscript_t*);
        if ( !strcmp(id, gff_id2string(args.gff, transcript, tr->id)) ) return tr;
    }
    return NULL;
}

// Split a line on tabs in place. Unlike ksplit, keeps the empty fields
static int split_tabs(char *str, char ***col, int *mcol)
{
    int n = 0;
    while ( 1 )
    {
        hts_expand(char*, n+1, *mcol, *col);
        (*col)[n++] = str;
        if ( !(str = strchr(str,'\t')) ) break;
        *str++ = 0;
    }
    return n;
}

static void phylocsf_init(const char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Failed to read %s\n", fname);

    kstring_t str = {0,0,0};
    char **col = NULL;
    int mcol = 0, nline = 1, i, n;
    int itr = -1, iex = -1, iorf = -1, imax = -1;
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Failed to read %s\n", fname);
    n = split_tabs(str.s, &col, &mcol);
    for (i=0; i<n; i++)
    {
        char *name = col[i];
        if ( !i && *name=='#' ) name++;
        if ( !strcmp(name,"transcript") ) itr = i;
        else if ( !strcmp(name,"exon") ) iex = i;
        else if ( !strcmp(name,"corresponding_orf_score") ) iorf = i;
        else if ( !strcmp(name,"max_score") ) imax = i;
    }
    if ( itr<0 || iex<0 || iorf<0 || imax<0 )
        error("The header of %s must name the columns transcript, exon, corresponding_orf_score and max_score\n", fname);

    args.phylo_idx = khash_str2int_init();
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        nline++;
        n = split_tabs(str.s, &col, &mcol);
        if ( n <= itr || n <= iex || n <= iorf || n <= imax ) error("Too few columns at line %d of %s\n", nline, fname);

        char *id = col[itr], *orf = col[iorf], *max = col[imax], *end;
        int iexon = strtol(col[iex], &end, 10);
        if ( *end || iexon < 1 ) error("Could not parse the exon number at line %d of %s\n", nline, fname);

        // A missing score is an empty, NULL or NA field
        double orf_score = strtod(orf, &end);
        if ( end==orf || *end ) continue;
        double max_score = strtod(max, &end);
        if ( end==max || *end ) max_score = NAN;

        id[strcspn(id,".")] = 0;
        int idx;
        if ( khash_str2int_get(args.phylo_idx, id, &idx) < 0 )
        {
            idx = args.nphylo++;
            hts_expand0(phylocsf_t, args.nphylo, args.mphylo, args.phylo);
            khash_str2int_set(args.phylo_idx, strdup(id), idx);
        }
        phylocsf_t *phylo = &args.phylo[idx];
        if ( iexon > phylo->nexon )
        {
            phylo->orf = realloc(phylo->orf, iexon * sizeof(*phylo->orf));
            phylo->max = realloc(phylo->max, iexon * sizeof(*phylo->max));
            for (i=phylo->nexon; i<iexon; i++) phylo->orf[i] = phylo->max[i] = NAN;
            phylo->nexon = iexon;
        }
        phylo->orf[iexon-1] = orf_score;
        phylo->max[iexon-1] = max_score;
    }
    if ( hts_close(fp)!=0 ) error("Close failed: %s\n", fname);
    free(str.s);
    free(col);
}

static void phylocsf_destroy(void)
{
    int i;
    for (i=0; i<args.nphylo; i++)
    {
        free(args.phylo[i].orf);
        free(args.phylo[i].max);
    }
    free(args.phylo);
    khash_str2int_destroy_free(args.phylo_idx);
}

// The GFF and fasta files may name chromosomes with or without the "chr" prefix, independently of the VCF
static int gff_has(void *gff, const char *chr) { return gff_has_seq((gff_t*)gff, chr); }
static int fai_has(void *fai, const char *chr) { return faidx_has_seq((faidx_t*)fai, chr); }

static const char *find_chr(const char *chr, int (*has_chr)(void*,const char*), void *db, kstring_t *str)
{
    str->l = 0;
    if ( !strcmp(chr,"MT") ) kputs("chrM", str);
    else if ( !strcmp(chr,"chrM") ) kputs("MT", str);
    else if ( !strncmp(chr,"chr",3) ) kputs(chr+3, str);
    else ksprintf(str, "chr%s", chr);

    if ( has_chr(db, chr) ) return chr;
    if ( has_chr(db, str->s) ) return str->s;
    return NULL;
}

static void set_chr(int rid)
{
    if ( args.rid==rid ) return;
    args.rid = rid;
    const char *chr = bcf_hdr_id2name(args.hdr_in, rid);
    args.chr_gff = find_chr(chr, gff_has, args.gff, &args.chr_gff_str);
    args.chr_ref = find_chr(chr, fai_has, args.ref, &args.chr_ref_str);
    args.chr_anc = find_chr(chr, fai_has, args.anc, &args.chr_anc_str);
}

static void variant_init(variant_t *var, bcf1_t *rec)
{
    const char *ref = rec->d.allele[0], *alt = rec->d.allele[1];
    int rlen = strlen(ref), alen = strlen(alt);

    var->ref = var->alt = 0;
    if ( rlen==1 && alen==1 ) var->ref = toupper((unsigned char)*ref), var->alt = toupper((unsigned char)*alt);

    var->beg = rec->pos;
    var->end = rec->pos + rlen - 1;
    if ( rlen!=alen && toupper((unsigned char)*ref)==toupper((unsigned char)*alt) )
    {
        if ( rlen > 1 ) var->beg++;     // skip the padding base of a deletion
        else var->end++;                // an insertion touches the bases on both sides
    }
    var->pos = rlen < alen ? var->end : var->beg;
}

static void add_tag(const char *tag, int is_filter)
{
    kstring_t *str = is_filter ? &args.csq_filters : &args.csq_flags;
    if ( str->l ) kputc('&', str);
    kputs(tag, str);
}

// The first exon in the transcript order that the variant overlaps, -1 if none
static int find_exon(tscript_t *tx, variant_t *var)
{
    int i;
    for (i=0; i<tx->nexon; i++)
        if ( var->beg <= tx->exon[i].end && var->end >= tx->exon[i].beg ) return i;
    return -1;
}

static int find_intron(gf_tscript_t *tr, tscript_t *tx, variant_t *var, uint32_t *beg, uint32_t *end)
{
    int i;
    for (i=0; i+1<tx->nexon; i++)
    {
        exon_t *a = &tx->exon[i], *b = &tx->exon[i+1];
        *beg = tr->strand==STRAND_FWD ? a->end + 1 : b->end + 1;
        *end = tr->strand==STRAND_FWD ? b->beg - 1 : a->beg - 1;
        if ( *beg > *end ) continue;
        if ( var->beg <= *end && var->end >= *beg ) return i;
    }
    return -1;
}

// Transcript offset of a position, or of the first base of the next exon if the position is
// intronic. Returns -1 if there is no exon downstream
static int exonic_offset(gf_tscript_t *tr, tscript_t *tx, uint32_t pos, int *iexon)
{
    int i, toff = tscript_offset(tr, tx, pos, iexon);
    if ( toff >= 0 ) return toff;
    for (i=0; i<tx->nexon; i++)
        if ( tr->strand==STRAND_FWD ? tx->exon[i].beg > pos : tx->exon[i].end < pos ) break;
    if ( i==tx->nexon ) return -1;
    *iexon = i;
    return tx->exon[i].tbeg;
}

/*
    The LOFTEE-2 distance from the start of the last coding exon, negative inside it. Like the
    distance from the start of the CDS, it loses a base for each exon boundary crossed, as in LOFTEE,
    from which LOFTEE-2 took them.
*/
static int dist_last_exon(tscript_t *tx, int toff, int iexon)
{
    int last = tx->cds_iend;
    int last_beg = tx->exon[last].tbeg > tx->cds_tbeg ? tx->exon[last].tbeg : tx->cds_tbeg;
    if ( iexon >= last ) return last_beg - toff;
    return last_beg - 1 - toff - (last - 1 - iexon);
}

static int cmp_indel(const void *a, const void *b)
{
    return ((const indel_t*)a)->off - ((const indel_t*)b)->off;
}

/*
    The transcript offset of the stop codon a frameshift creates, -1 if unknown. The protein predicted
    by csq ends at the stop, e.g. "60DMGQK*>60DWDSLKG*" or "60DMGQK*>60DW..67" with --brief-predictions.
    The translation goes no further than the end of the CDS, so a protein without the stop means the
    frame reads through into the 3' UTR. The stop is located in the CDS of the haplotype, then moved
    back by the indels upstream of it; a haplotype-aware compound can list indels past the stop too.
*/
static int frameshift_stop(gf_tscript_t *tr, tscript_t *tx, bcsq_t *bcsq)
{
    char *str = bcsq->aa ? strchr(bcsq->aa,'>') : NULL, *end;
    if ( !str ) return -1;
    int stop = strtol(str+1, &end, 10);
    if ( end==str+1 ) return -1;
    if ( (str = strstr(end,"..")) ) stop = atoi(str+2);
    else stop += strcspn(end,"*");

    // dna_change, e.g. 1480CA>C+1501G>T
    args.nindel = 0;
    str = bcsq->dna;
    while ( str && *str )
    {
        uint32_t pos = strtoul(str, &end, 10);
        size_t rlen = strcspn(end, ">");
        if ( end==str || end[rlen]!='>' ) break;
        char *alt = end + rlen + 1;
        size_t alen = strcspn(alt, "+");
        int iexon, toff = alen!=rlen ? exonic_offset(tr, tx, pos - 1, &iexon) : -1;
        if ( toff >= 0 )
        {
            hts_expand(indel_t, args.nindel+1, args.mindel, args.indel);
            args.indel[args.nindel].off  = toff - tx->cds_tbeg;
            args.indel[args.nindel].dlen = (int)alen - (int)rlen;
            args.nindel++;
        }
        str = alt + alen;
        if ( *str=='+' ) str++;
    }
    qsort(args.indel, args.nindel, sizeof(*args.indel), cmp_indel);

    int i, cds_off = 3*(stop - 1);
    for (i=0; i<args.nindel && args.indel[i].off < cds_off; i++)
        cds_off -= args.indel[i].dlen;

    int toff = tx->cds_tbeg + cds_off;
    return toff < tx->tlen ? toff : tx->tlen - 1;
}

/*
    Truncations by the LOFTEE-2 distances from the start of the CDS and from the start of the last
    coding exon. A variant starting in an intron is measured from the first base of the next exon.

    What matters for a frameshift is the stop codon it creates: whether it triggers NMD, and whether
    it is close enough to the start for translation to reinitiate downstream. The stop is used when
    csq gives it; END_TRUNC still refers to the frameshift itself and STOP_END_TRUNC to the stop.
    There is no NMD to escape in single-exon transcripts, hence no STOP_END_TRUNC there.
*/
static void check_truncation(gf_tscript_t *tr, tscript_t *tx, variant_t *var, bcsq_t *bcsq)
{
    int iexon, toff = exonic_offset(tr, tx, var->pos, &iexon);
    if ( toff < 0 ) return;

    int stop = toff, istop = iexon;
    if ( bcsq->type & LOF_FRAMESHIFT )
    {
        int fs_stop = frameshift_stop(tr, tx, bcsq);
        if ( fs_stop > toff ) stop = fs_stop, istop = tscript_exon(tx, stop);   // or in the inserted bases
    }
    int from_start = stop - tx->cds_tbeg - (istop - tx->cds_ibeg);

    if ( tx->nexon==1 )
    {
        uint32_t tend = tr->strand==STRAND_FWD ? tr->end : tr->beg;
        uint32_t to_end = var->pos > tend ? var->pos - tend : tend - var->pos;
        if ( to_end <= MIN_DIST_TSCRIPT_END ) add_tag("END_TRUNC", FILTER);
        if ( from_start <= MIN_DIST_CDS_START_1EX ) add_tag("START_TRUNC", FILTER);
        return;
    }

    if ( dist_last_exon(tx, toff, iexon) <= MIN_DIST_LAST_EXON ) add_tag("END_TRUNC", FILTER);
    else if ( dist_last_exon(tx, stop, istop) <= MIN_DIST_LAST_EXON ) add_tag("STOP_END_TRUNC", FILTER);
    if ( from_start <= MIN_DIST_CDS_START ) add_tag("START_TRUNC", FILTER);
}

static void check_phylocsf(gf_tscript_t *tr, int iexon)
{
    const char *id = gff_id2string(args.gff, transcript, tr->id);
    args.tmp.l = 0;
    kputsn(id, strcspn(id,"."), &args.tmp);

    int idx;
    if ( khash_str2int_get(args.phylo_idx, args.tmp.s, &idx) < 0 ) return;
    phylocsf_t *phylo = &args.phylo[idx];
    if ( iexon >= phylo->nexon || isnan(phylo->orf[iexon]) || phylo->orf[iexon] >= 0 ) return;
    add_tag(phylo->max[iexon] > 0 ? "PHYLOCSF_UNLIKELY_ORF" : "PHYLOCSF_WEAK", FLAG);
}

// The 9 bases around an SNV in an acceptor contain another acceptor in frame
static int is_nagnag(gf_tscript_t *tr, variant_t *var)
{
    char seq[10];
    if ( var->pos < 4 || fetch_seq(args.ref, args.chr_ref, var->pos - 4, var->pos + 4, tr->strand, seq) < 0 ) return 0;
    int i;
    for (i=0; i+4<9; i++)
        if ( seq[i]=='A' && seq[i+1]=='G' && seq[i+3]=='A' && seq[i+4]=='G' ) return 1;
    return 0;
}

static void check_intron(gf_tscript_t *tr, tscript_t *tx, variant_t *var, int csq, uint32_t beg, uint32_t end)
{
    int fwd = tr->strand==STRAND_FWD;
    if ( end - beg + 1 < MIN_INTRON_SIZE ) add_tag("SMALL_INTRON", FILTER);
    if ( csq & (LOF_DONOR|LOF_ACCEPTOR) )
    {
        char donor[3], acceptor[3];
        if ( fetch_seq(args.ref, args.chr_ref, fwd ? beg : end - 1, fwd ? beg + 1 : end, tr->strand, donor) < 0 ||
             fetch_seq(args.ref, args.chr_ref, fwd ? end - 1 : beg, fwd ? end : beg + 1, tr->strand, acceptor) < 0 )
            error("Failed to fetch the intron %s:%"PRIu32"-%"PRIu32" from the reference\n", args.chr_ref, beg+1, end+1);

        // The alleles on the transcript strand
        char ref = fwd ? var->ref : complement(var->ref);
        char alt = fwd ? var->alt : complement(var->alt);
        if ( csq & LOF_DONOR && var->ref && !strcmp(donor,"GC") && ref=='C' && alt=='T' ) add_tag("GC_TO_GT_DONOR", FILTER);
        if ( strcmp(donor,"GT") || strcmp(acceptor,"AG") ) add_tag("NON_CAN_SPLICE", FLAG);
        if ( fwd ? var->end < tx->cds_beg : var->beg > tx->cds_end ) add_tag("5UTR_SPLICE", FILTER);
        if ( fwd ? var->beg > tx->cds_end : var->end < tx->cds_beg ) add_tag("3UTR_SPLICE", FILTER);
    }
    if ( csq & LOF_ACCEPTOR && var->ref && is_nagnag(tr, var) ) add_tag("NAGNAG_SITE", FLAG);
}

// The ancestral allele is the one of the purported loss of function. SNVs only, as LOFTEE
static void check_ancestral(variant_t *var)
{
    if ( !args.chr_anc ) { args.nno_ancestral++; return; }
    char base[2];
    if ( fetch_seq(args.anc, args.chr_anc, var->pos, var->pos, STRAND_FWD, base) < 0 ) return;
    if ( base[0]==var->alt ) add_tag("ANC_ALLELE", FILTER);
}

// The filters and flags in the order LOFTEE runs them
static const char *assess(gf_tscript_t *tr, tscript_t *tx, variant_t *var, bcsq_t *bcsq)
{
    int csq = bcsq->type;
    args.csq_filters.l = args.csq_flags.l = 0;

    if ( csq & (LOF_STOP_GAINED|LOF_FRAMESHIFT) )
    {
        int iexon = find_exon(tx, var);
        if ( iexon < 0 ) add_tag("NO_EXON_NUMBER", FLAG);
        else
        {
            check_truncation(tr, tx, var, bcsq);
            if ( tx->nexon==1 ) add_tag("SINGLE_EXON", FLAG);
            else if ( tx->incomplete ) add_tag("INCOMPLETE_CDS", FILTER);
            check_phylocsf(tr, iexon);
        }
    }
    uint32_t beg, end;
    if ( find_intron(tr, tx, var, &beg, &end) >= 0 ) check_intron(tr, tx, var, csq, beg, end);
    if ( var->ref ) check_ancestral(var);

    return args.csq_filters.l ? "LC" : "HC";
}

// The consequence types of a BCSQ entry relevant to LOFTEE, e.g. "frameshift&splice_acceptor|..."
static int parse_csq(const char *str)
{
    int csq = 0;
    while ( *str && *str!='|' )
    {
        size_t len = strcspn(str, "&|");
        if ( len==11 && !strncmp(str,"stop_gained",len) ) csq |= LOF_STOP_GAINED;
        else if ( len==10 && !strncmp(str,"frameshift",len) ) csq |= LOF_FRAMESHIFT;
        else if ( len==12 && !strncmp(str,"splice_donor",len) ) csq |= LOF_DONOR;
        else if ( len==15 && !strncmp(str,"splice_acceptor",len) ) csq |= LOF_ACCEPTOR;
        str += len;
        if ( *str=='&' ) str++;
    }
    return csq;
}

// Split a BCSQ entry in place
static void parse_bcsq(char *str, bcsq_t *bcsq)
{
    char *field[7] = {NULL};
    int i;
    for (i=0; i<7 && str; i++)
    {
        field[i] = str;
        if ( (str = strchr(str,'|')) ) *str++ = 0;
    }
    bcsq->type    = parse_csq(field[0]);
    bcsq->tscript = field[2];
    bcsq->aa      = field[5];
    bcsq->dna     = field[6];
}

/*
    The confidence of one BCSQ entry, NULL if it is not a loss of function. Haplotype-aware csq prints
    a compound consequence at one site only and "@POS" at the others, which are left out, and marks
    with '*' consequences downstream of a stop on the same haplotype, which are assessed as any other
*/
static const char *assess_bcsq(char *str, variant_t *var)
{
    if ( *str=='@' ) return NULL;
    if ( *str=='*' ) str++;

    bcsq_t bcsq;
    parse_bcsq(str, &bcsq);
    if ( !bcsq.type || !bcsq.tscript ) return NULL;

    gf_tscript_t *tr = find_tscript(bcsq.tscript, var);
    if ( !tr ) { args.nunknown_tscript++; return NULL; }
    if ( tr->type!=GF_PROTEIN_CODING ) return NULL;

    tscript_t *tx = tscript_get(tr);
    if ( !tx->ok ) return NULL;
    return assess(tr, tx, var, &bcsq);
}

static int has_lof_csq(const char *bcsq)
{
    return strstr(bcsq,"stop_gained") || strstr(bcsq,"frameshift") || strstr(bcsq,"splice_donor") || strstr(bcsq,"splice_acceptor");
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    char *ref_fname = NULL, *gff_fname = NULL, *anc_fname = NULL, *phylo_fname = NULL;
    int c, force = 0;

    memset(&args, 0, sizeof(args));
    args.hdr_in  = in;
    args.hdr_out = out;
    args.csq_tag = "BCSQ";
    args.rid     = -1;

    static struct option loptions[] =
    {
        {"ancestral",required_argument,NULL,'a'},
        {"custom-tag",required_argument,NULL,'c'},
        {"fasta-ref",required_argument,NULL,'f'},
        {"force",no_argument,NULL,1},
        {"gff-annot",required_argument,NULL,'g'},
        {"phylocsf",required_argument,NULL,'p'},
        {NULL,0,NULL,0}
    };
    while ((c = getopt_long(argc, argv, "a:c:f:g:p:h?", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'a': anc_fname = optarg; break;
            case 'c': args.csq_tag = optarg; break;
            case 'f': ref_fname = optarg; break;
            case 'g': gff_fname = optarg; break;
            case 'p': phylo_fname = optarg; break;
            case  1 : force = 1; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !ref_fname ) error("Missing the -f option\n");
    if ( !gff_fname ) error("Missing the -g option\n");
    if ( !anc_fname ) error("Missing the -a option\n");
    if ( !phylo_fname ) error("Missing the -p option\n");

    if ( !bcf_hdr_idinfo_exists(in, BCF_HL_INFO, bcf_hdr_id2int(in, BCF_DT_ID, args.csq_tag)) )
        error("The INFO/%s tag is not in the header, run bcftools csq first\n", args.csq_tag);
    const char *tags[] = { "LOF", "LOF_FILTERS", "LOF_FLAGS" };
    for (c=0; c<3; c++)
        if ( bcf_hdr_idinfo_exists(in, BCF_HL_INFO, bcf_hdr_id2int(in, BCF_DT_ID, tags[c])) )
            error("The input has INFO/%s already, remove it with \"bcftools annotate -x INFO/LOF,INFO/LOF_FILTERS,INFO/LOF_FLAGS\"\n", tags[c]);

    if ( !(args.ref = fai_load(ref_fname)) ) error("Failed to load the fasta index of %s\n", ref_fname);
    if ( !(args.anc = fai_load(anc_fname)) ) error("Failed to load the fasta index of %s\n", anc_fname);
    phylocsf_init(phylo_fname);

    args.gff = gff_init(gff_fname);
    gff_set(args.gff, force_out_of_phase, force);
    gff_parse(args.gff);
    args.idx_tscript = gff_get(args.gff, idx_tscript);
    args.idx_exon    = gff_get(args.gff, idx_exon);
    args.itr_tscript = regitr_init(args.idx_tscript);
    args.itr_exon    = regitr_init(args.idx_exon);

    bcf_hdr_printf(out, "##INFO=<ID=LOF,Number=.,Type=String,Description=\"LOFTEE confidence of each %s consequence: "
        "HC high, LC low (LOF_FILTERS set), '.' not a loss of function\">", args.csq_tag);
    bcf_hdr_printf(out, "##INFO=<ID=LOF_FILTERS,Number=.,Type=String,Description=\"LOFTEE filters that make each %s consequence LC, "
        "separated by '&', '.' if none\">", args.csq_tag);
    bcf_hdr_printf(out, "##INFO=<ID=LOF_FLAGS,Number=.,Type=String,Description=\"LOFTEE flags of each %s consequence, informative "
        "only, separated by '&', '.' if none\">", args.csq_tag);
    return 0;
}

// Append the value of one consequence to an INFO list
static void push_value(kstring_t *str, int i, const char *value)
{
    if ( i ) kputc(',', str);
    kputs(value ? value : ".", str);
}

bcf1_t *process(bcf1_t *rec)
{
    if ( bcf_get_info_string(args.hdr_in, rec, args.csq_tag, &args.bcsq, &args.mbcsq) <= 0 ) return rec;
    if ( !has_lof_csq(args.bcsq) ) return rec;
    if ( rec->n_allele > 2 ) { args.nmultiallelic++; return rec; }
    if ( rec->n_allele < 2 || rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*' ) return rec;

    set_chr(rec->rid);
    if ( !args.chr_ref ) error("The sequence \"%s\" is not in the reference fasta\n", bcf_seqname(args.hdr_in, rec));

    variant_t var;
    variant_init(&var, rec);

    int nlof = 0, i;
    char *bcsq = args.bcsq;
    args.lof.l = args.filters.l = args.flags.l = 0;
    for (i=0; bcsq; i++)
    {
        char *next = strchr(bcsq, ',');
        if ( next ) *next++ = 0;
        const char *lof = assess_bcsq(bcsq, &var);
        if ( lof ) nlof++;
        push_value(&args.lof, i, lof);
        push_value(&args.filters, i, lof && args.csq_filters.l ? args.csq_filters.s : NULL);
        push_value(&args.flags, i, lof && args.csq_flags.l ? args.csq_flags.s : NULL);
        bcsq = next;
    }
    if ( !nlof ) return rec;

    if ( bcf_update_info_string(args.hdr_out, rec, "LOF", args.lof.s) < 0 ||
         bcf_update_info_string(args.hdr_out, rec, "LOF_FILTERS", args.filters.s) < 0 ||
         bcf_update_info_string(args.hdr_out, rec, "LOF_FLAGS", args.flags.s) < 0 )
        error("Failed to update INFO/LOF at %s:%"PRIhts_pos"\n", bcf_seqname(args.hdr_in, rec), rec->pos + 1);
    return rec;
}

void destroy(void)
{
    if ( args.nmultiallelic )
        fprintf(stderr,"Warning: skipped %d multiallelic records, split them with \"bcftools norm -m-\" before running bcftools csq\n", args.nmultiallelic);
    if ( args.nunknown_tscript )
        fprintf(stderr,"Warning: %d loss-of-function consequences were left out, their transcripts are not in the GFF\n", args.nunknown_tscript);
    if ( args.nno_ancestral )
        fprintf(stderr,"Warning: the ancestral allele of %d consequences was not checked, their chromosome is not in the ancestral sequence\n", args.nno_ancestral);

    tscript_destroy();
    regitr_destroy(args.itr_tscript);
    regitr_destroy(args.itr_exon);
    gff_destroy(args.gff);
    fai_destroy(args.ref);
    fai_destroy(args.anc);
    phylocsf_destroy();
    free(args.bcsq);
    free(args.chr_gff_str.s);
    free(args.chr_ref_str.s);
    free(args.chr_anc_str.s);
    free(args.lof.s);
    free(args.filters.s);
    free(args.flags.s);
    free(args.csq_filters.s);
    free(args.csq_flags.s);
    free(args.tmp.s);
    free(args.indel);
}
