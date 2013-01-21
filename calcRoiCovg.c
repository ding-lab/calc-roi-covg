/// Author: Cyriac Kandoth
/// Date: 01.19.2011
/// Description: Counts bases with sufficient read-depth in regions of interest within two BAMs
/// Notes:
/// - If ROIs of the same gene overlap, they will not be merged. Use BEDtools' mergeBed if needed
/// - The totals written at the end count each base only once, even if it is in multiple ROIs

#define _GNU_SOURCE
#include <stdio.h>
#include <stdbool.h>
#include "sam.h"
#include "faidx.h"
#include "khash.h"

KHASH_MAP_INIT_STR( s, int )

// Initializes the header hash in a given bam_header_t. Defined in samtools/bam_aux.c
void bam_init_header_hash( bam_header_t *header );

// Create some values that represent the different refseq bp-classes
enum bp_class_t { AT, CG, CpG, IUB, UNKNOWN };

typedef struct
{
  int beg, end; // The start and stop of a region of interest
  int min_mapq; // Minimum mapping quality of the reads to pileup
  int min_depth_bam1, min_depth_bam2; // Minimum read depth required per bam
  bool *bam1_cvg; // Tags bases in a region of bam1 with the minimum required read-depth
  int covd_bases; // Counts bases in a region that has the minimum read depth in both bams
  int base_cnt[4]; // Counts covered bases in an ROI of 4 bp-classes AT, CG, CpG, IUB
  int tot_covd_bases; // Counts bases in all ROIs that have the minimum read depth in both bams
  int tot_base_cnt[4]; // Counts covered bases in all ROIs of 4 bp-classes AT, CG, CpG, IUB
  char *ref_seq; // Contains the reference sequence for the entire chromosome a region lies in
  char *bp_class; // This prevents counting the same base twice when ROIs overlap in a chromosome
  int ref_id, ref_len; // A chromosome's ID in the the BAM header hash, and it's length
  samfile_t *sam1, *sam2; // The two bam files that need to be piled-up
} pileup_data_t;

// Callback for bam_fetch() pushed only alignments that pass the minimum mapping quality
static int fetch_func( const bam1_t *b, void *data )
{
  bam_plbuf_t *buf = (bam_plbuf_t*)data;
  bam_plbuf_push( b, buf ); // Invokes the callback function specified in bam_plbuf_init()
  return 0;
}

// Callback for bam_plbuf_init() when running a pileup on bam1
static int pileup_func_1( uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data )
{
  pileup_data_t *tmp = (pileup_data_t*)data;
  if( (int)pos >= tmp->beg && (int)pos < tmp->end )
  {
    // Count the number of reads that pass the mapping quality threshold across this base
    int i, mapq_n = 0;
    for( i = 0; i < n; ++i )
    {
      const bam_pileup1_t *base = pl + i;
      if( !base->is_del && base->b->core.qual >= tmp->min_mapq )
        mapq_n++;
    }
    tmp->bam1_cvg[(int)pos - tmp->beg] = ( mapq_n >= tmp->min_depth_bam1 );
  }
  return 0;
}

// Callback for bam_plbuf_init() when running a pileup on bam2
static int pileup_func_2( uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data )
{
  pileup_data_t *tmp = (pileup_data_t*)data;
  if( (int)pos >= tmp->beg && (int)pos < tmp->end && tmp->bam1_cvg[(int)pos - tmp->beg] )
  {
    // Count the number of reads that pass the mapping quality threshold across this base
    int i, mapq_n = 0;
    for( i = 0; i < n; ++i )
    {
      const bam_pileup1_t *base = pl + i;
      if( !base->is_del && base->b->core.qual >= tmp->min_mapq )
        mapq_n++;
    }
    if( mapq_n >= tmp->min_depth_bam2 )
    {
      int class = (int)(tmp->bp_class[pos]);
      if( class == UNKNOWN )
      {
        char base = tmp->ref_seq[pos];
        char prev_base = tmp->ref_seq[pos-1];
        char next_base = tmp->ref_seq[pos+1];
        if( base == 'a' || base == 'A' || base == 't' || base == 'T' )
          class = AT;
        else if((( base == 'c' || base == 'C' ) && ( next_base == 'g' || next_base == 'G' )) ||
                (( base == 'g' || base == 'G' ) && ( prev_base == 'c' || prev_base == 'C' )))
          class = CpG;
        else if( base == 'c' || base == 'C' || base == 'g' || base == 'G' )
          class = CG;
        else
          class = IUB;
        ++tmp->covd_bases;
        ++tmp->base_cnt[class];
        ++tmp->tot_covd_bases;
        ++tmp->tot_base_cnt[class];
        tmp->bp_class[pos] = (char)class; // Tag this as seen and save its class for an overlapping ROI
      }
      else
      {
        ++tmp->covd_bases;
        ++tmp->base_cnt[class];
      }
    }
  }
  return 0;
}

int main( int argc, char *argv[] )
{
  // This data is shared across functions
  pileup_data_t data;
  data.ref_id = -1;
  data.ref_seq = NULL;
  data.bp_class = NULL;

  // Set default min_mapq and min_depths
  data.min_mapq = 20;
  data.min_depth_bam1 = 6;
  data.min_depth_bam2 = 8;
  if( argc != 6 && argc != 9 )
  {
    fprintf( stderr, "\nUsage: %s <bam1> <bam2> <roi_file> <ref_seq_fasta> <output_file> [min_depth_bam1 min_depth_bam2 min_mapq]", argv[0] );
    fprintf( stderr, "\nDefaults: min_depth_bam1 = %i, min_depth_bam2 = %i, min_mapq = %i", data.min_depth_bam1, data.min_depth_bam2, data.min_mapq );
    fprintf( stderr, "\nNOTE: ROI file *must* be sorted by chromosome/contig names\n\n" );
    return 1;
  }
  // Set user-defined min_depths if specified
  if( argc == 9 )
  {
    data.min_depth_bam1 = atoi( argv[6] );
    data.min_depth_bam2 = atoi( argv[7] );
    data.min_mapq = atoi( argv[8] );
  }

  // Open both BAM files and load their index files
  data.sam1 = samopen( argv[1], "rb", 0 );
  if( !data.sam1 )
    fprintf( stderr, "Failed to open BAM file %s\n", argv[1] );
  bam_index_t *idx1 = bam_index_load( argv[1] );
  if( !idx1 )
    fprintf( stderr, "BAM index file is not available for %s\n", argv[1] );
  data.sam2 = samopen( argv[2], "rb", 0 );
  if( !data.sam2 )
    fprintf( stderr, "Failed to open BAM file %s\n", argv[2] );
  bam_index_t *idx2 = bam_index_load( argv[2] );
  if( !idx2 )
    fprintf( stderr, "BAM index file is not available for %s\n", argv[2] );

  // Open the file with the annotated regions of interest
  FILE *roiFp = fopen( argv[3], "r" );
  if( !roiFp )
    fprintf( stderr, "Failed to open ROI file %s\n", argv[3] );

  // Load an index to the reference sequence fasta file
  faidx_t *ref_fai = fai_load( argv[4] );
  if( !ref_fai )
    fprintf( stderr, "Failed to open reference fasta file %s\n", argv[4] );

  // Open the output file to write to
  FILE* outFp = fopen( argv[5], "w" );
  if( !outFp )
    fprintf( stderr, "Failed to open output file %s\n", argv[5] );

  // Show the user any and all errors they need to fix above before quitting the program
  if( !data.sam1 || !idx1 || !data.sam2 || !idx2 || !roiFp || !ref_fai || !outFp )
    return 1;

  // Write a header with column titles for the output file
  fprintf( outFp, "#NOTE: Last line in file shows non-overlapping totals across all ROIs\n" );
  fprintf( outFp, "#Gene\tROI\tLength\tCovered\tATs_Covered\tCGs_Covered\tCpGs_Covered\n" );

  // Initialize a header hash to check for valid ref_names in bam1
  bam_init_header_hash( data.sam1->header );
  khash_t(s) *hdr_hash = (khash_t(s)*)data.sam1->header->hash;

  // Initialize the counters for the total number of non-overlapping bases in all ROIs
  data.tot_covd_bases = data.tot_base_cnt[AT] = data.tot_base_cnt[CG] = data.tot_base_cnt[CpG] = 0;

  size_t length;
  char *line = NULL;
  line = (char*)malloc( 200 );
  while( getline( &line, &length, roiFp ) != -1 )
  {
    char ref_name[50], gene_name[100];
    if( sscanf( line, "%s %i %i %s", ref_name, &data.beg, &data.end, gene_name ) == 4 )
    {
      int ref_id;
      // If this region is valid in bam1, we'll assume it's also valid in bam2
      khiter_t iter = kh_get( s, hdr_hash, ref_name );
      if( iter == kh_end( hdr_hash ) || data.beg > data.end )
      {
        fprintf( stderr, "Skipping invalid ROI: %s", line );
      }
      else
      {
        --data.beg; // Make the start locus a 0-based coordinate
        ref_id = kh_value( hdr_hash, iter );
        int bases = data.end - data.beg;
        data.covd_bases = data.base_cnt[AT] = data.base_cnt[CG] = data.base_cnt[CpG] = 0;
        data.bam1_cvg = (bool*)calloc( bases, sizeof( bool )); // calloc also sets them to zero

        // Load this whole chromosome's refseq unless already loaded for the previous ROI
        if( data.ref_seq == NULL || ref_id != data.ref_id )
        {
          if( data.ref_seq )
            free( data.ref_seq );
          if( data.bp_class )
            free( data.bp_class );
          data.ref_seq = fai_fetch( ref_fai, data.sam1->header->target_name[ref_id], &data.ref_len );
          data.bp_class = (char*)malloc( data.ref_len * sizeof( char ));
          memset( data.bp_class, UNKNOWN, data.ref_len );
          data.ref_id = ref_id;
        }

        // If the ROI is at a chromosome tip, edit it so we can look for CpGs without segfaulting
        if( data.beg == 0 )
          ++data.beg;
        if( data.end == data.ref_len )
          --data.end;

        // Pileup bam1 and tag all the bases which have sufficient read depth
        bam_plbuf_t *buf1 = bam_plbuf_init( pileup_func_1, &data ); // Initialize pileup
        bam_fetch( data.sam1->x.bam, idx1, ref_id, data.beg, data.end, buf1, fetch_func );
        bam_plbuf_push( 0, buf1 );
        bam_plbuf_destroy( buf1 );

        // Pileup bam2 and count bases with sufficient read depth, and tagged earlier in bam1
        bam_plbuf_t *buf2 = bam_plbuf_init( pileup_func_2, &data ); // Initialize pileup
        bam_fetch( data.sam2->x.bam, idx2, ref_id, data.beg, data.end, buf2, fetch_func );
        bam_plbuf_push( 0, buf2 );
        bam_plbuf_destroy( buf2 );

        fprintf( outFp, "%s\t%s:%i-%i\t%i\t%i\t%i\t%i\t%i\n", gene_name, ref_name, data.beg+1, data.end,
                 bases, data.covd_bases, data.base_cnt[AT], data.base_cnt[CG], data.base_cnt[CpG] );
        free( data.bam1_cvg );
      }
    }
    else
    {
      fprintf( stderr, "Badly formatted ROI: %s", line );
      fprintf( stderr, "\nROI file should be a tab-delimited list of [chrom, start, stop, annotation]" );
      fprintf( stderr, "\nwhere start and stop are both 1-based chromosomal loci" );
      fprintf( stderr, "\nFor example:\n20\t44429404\t44429608\tELMO2\nMT\t5903\t7445\tMT-CO1\n" );
      fprintf( stderr, "\nNOTE: ROI file *must* be sorted by chromosome/contig names\n\n" );
      return 1;
    }
  }

  // The final line in the file contains the non-overlapping base counts across all ROIs
  fprintf( outFp, "#NonOverlappingTotals\t\t\t%i\t%i\t%i\t%i\n", data.tot_covd_bases,
           data.tot_base_cnt[AT], data.tot_base_cnt[CG], data.tot_base_cnt[CpG] );

  // Cleanup
  if( line )
    free( line );
  if( data.ref_seq )
    free( data.ref_seq );
  if( data.bp_class )
    free( data.bp_class );
  bam_index_destroy( idx1 );
  bam_index_destroy( idx2 );
  samclose( data.sam1 );
  samclose( data.sam2 );
  fai_destroy( ref_fai );
  fclose( roiFp );
  fclose( outFp );
  return 0;
}
