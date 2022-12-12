// We can't change global parameters inside this scope, so we build the ones we need locally
def n_haps = 0
if (!params.smoothxg_haplotypes_smooth) {
  n_haps = params.n_haplotypes
}

def wfmash_merge_cmd = params.wfmash_merge_segments ? "-M" : ""
def wfmash_exclude_cmd = params.wfmash_exclude_delim ? "-Y ${params.wfmash_exclude_delim}" : "-X"
def wfmash_split_cmd = params.wfmash_no_splits ? "-N" : ""
def wfmash_block_length = params.wfmash_segment_length*5
def wfmash_block_length_cmd = "-l ${wfmash_block_length}"
def wfmash_mash_kmer_cmd = "-k ${params.wfmash_mash_kmer}"
def wfmash_kmer_thres_cmd = "-H ${params.wfmash_mash_kmer_thres}"
def wfmash_n_mappings_minus_1 = params.n_haplotypes - 1
def wfmash_sparse_map_cmd = ""
if (params.wfmash_sparse_map == "auto") {
  n = n_haps
  x = Math.log(n)/n * 10
  wfmash_sparse_map_frac = 1
  if (x >= 1) {
    wfmash_sparse_map_frac = x
  }
  wfmash_sparse_map_cmd = "-x${wfmash_sparse_map_frac}"
} else {
  if (params.wfmash_sparse_map != null) {
    wfmash_sparse_map_cmd = "-x${params.wfmash_sparse_map}"
  }
}
def wfmash_temp_dir = params.wfmash_temp_dir ? "-B${params.wfmash_temp_dir}" : ""

def wfmash_prefix = "wfmash"

def make_file_prefix = { f -> """\
${f.getName()}\
""" }

workflow COMMUNITY {
  take:
  ch_fasta

  main:

  fasta = ch_fasta.map { f -> tuple(make_file_prefix(f), f) }
  fasta_file_name = ch_fasta.map {it.getName()}

  ch_empty = Channel.empty()
  // we currently don't want to emit anything

  emit:
  ch_fasta

}