nextflow.enable.dsl=2

// ---- require runs ----
if( !params.runs ) error "config.yaml must define a 'runs' list (sample_id, type, fastq1, fastq2)."

// ---- build simple lists (no optional syntax) ----
def PE_RNA_LIST = []
def SE_RNA_LIST = []
def PE_ALL_LIST = []
def SE_ALL_LIST = []

params.runs.each { m ->
  def sid  = (m.sample_id ?: new File(m.fastq1).getName())
  def typ  = m.type?.toString()?.toLowerCase()
  def r1   = file(m.fastq1)
  def isSE = (!m.fastq2) || (m.fastq2 == 'NA')

  if( isSE ) {
    // SE tuples: (r1, dummy=null, type, sample)
    if( typ == 'rna' ) SE_RNA_LIST << [ r1, null, typ, sid ]
    SE_ALL_LIST << [ r1, null, typ, sid ]
  } else {
    def r2 = file(m.fastq2)
    // PE tuples: (r1, r2, type, sample)
    if( typ == 'rna' ) PE_RNA_LIST << [ r1, r2, typ, sid ]
    PE_ALL_LIST << [ r1, r2, typ, sid ]
  }
}

// ---- channels ----
def PE_RNA_CH = Channel.from(PE_RNA_LIST)
def SE_RNA_CH = Channel.from(SE_RNA_LIST)
def PE_ALL_CH = Channel.from(PE_ALL_LIST)
def SE_ALL_CH = Channel.from(SE_ALL_LIST)

// constant value channel for the PHLAT index (reused for every run)
def PHLAT_INDEX_VAL = Channel.value( file(params.phlat_index) )

// ======================= arcasHLA (RNA only) =======================

// PE RNA
process ARCASHLA_PE {
  publishDir { "${params.outdir}/arcashla/${sample}" }, mode: 'copy'
  cpus   { params.cpus as int }
  memory { "${params.memory_gb as int} GB" }
  container 'docker://jfx319/arcashla'

  input:
  tuple path(r1), path(r2), val(type), val(sample)

  output:
  path "arcashla_out"

  """
  set -euo pipefail
  mkdir -p arcashla_out
  echo "arcasHLA PE | sample=${sample}"
  arcasHLA genotype \
    -t ${params.cpus} \
    -g "${params.arcashla_gene_list}" \
    "${r1}" "${r2}" \
    -o arcashla_out
  """
}

// SE RNA
process ARCASHLA_SE {
  publishDir { "${params.outdir}/arcashla/${sample}" }, mode: 'copy'
  cpus   { params.cpus as int }
  memory { "${params.memory_gb as int} GB" }
  container 'docker://jfx319/arcashla'

  input:
  tuple path(r1), val(dummy), val(type), val(sample)

  output:
  path "arcashla_out"

  """
  set -euo pipefail
  mkdir -p arcashla_out
  echo "arcasHLA SE | sample=${sample}"
  arcasHLA genotype \
    -t ${params.cpus} \
    -g "${params.arcashla_gene_list}" \
    "${r1}" \
    -o arcashla_out
  """
}

// ======================= OptiType (RNA or DNA) =======================

process OPTITYPE_PE {
  publishDir { "${params.outdir}/optitype/${sample}" }, mode: 'copy'
  cpus   { params.cpus as int }
  memory { "${params.memory_gb as int} GB" }
  container 'docker://fred2/optitype'

  input:
  tuple path(r1), path(r2), val(type), val(sample)

  output:
  path "optitype_out"

  """
  set -euo pipefail
  mkdir -p optitype_out
  echo "OptiType PE | sample=${sample} | type=${type}"
  OptiTypePipeline.py \
    -i "${r1}" "${r2}" \
    --${type} \
    -o optitype_out \
    -p ${sample} \
    -v
  """
}

process OPTITYPE_SE {
  publishDir { "${params.outdir}/optitype/${sample}" }, mode: 'copy'
  cpus   { params.cpus as int }
  memory { "${params.memory_gb as int} GB" }
  container 'docker://fred2/optitype'

  input:
  tuple path(r1), val(dummy), val(type), val(sample)

  output:
  path "optitype_out"

  """
  set -euo pipefail
  mkdir -p optitype_out
  echo "OptiType SE | sample=${sample} | type=${type}"
  OptiTypePipeline.py \
    -i "${r1}" \
    --${type} \
    -o optitype_out \
    -p ${sample} \
    -v
  """
}

// ======================= PHLAT (RNA or DNA) =======================

process PHLAT_PE {
  publishDir { "${params.outdir}/phlat/${sample}" }, mode: 'copy'
  cpus   { params.cpus as int }
  memory { "${params.memory_gb as int} GB" }
  container 'docker://stevenjpastor/phlat_env'

  input:
  tuple path(r1), path(r2), val(type), val(sample)
  path index_dir

  output:
  path "phlat_out"

  """
  set -euo pipefail
  mkdir -p phlat_out

  # per-task writable tmp to avoid /tmp FIFO errors
  export TMPDIR="\$PWD/tmp_phlat"; export TMP="\$TMPDIR"; export TEMP="\$TMPDIR"
  mkdir -p "\$TMPDIR"

  echo "PHLAT PE | sample=${sample} | type=${type}"
  python2 -O /usr/bin/phlat-release/dist/PHLAT.py \
    -1 "${r1}" \
    -2 "${r2}" \
    -index "${index_dir}" \
    -b2url /usr/bin/bowtie2 \
    -tag ${sample} \
    -e /usr/bin/phlat-release \
    -o phlat_out \
    -p ${params.cpus}
  """
}

process PHLAT_SE {
  publishDir { "${params.outdir}/phlat/${sample}" }, mode: 'copy'
  cpus   { params.cpus as int }
  memory { "${params.memory_gb as int} GB" }
  container 'docker://stevenjpastor/phlat_env'

  input:
  tuple path(r1), val(dummy), val(type), val(sample)
  path index_dir

  output:
  path "phlat_out"

  """
  set -euo pipefail
  mkdir -p phlat_out

  # per-task writable tmp to avoid /tmp FIFO errors
  export TMPDIR="\$PWD/tmp_phlat"; export TMP="\$TMPDIR"; export TEMP="\$TMPDIR"
  mkdir -p "\$TMPDIR"

  echo "PHLAT SE | sample=${sample} | type=${type}"
  python2 -O /usr/bin/phlat-release/dist/PHLAT.py \
    -1 "${r1}" \
    -index "${index_dir}" \
    -b2url /usr/bin/bowtie2 \
    -tag ${sample} \
    -e /usr/bin/phlat-release \
    -o phlat_out \
    -p ${params.cpus}
  """
}

// ======================= workflow =======================
workflow {
  // arcasHLA only for RNA
  ARCASHLA_PE(PE_RNA_CH)
  ARCASHLA_SE(SE_RNA_CH)

  // OptiType for all
  OPTITYPE_PE(PE_ALL_CH)
  OPTITYPE_SE(SE_ALL_CH)

  // PHLAT for all (broadcast index to each run)
  PHLAT_PE(PE_ALL_CH, PHLAT_INDEX_VAL)
  PHLAT_SE(SE_ALL_CH, PHLAT_INDEX_VAL)
}

