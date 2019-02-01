#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2018 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  variant-filtering-nf : Nextflow pipeline for variant confidence scoring"
log.info "                         based on deep/machine learning algorithm        "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

params.help = null

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run main.nf --mode learning --truthVCF GIAB.vcf --newSeqVCF our_GIAB_seq.vcf'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --ref                  FILE (with index)     Reference fasta file indexed.'
    log.info '    --mode                 STRING                Indicate which mode (learning or scoring) to run.'
    log.info '    if --mode learning:'
    log.info '      --toTrainVCF         VCF                   Variant calls used to train the model, after assigning a TP/FP status.'
    log.info '      --truthVCF           VCF                   Validated variant calls, used to assign status.'
    log.info '      OR '
    log.info '      --duplicatedSeqVCF   VCF                   Variant calls from an other sequencing, used to assign status.'
    log.info '        AND --duplicatedSeqCov  TXT              File containing for each sequenced position in toTrainVCF the coverage in duplicatedSeq data'
    log.info '                                                 for quality check. Use iarcbioinfo/mpileup-nf pipeline on BAM/BED used to generate --toTrainVCF.'
    log.info '    if --mode scoring:'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --output_folder      FOLDER                  Output folder (default: variant-filtering-output).'
    log.info '    --nsplit             INTEGER                 Split the input toTrainVCF to transform it into table in parallel.'
    log.info ''
    exit 0
}

params.mode = null
params.toTrainVCF = null
params.truthVCF = null
params.duplicatedSeqVCF = null
params.duplicatedSeqCov = null
params.cpu = 1
params.output_folder = 'variant-filtering-output'
params.nsplit = 1
params.caller = 'needlestack'

if(params.mode == null){
  exit 1, "Please specify a mode: --mode [learning,scoring]"
}

if(params.mode == "learning"){
    if (params.toTrainVCF == null){
    exit 1, "Please specify a VCF to train the model: --toTrainVCF myfile.vcf"
    }
    if(params.truthVCF == null && params.duplicatedSeqVCF == null){
    exit 1, "Please specify a VCF to assign TP/FP status: --truthVCF myfile.vcf OR --duplicatedSeqVCF myfile.vcf"
    }
    if(params.truthVCF != null && params.duplicatedSeqVCF != null){
    exit 1, "Please specify either --truthVCF myfile.vcf or --duplicatedSeqVCF, not both"
    }
    if(params.duplicatedSeqVCF != null && params.duplicatedSeqCov == null){
    exit 1, "Please provide --duplicatedSeqCov input using iarcbioinfo/mpileup-nf pipeline on BAM/BED used to generate --toTrainVCF"
    }
}

truthVCF = file(params.truthVCF)
toTrainVCF = file(params.toTrainVCF)
duplicatedSeqVCF = file(params.duplicatedSeqVCF)
duplicatedSeqCov = file(params.duplicatedSeqCov)

if(params.mode == "learning"){

  if(params.truthVCF != null){

  }

  if(params.duplicatedSeqVCF != null){

      process splitVCF {
        input:
        file toTrainVCF

        output:
        file 'split*' into splitted_vcf mode flatten

        shell:
        '''
        zcat !{toTrainVCF} | sed '/^#CHROM/q' | grep -v "<redacted>" > header
        ((nb_total_lines= $((`zcat !{toTrainVCF} | wc -l`)) ))
        ((core_lines = $nb_total_lines - $((`cat header | wc -l`)) ))
        ((lines_per_file = ( $core_lines + !{params.nsplit} - 1) / !{params.nsplit}))
        ((start=( $((`cat header | wc -l`)) +1 ) ))
        zcat !{toTrainVCF} | tail -n+$start | split -l $lines_per_file -a 10 --filter='{ cat header; cat; } | bgzip > $FILE.gz' - split_
        '''
      }

      process addCoverage {
        publishDir params.output_folder+'/COVERAGE/', mode: 'move', pattern: '*.pdf'

        input:
        file duplicatedSeqCov
        file svcf from splitted_vcf

        output:
        file "*addCoverage.vcf" into toTrainVCFCov mode flatten
        file  "*.pdf" into PDFoutput

        shell:
        remove_bc_arg = "${params.caller}" == "needlestack" ? "--remove_bc" : ""
        '''
        !{baseDir}bin/add_coverage.r --coverage=!{duplicatedSeqCov} --VCFToAddCoverage=!{svcf} !{remove_bc_arg}
        '''
      }

      process VCFToTable {
        publishDir params.output_folder+'/INFOGENO/', mode: 'move', pattern: '*.pdf'

        input:
        file toTrainVCFCovFile from toTrainVCFCov

        output:
        file "*.txt" into toTrainTable
        file  "*.pdf" into PDFoutput

        shell:
        '''
        !{baseDir}bin/VCF_to_table.r --vcf=!{toTrainVCFCovFile}
        '''
      }

if(params.caller == 'needlestack'){
  process AddTableFeatures {
    publishDir params.output_folder+'/INFOGENO/', mode: 'move', pattern: '*.pdf'

    input:
    file table from toTrainTable

    output:
    file "*.txt" into toTrainTableFeatures

    shell:
    '''
    !{baseDir}bin/add_calling_features.r --table=!{table}
    '''
  }
} else { toTrainTable.into(toTrainTableFeatures) }

      process mergeTables {

        input:
        file all_table from toTrainTableFeatures.collect()

        output:
        file "*table.txt" into merged_table

        shell:
        '''
        head -n1 !{all_table[0]} > !{params.toTrainVCF.baseName}_table.txt
        awk 'FNR>1' !{all_table} >> !{params.toTrainVCF.baseName}_table.txt
        '''
      }

    }
  }
}
