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
log.info "  vf-nf : Nextflow pipeline for variant confidence scoring"
log.info "          based on the random forest algorithm            "
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
    log.info 'nextflow run vf.nf --mode learning --truthVCF GIAB.vcf --newSeqVCF our_GIAB_seq.vcf'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --learning AND/OR --scoring        STRING               Indicate which mode (learning or scoring) to run.'
    log.info ''
    log.info '    if --learning:'
    log.info '      --trainingTable                  TXT                  File containing variant calls used to train the model. Must contain a column "status" if'
    log.info '                                                            options --duplicatedSeqVCF and --duplicatedSeqCov are not provided.'
    log.info '      --features                       LIST                 List of features used to train the model separated by commas, e.g. "--features AF,DP,RVSB".'
    log.info '                                                            Features are predifined when running with option --needlestack.'
    log.info '      if --replication'
    log.info '        --duplicatedSeqVCF             VCF                  Variant calls from an other sequencing, used to assign status to trainingTable.'
    log.info '        --duplicatedSeqCov             TXT                  File containing for each sequenced position in toTrainVCF the coverage in duplicatedSeq data'
    log.info '                                                            for quality check. Use iarcbioinfo/mpileup-nf pipeline on BAM/BED used to generate --toTrainVCF.'
    log.info ''
    log.info '    if --scoring:'
    log.info '      --targetTable                    TXT                  File containing variant calls on which models would be apply.'
    log.info '    if --scoring without --learning:'
    log.info '      --modelSNV                       RDATA                Variant calls from an other sequencing, used to assign status to toTrainTable.'
    log.info '      --modelINDEL                     RDATA                File containing for each sequenced position in toTrainVCF the coverage in duplicatedSeq data'
    log.info ''
    log.info 'Optional arguments:'
    log.info '    --output_folder                    FOLDER               Output folder (default: vf_output).'
    log.info '    --nsplit                           INTEGER              Split the input toTrainVCF to transform it into table in parallel.'
    log.info '    --needlestack                      FLAG                 Specify that calling was launched with needlestack, to use predifined features.'
    log.info ''
    exit 0
}

params.learning = null
params.scoring = null
params.trainingTable = null
params.targetTable = null
params.duplicatedSeqVCF = null
params.duplicatedSeqCov = null
params.modelSNV = null
params.modelINDEL = null
params.features = null
params.output_folder = 'vf_output'
params.nsplit = 1
params.needlestack = null

if(params.learning == null & params.scoring == null){
  exit 1, "Please specify a mode: --learning and/or --scoring"
}

if(params.learning != null){
    if (params.trainingTable == null){
    exit 1, "Please specify a table file to train the model: --trainingTable myfile.txt"
    }
    if (params.features == null & params.needlestack == null){
    exit 1, "Please specify the features if calling was not done with needlestack (--needlestack option)"
    }
    if((params.duplicatedSeqVCF != null && params.duplicatedSeqCov == null) || (params.duplicatedSeqVCF == null && params.duplicatedSeqCov != null)){
    exit 1, "Please provide both --duplicatedSeqVCF and --duplicatedSeqCov"
    }
    if(params.duplicatedSeqVCF == null && params.duplicatedSeqCov == null){
    "Caution: please make sure that --trainingTable file contains a column 'status'"
    }
}

if(params.scoring != null){
    if (params.targetTable == null){
      exit 1, "Please specify a target table file to predict the status: --targetTable myfile.txt"
    }
    if (params.learning == null && (params.modelSNV == null || params.modelINDEL == null)) {
      exit 1, "Please provide trained models when --learning mode is inactivated (e.g. --modelSNV msnv.Rdata --modelINDEL mindel.Rdata)"
    }
}

if(params.trainingTable != null)trainingTable = file(params.trainingTable)
if(params.targetTable != null) targetTable = file(params.targetTable)
if(params.duplicatedSeqVCF != null) duplicatedSeqVCF = file(params.duplicatedSeqVCF)
if(params.duplicatedSeqCov != null) duplicatedSeqCov = file(params.duplicatedSeqCov)
if(params.modelSNV != null) modelSNV = file(params.modelSNV)
if(params.modelINDEL != null) modelINDEL = file(params.modelINDEL)
if(params.needlestack != null) params.features = "RVSB,QVAL,AF,ERR_INFO,DP,medianDP_INFO,FS,MIN_DIST,AO,QUAL,MaxRatioWin,NbVarWin,IoD,HpLength,N_QVAL_20_50_INFO"

if(params.mode == "learning"){

  if(params.duplicatedSeqVCF != null){

      process splitTrainingTable {
        input:
        file trainingTable

        output:
        file 'split*' into splitted_trainingTable mode flatten

        shell:
        '''
        head -n1 !{trainingTable} | sed '/^#CHROM/q' | grep -v "<redacted>" > header
        ((nb_total_lines= $((`zcat !{trainingTable} | wc -l`)) ))
        ((core_lines = $nb_total_lines - $((`cat header | wc -l`)) ))
        ((lines_per_file = ( $core_lines + !{params.nsplit} - 1) / !{params.nsplit}))
        ((start=( $((`cat header | wc -l`)) +1 ) ))
        zcat !{trainingTable} | tail -n+$start | split -l $lines_per_file -a 10 --filter='{ cat header; cat; } | bgzip > $FILE.gz' - split_
        '''
      }

      process addCoverage {
        publishDir params.output_folder+'/COVERAGE/', mode: 'move', pattern: '*.pdf'

        input:
        file duplicatedSeqCov
        file strainingTable from splitted_trainingTable

        output:
        file "*.txt" into trainingTableCov mode flatten
        file  "*.pdf" into PDFoutput

        shell:
        remove_bc_arg = "${params.needlestack}" != null ? "--remove_bc" : ""
        '''
        !{baseDir}/bin/add_coverage_to_annot.r --coverage=!{duplicatedSeqCov} --table=!{strainingTable} !{remove_bc_arg} --plot_coverage
        '''
      }

      if(params.caller == 'needlestack'){
        process splitInfoGeno{
          publishDir params.output_folder+'/INFOGENO/', mode: 'move', pattern: '*.pdf'

          input:
          file table from trainingTableCov

          output:
          file "*.txt" into trainingTableCov_ready

          shell:
          '''
          !{baseDir}/bin/split_INFO_GENOTYPE.r --table=!{table} --split_info --split_geno --plots
          '''
        }
      } else { trainingTableCov.into(trainingTableCov_ready) }

      process addStatus {

      input:
      file table from trainingTableCov_ready
      file duplicatedSeqVCF

      output:
      file "*.txt" into trainingTableCov_ready_status

      shell:
      to_log10_arg = "${params.needlestack}" != null ? "--to_log10=ERR_INFO" : ""
      variable_to_plot_arg = "${params.needlestack}" != null ? "--variable_to_plot=ERR_INFO,RVSB_INFO,AF,FS_INFO,REVEL,PopFreqMax" : ""
      reformat_indels_arg = "${params.needlestack}" != null ? "--reformat_indels" : ""
      '''
      !{baseDir}/bin/add_positive_status.r --table=!{table} --vcf=!{duplicatedSeqVCF} !{to_log10_arg} !{variable_to_plot_arg} !{reformat_indels_arg}
      '''

      }

      if(params.caller == 'needlestack'){
        process addFeatures{

          input:
          file table from trainingTableCov_ready_status

          output:
          file "*.txt" into trainingTableCov_ready_status_features

          shell:
          '''
          !{baseDir}/bin/add_calling_features.r --table=!{table}
          '''
        }
      } else { trainingTableCov_ready_status.into(trainingTableCov_ready_status_features) }

      process mergeTables {

      input:
      file all_table from trainingTableCov_ready_status_features.collect()

      output:
      file "*table.txt" into merged_table

      shell:
      '''
      head -n1 !{all_table[0]} > !{params.trainingTable.baseName}_processed_table.txt
      awk 'FNR>1' !{all_table} >> !{params.trainingTable.baseName}_processed_table.txt
      '''
      }
  }
}
