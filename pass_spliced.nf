#!/usr/bin/env nextflow

/*
*Bug:
*
*/



/*
*params input 
*/
params.csfasta_1 = "$baseDir/data/data_pass/*F3.csfasta"
params.csfasta_2 = "$baseDir/data/data_pass/*F5-BC.csfasta"
params.qual_1 = "$baseDir/data/data_pass/*F3.QV.qual"
params.qual_2 = "$baseDir/data/data_pass/*F5-BC.QV.qual"
params.genome = "$baseDir/data/data_pass/1M_hg19.fasta"
params.path_Pass  = "/usr/local/bin/pass" 
params.cpu = "4" 
params.index = null
params.help = false

//print usage
if (params.help) {
    log.info ''
    log.info 'Pass pipeline RNA SOLID'
    log.info '-----------------------'
    log.info '.'
    log.info ''
    log.info 'Usage: '
    log.info '    pass.nf --csfasta_1  .../*1.{csfasta,QV.qual} --csfasta_2 .../2.{csfasta,QV.qual} '
    log.info '            --genome GENOME_FILE '
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --csfasta_1                         First read cs.fasta [ex :".../*1.csfasta."'
    log.info '    --csfasta_2                         Second read cs.fasta [ex :".../*2.csfasta"'
    log.info '    --qual_1                            First File QV.qual [ex :.../*1.qual'
    log.info '    --qual_2                            Second File QV.qual [ex :.../*2.qual'
    log.info '    --index GENOME_INDEX_FILE           Index file.[Optional but must faster].'
    log.info '    --genome GENOME_FILE                Reference genome file(s).'
    log.info '    --path_Pass                         path of tool Pass (by default /usr/local/bin/pass).'  
    log.info '    --cpu                               Numbers of core use (CPU)by default [4].'  
    exit 1
}

/*
*create a read_pairs by params.reads and genome ref 
*/
genome_file = file(params.genome)

read_1 = Channel.fromPath(params.csfasta_1).map { file -> tuple(file.baseName[0..-3], file) }
read_2 = Channel.fromPath(params.csfasta_2).map { file -> tuple(file.baseName[0..-6], file) }
qual_1 = Channel.fromPath(params.qual_1).map { file -> tuple(file.baseName[0..-6], file) }
qual_2 = Channel.fromPath(params.qual_2).map { file -> tuple(file.baseName[0..-9], file) }

csfasta_qual_1 = read_1.combine(qual_1,by:0)
csfasta_qual_2 = read_2.combine(qual_2,by:0)

path_Pass= params.path_Pass
cpu = params.cpu


//csfasta_qual_1.subscribe { println "$it"}
//csfasta_qual_2.subscribe { println "$it"}

if(params.index == null){


    process buildIndex{

        input:
        file genome from genome_file
    
        output: 
        file '*.bin' into genome_index

        """
        ${path_Pass} -p 11111100111111 -d  $genome \
             -cpu ${cpu} -solidCS -D ${genome}.bin
        """
    }
}
else{
    genome_index = file(params.index)
}


/*
* Process for mapping read 1 
*
* Strategy to optimize spliced mapping that requires 3 steps :
*   -process passRead_1_step_1 : Produce global alignments and recover the not aligned
*   -process passRead_1_step_2 : Mapping coordinates  the "not aligned reads" using the local alignment mapping. 
*   -process passRead_1_step_3 : Focusing the spliced alignments at the found coordinates
*OPTION:
*-flc (INT) Low complexity filter for seed words indexing :
*           0: no seed word filter; 
*           1: it filters the seed word if there are at least two homopolymeric hexamers;
*           >=n: it filters the seed word if there is a number of same hexamers;
*-fid       Filter alignments that have the percent identities less than this value [90];
*-sam       Format output;
*-b         Best hits
*-not_aligned       It appends to not_aligned.fa file unaligned sequences in the same format as the input.
*-pst_word_range    The PST word range determines number and type of PST processed during the initial 
                    elaboration of the structures
*/

process spliced_Alignment_Read_1_step_1{
    tag{id}

    input: 
    set id ,file(csfasta), file(qual) from csfasta_qual_1
    file genome from genome_index

    output: 
    set id , file ("*global_result.sam") into global_result_1
    set id , file ("*not_aligned.csfasta"), file ("*not_aligned.qual") into not_aligned_1

    
    """
${path_Pass} -R $genome  \
    -cpu ${cpu} -flc 1 -fid 90 -sam -b -pst_word_range 4 4 \
    -csfasta $csfasta \
    -qual $qual \
    -not_aligned -na_file ${id}1_not_aligned \
    > ${id}1_global_result.sam
    """
}

/*
*OPTION:
*-l     Enable the local alignment function. It try to recognize local small alignments of each read;
*-fle   length of the minimal alignment size.
*-focus_check   Save candidates spliced reads Into 'filename' in order to map them again, 
*               using short seed words together the parameter -focus.
*/

process spliced_Alignment_Read_1_step_2{
    tag{id}

    input: 
    set id ,file(csfasta), file(qual) from not_aligned_1
    file genome from genome_index

    output: 
     set id , file ("*_focus.csfasta"), file ("*_focus.qual") into focus_1

    """
${path_Pass} -R $genome  \
    -cpu ${cpu} -flc 1 -fid 90 -sam\
    -csfasta $csfasta \
    -qual $qual \
    -l -fle 10 -focus_check ${id}1_focus \
    > /dev/null 
    """
}


/*
*OPTION:
*-l     Enable the local alignment function. It try to recognize local small alignments of each read.
*-fle   length of the minimal alignment size.
*-percent_tolerance     The tolerated percent of read size not covered by the alignment.
*-spliced   It Outputs spliced alignments using the user intron scoring system.
*/

process spliced_Alignment_Read_1_step_3{
    tag{id}

    input: 
    set id ,file(csfasta), file(qual) from focus_1
    file genome from genome_index

    output: 
    set id , file ("*_spliced_result.sam") into spliced_result_1


    """
${path_Pass} -R $genome  \
    -cpu ${cpu} -flc 1 -fid 90 -sam \
    -csfasta $csfasta \
    -qual $qual \
    -spliced dna -percent_tolerance 30 -fle 10 \
    -not_aligned \
    > ${id}1_spliced_result.sam
    """
}


/*
* Process for mapping read 2 
*
*   strategy to optimize spliced mapping that requires 3 steps :
*   -process passRead_1_step_1 : Produce global alignments and recover the not aligned
*   -process passRead_1_step_2 : Mapping coordinates  the "not aligned reads" using the local alignment mapping. 
*   -process passRead_1_step_3 : Focusing the spliced alignments at the found coordinates
*OPTION: 
* Same the process for read 1
*/

process spliced_Alignment_Read_2_step_1{
    tag{id}

    input: 
    set id , file(csfasta), file(qual) from csfasta_qual_2
    file genome from genome_index

    output: 
    set id ,file ("*_global_result.sam") into global_result_2
    set id , file ("*not_aligned.csfasta"), file ("*not_aligned.qual") into not_aligned_2


    """
${path_Pass} -R $genome  \
    -cpu ${cpu} -flc 1 -fid 90 -sam -b -pst_word_range 4 4 \
    -csfasta $csfasta \
    -qual $qual \
    -not_aligned -na_file $id"2_not_aligned" \
    > ${id}2_global_result.sam
    """
}

process spliced_Alignment_Read_2_step_2{
    tag{id}

    input: 
    set id ,file(csfasta), file(qual) from not_aligned_2
    file genome from genome_index

    output: 
   set id , file ("*_focus.csfasta"), file ("*_focus.qual") into focus_2


    """
${path_Pass} -R $genome  \
    -cpu ${cpu} -flc 1 -fid 90 -sam \
    -csfasta $csfasta \
    -qual $qual \
    -l -fle 10 -focus_check $id"2_focus" \
    > /dev/null 
    """
}

process spliced_Alignment_Read_2_step_3{
    tag{id}

    input: 
    set id ,file(csfasta), file(qual) from focus_2
    file genome from genome_index

    output: 
    set id, file ("*spliced_result.sam") into spliced_result_2


    """
${path_Pass} -R $genome  \
    -cpu ${cpu} -flc 1 -fid 90 -sam \
    -csfasta $csfasta \
    -qual $qual \
    -spliced dna -percent_tolerance 30 -fle 10 \
    -not_aligned \
    > ${id}2_spliced_result.sam
    """
}


global = global_result_1.combine (global_result_2,by:0)
spliced = spliced_result_1.combine (spliced_result_2,by:0)

/*
*This program executes a refined paired-end with two files SAM  
*
*OPTION :
* -range (int , int, int)  (Estimated_Min_Distance) (Estimated_Max_Distance) (Tolerated_Max_distance)  [0 4000 100000]
* -stdout  Instead to generate a file for each class of paired-ends it outputs all information to the stdout
* -pe_type (int)     if 0    [ ---> <--- ] : typically for BAC-END or FOSMID-END or ILLUMINA paired-end library
*                    if 1    [ ---> ---> or <--- <--- ] : typically for SOLiD mate-pair [default] 
*                    if 2    [ <--- ---> ] : typically for ILLUMINA paired-end library
* -tags (string , string) tag for each paired (tag1_end tag2_end) (example: R3 F3 default)
* -unique_pair (int)     if 1 save unique pair
* -unique_single (int)   if 1 save unique single
*/
process pairing_global_result{
    tag{id}
    errorStrategy 'ignore'

    input: 
    set id , file (read_1_map) , file (read_2_map) from global

    output: 
    set id , file ("*_global_paired_alignments.sam") into global_result_paired


    """
${path_Pass} -program pairing \
    -sam1 $read_1_map \
    -sam2 $read_2_map \
    -range 0 3000 10000 \
    -unique_pair 1 \
    -unique_single 1 \
    -stdout -cpu ${cpu} -pe_type 0 -tags F3 F5-BC \
    > ${id}_global_paired_alignments.sam


    """
    //attention au tags -tags F3 F5-BC
}


process pairing_spliced_result{
    tag{id}
    errorStrategy 'ignore'    

    input: 
    set id , file (read_1_map), file (read_2_map) from spliced


    output: 
    set id , file ("*_spliced_paired_alignments.sam") into spliced_result_paired


    """
${path_Pass} -program pairing \
    -sam1 $read_1_map \
    -sam2 $read_2_map \
    -range 0 3000 10000 \
    -unique_pair 1 \
    -unique_single 1 \
    -stdout -cpu ${cpu} -pe_type 0 -tags F3 F5-BC \
    > ${id}_spliced_paired_alignments.sam

    """

        //attention au tags -tags F3 F5-BC
}

/*
*Convert file SAM to BAM with samtools
*
*OPTION :
*-S       input is SAM
*-b       output BAM
*-o FILE  output file name [stdout]
*/

process convert_Sam_to_Bam_spliced{
    tag{id}
    errorStrategy 'ignore'
    
    publishDir "result/pass/${id}"
    input:
    set id , file (spliced_paired_alignement) from spliced_result_paired

    output:
    set id ,file ("${id}pass_spliced/*_spliced_paired_alignments.bam") into spliced_paired_bam


    """
    samtools view -Sb -o ${id}_spliced_paired_alignments.bam $spliced_paired_alignement
    samtools flagstat ${id}_spliced_paired_alignments.bam > stat

    mkdir ${id}pass_spliced
    mv ${id}*spliced_alignments* ${id}pass_spliced/.
    mv stat ${id}pass_spliced/.
    """

}
process convert_Sam_to_Bam_global{
    tag{id}
    errorStrategy 'ignore'
    publishDir "result/pass/${id}"
    input:
    set id , file (global_paired_alignement) from global_result_paired

    output:
    set id ,file ("${id}pass_global/*_global_paired_alignments.bam") into global_paired_bam

    """
    samtools view -Sb -o ${id}_global_paired_alignments.bam $global_paired_alignement
    samtools flagstat ${id}_global_paired_alignments.bam > stat 

    mkdir ${id}pass_global
    mv ${id}*paired_alignments*.bam ${id}pass_global/.
    mv stat ${id}pass_global/.
    """

}


/*
*Merge the two file bam spliced and global alignement , and sort by read names 
*
*OPTION:
*-n       sort by read names
*/

process mergeBam{
    tag{id}
    errorStrategy 'ignore'
    publishDir "result/pass/${id}"

    input:
    set id , file (spliced_paired_alignement) from spliced_paired_bam
    set id , file (global_paired_alignement) from global_paired_bam

    output:
    set id ,file ("${id}pass_merge/*_all_paired_alignments.bam") into all_paired_bam

    """
    samtools merge -n ${id}_all_paired_alignments.bam $spliced_paired_alignement $global_paired_alignement 
    samtools flagstat ${id}_all_paired_alignments.bam > ${id}stat 

    mkdir ${id}pass_merge
    mv ${id}*paired_alignments*.bam ${id}pass_merge/.
    mv stat ${id}pass_merge/.
    """
}

