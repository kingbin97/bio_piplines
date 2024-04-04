#!/bin/bash
########################################################
usage() {
  echo -e "\e[32m Usage: run_reads_mapping <-r> <-m> <-c> [-t] [-o] [-p]                     \e[0m"
  echo -e "\e[32m        -t       threads           [Default: 60]                                 \e[0m"
  echo -e "\e[32m        -o       work diectory     [Default: $(pwd)]                             \e[0m"
  echo -e "\e[32m        -r       reference index                                                  \e[0m"
  echo -e "\e[32m        -m       mapping soft      [hisat2 | star | bwamem | bwamem2 | bowtie2]   \e[0m"
  echo -e "\e[32m        -c       sample configure  [3 col: Sample fq1 fq2]                        \e[0m"
  echo -e "\e[32m        -p       out prefix        [Default: the same as -m]                      \e[0m"
  echo -e "\e[34m                              _  __  _                   _       _                \e[0m"
  echo -e "\e[34m                             | |/ / (_)                 | |     (_)               \e[0m"
  echo -e "\e[34m                             | ' /   _   _ __     __ _  | |__    _   _ __         \e[0m"
  echo -e "\e[34m                             |  <   | | | '_ \   / _  | | '_ \  | | | '_ \        \e[0m"
  echo -e "\e[34m                             | . \  | | | | | | | (_| | | |_) | | | | | | |       \e[0m"
  echo -e "\e[34m                             |_|\_\ |_| |_| |_|  \__, | |_.__/  |_| |_| |_|       \e[0m"
  echo -e "\e[34m                                                  __/ |                           \e[0m"
  echo -e "\e[34m                                                 |___/                            \e[0m"
}
if [ $# -eq 0 ] || [ "$1" = "-h" ]; then
  usage
  exit 0
fi

############################################################
### Function                                             ###
############################################################

########################################################
# 检查目录是否可写入
########################################################
function check_directory_permissions() {
    local directory="$1"
    if [ ! -d "$directory" ]; then
        echo -e "\e[31m error: 目录 $directory 不存在 \e[0m" >&2
        exit 1
    fi
    if [ ! -r "$directory" ]; then
        echo -e "\e[31m error: 目录 $directory 不可读 \e[0m" >&2
        exit 1
    fi
    if [ ! -w "$directory" ]; then
        echo -e "\e[31m error: 目录 $directory 不可写 \e[0m" >&2
        exit 1
    fi
    echo -e "\e[32m Successful: 目录 $directory 存在且可读写 \e[0m"
}

########################################################
# 检查变量
########################################################
function check_variable() {
  local var_name=$1
  local default_value=$2
  if [ -z "${!var_name}" ]; then
    echo "$var_name 没有值，将被指定为默认值：$default_value"
    eval "$var_name=$default_value"
  else
    echo "$var_name 已经有值：${!var_name}"
  fi
}
########################################################
# 检查软件
########################################################
function check_soft() {
  local soft_name=$1
  echo "检查 $soft_name "
  if ! command -v "$soft_name" &> /dev/null; then
    echo "未找到 $soft_name 命令，请先安装 $soft_name "
    exit 1
  fi
  echo "$soft_name"
}

########################################################
# 检查文件夹
########################################################
function check_and_create_directory() {
  local directory="$1"
  if [ -d "$directory" ]; then
    echo "目录已存在：$directory"
  else
    echo "正在创建目录：$directory"
    if mkdir -p "$directory"; then
      echo "目录创建成功：$directory"
    else
      echo "目录创建失败：$directory，请检查权限或磁盘空间。" >&2
      exit 1
    fi
  fi
}

########################################################
# 检查文件
########################################################
function check_file_readable() {
    local input_reads="$1"
    echo -e "\e[32m  检测文件: ${input_reads}\e[0m"
    if [ -f "$input_reads" ]; then
        if [ -r "$input_reads" ]; then
            echo -e "\e[32m  文件已存在且可读：$input_reads\e[0m"
        else
            echo -e "\e[31m  错误：文件存在但不可读：$input_reads\e[0m" >&2
            return 1  # 或使用 exit 1 来终止脚本
        fi
    else
        echo -e "\e[31m  错误：文件不存在：$input_reads\e[0m" >&2
        return 1  # 或使用 exit 1 来终止脚本
    fi
}
########################################################
# 函数计算指定百分比的内存值
########################################################
function calculate_memory_percentage() {
  local percentage="$1"
  local total_memory_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
  local total_memory_gb=$(echo "scale=2; $total_memory_kb / 1024 / 1024" | bc)
  local percentage_gb=$(echo "scale=2; $total_memory_gb * $percentage / 100" | bc)
  echo "$percentage% 的内存值为: $percentage_gb GB"
}
# calculate_memory_percentage 10
# calculate_memory_percentage 20
########################################################
# run hisat2
########################################################
function run_hisat2() {
  echo -e "\e[32m   运行hisat2: $SAMPLE_RG  \e[0m"
  ###### 1.检查软件
  check_soft hisat2    # 检查 hisat2  软件
  check_soft samtools  # 检查 samtools软件
  ###### 2.输出配置
  echo "INDEX	$REF_INDEX
WORK_DIR	${WORK_DIR}
THREAD	$THREAD
log_DIR	${WORK_DIR}/${OUT_PREFIX}/log
stat_DIR	${WORK_DIR}/${OUT_PREFIX}/stat" >${WORK_DIR}/${OUT_PREFIX}/configure.hisat2.info.txt
  ###### 3. 输出命令
  echo "hisat2 -p $THREAD \
    --dta -x $REF_INDEX \
    --rg-id ${SAMPLE_RG} \
    --rg SM:${SAMPLE_RG} \
    --rg PL:Illumina \
    --summary-file ${WORK_DIR}/${OUT_PREFIX}/stat/${SAMPLE_RG}.hisat2_summary.txt \
    -1 $fq1 -2 $fq2 -S /dev/stdout | \
    samtools sort -m 5G -@10 -o ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.hisat2.sort.bam - " >>${WORK_DIR}/${OUT_PREFIX}/run_hisat2.command.txt
  ###### 4.比对
  hisat2 -p $THREAD \
    --dta -x $REF_INDEX \
    --rg-id ${SAMPLE_RG} \
    --rg SM:${SAMPLE_RG} \
    --rg PL:Illumina \
    --summary-file ${WORK_DIR}/${OUT_PREFIX}/stat/${SAMPLE_RG}.hisat2_summary.txt \
    -1 $fq1 -2 $fq2 -S /dev/stdout | \
    samtools sort -m 5G -@15 -o ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.hisat2.sort.bam -
  ###### 5.建索引
  echo "samtools index"
  samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.hisat2.sort.bam
  samtools flagstat ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.hisat2.sort.bam >${WORK_DIR}/${OUT_PREFIX}/stat/${SAMPLE_RG}.hisat2.sort.bam.flagstat.txt
}

########################################################
# run STAR
########################################################
function run_star() {
  echo -e "\e[32m   运行STAR: $SAMPLE_RG  \e[0m"
  ###### 1.检查软件
  check_soft STAR      # 检查 STAR    软件
  check_soft samtools  # 检查 samtools软件
  check_and_create_directory ${WORK_DIR}/${OUT_PREFIX}/mapping/run/${SAMPLE_RG}
  ###### 2.输出配置
  echo "INDEX	$REF_INDEX
WORK_DIR	${WORK_DIR}
THREAD	$THREAD
log_DIR	${WORK_DIR}/${OUT_PREFIX}/log
stat_DIR	${WORK_DIR}/${OUT_PREFIX}/stat" >${WORK_DIR}/${OUT_PREFIX}/configure.STAR.info.txt
  ###### 3. 输出命令
  echo "STAR \
    --genomeDir $REF_INDEX \
    --runThreadN $THREAD \
    --readFilesIn $fq1 $fq2 \
    --outFileNamePrefix ${WORK_DIR}/${OUT_PREFIX}/mapping/run/${SAMPLE_RG}/${SAMPLE_RG}_STAR. \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --sjdbOverhang 100 \
    --outSAMattrRGline ID:${SAMPLE_RG} PL:Illumina SM:${SAMPLE_RG} 2>> ${WORK_DIR}/${OUT_PREFIX}/log/STAR.${SAMPLE_RG}.log" >>${WORK_DIR}/${OUT_PREFIX}/run_STAR.command.txt
  ###### 4.比对
  STAR \
    --genomeDir $REF_INDEX \
    --runThreadN $THREAD \
    --readFilesIn $fq1 $fq2 \
    --outFileNamePrefix ${WORK_DIR}/${OUT_PREFIX}/mapping/run/${SAMPLE_RG}/${SAMPLE_RG}_STAR. \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --sjdbOverhang 100 \
    --outSAMattrRGline ID:${SAMPLE_RG} PL:Illumina SM:${SAMPLE_RG} 2>> ${WORK_DIR}/${OUT_PREFIX}/log/STAR.${SAMPLE_RG}.log
  mv ${WORK_DIR}/${OUT_PREFIX}/mapping/run/${SAMPLE_RG}/${SAMPLE_RG}_STAR.Aligned.sortedByCoord.out.bam ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}_STAR.Aligned.sortedByCoord.out.bam
  ###### 5.建索引
  echo "samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}_STAR.Aligned.sortedByCoord.out.bam"
  samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}_STAR.Aligned.sortedByCoord.out.bam
  samtools flagstat ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}_STAR.Aligned.sortedByCoord.out.bam >${WORK_DIR}/${OUT_PREFIX}/stat/${SAMPLE_RG}_STAR.Aligned.sortedByCoord.out.bam.flagstat.txt
}

########################################################
# run bowtie2
########################################################
function run_bowtie2() {
  echo -e "\e[32m   运行bowtie2: $SAMPLE_RG  \e[0m"
  ###### 1.检查软件
  check_soft bowtie2   # 检查 bowtie2 软件
  check_soft samtools  # 检查 samtools软件
  ###### 2.输出配置
  echo "INDEX	$REF_INDEX
WORK_DIR	${WORK_DIR}
THREAD	$THREAD
log_DIR	${WORK_DIR}/${OUT_PREFIX}/log
stat_DIR	${WORK_DIR}/${OUT_PREFIX}/stat" >${WORK_DIR}/${OUT_PREFIX}/configure.STAR.info.txt
  ###### 3. 输出命令
  echo "bowtie2 \
    --threads $THREAD \
    -x $REF_INDEX \
    --end-to-end \
    -1 $fq1 -2 $fq2 \
    --rg-id ${SAMPLE_RG} \
    --rg ID:${SAMPLE_RG} \
    --rg SM:${SAMPLE_RG} \
    --rg PL:Illumina 2>> ${WORK_DIR}/${OUT_PREFIX}/log/bowtie2.${SAMPLE_RG}.log | \
    samtools sort -@10 -m 5G -o ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bowtie2.sort.bam - " >>${WORK_DIR}/${OUT_PREFIX}/run_bowtie2.command.txt
  ###### 4.比对
  bowtie2 \
    --threads $THREAD \
    -x $REF_INDEX \
    --end-to-end \
    -1 $fq1 -2 $fq2 \
    --rg-id ${SAMPLE_RG} \
    --rg ID:${SAMPLE_RG} \
    --rg SM:${SAMPLE_RG} \
    --rg PL:Illumina 2>> ${WORK_DIR}/${OUT_PREFIX}/log/bowtie2.${SAMPLE_RG}.log | \
    samtools sort -@15 -m 5G -o ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bowtie2.sort.bam -
  ###### 5.建索引
  echo "samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bowtie2.sort.bam"
  samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bowtie2.sort.bam
  samtools flagstat ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bowtie2.sort.bam >${WORK_DIR}/${OUT_PREFIX}/stat/${SAMPLE_RG}.bowtie2.sort.bam.flagstat.txt

}

########################################################
# run bwamem
########################################################

function run_bwamem() {
  echo -e "\e[32m   运行bwa: $SAMPLE_RG  \e[0m"
  ###### 1.检查软件
  check_soft bwa       # 检查 bwa 软件
  check_soft samtools  # 检查 samtools 软件
  ###### 2.输出配置
  echo "INDEX	$REF_INDEX
WORK_DIR	${WORK_DIR}
THREAD	$THREAD
log_DIR	${WORK_DIR}/${OUT_PREFIX}/log
stat_DIR	${WORK_DIR}/${OUT_PREFIX}/stat" >${WORK_DIR}/${OUT_PREFIX}/configure.STAR.info.txt
  ###### 3. 输出命令
  echo "bwa mem \
    -t $THREAD \
    -R '@RG\tID:${SAMPLE_RG}\tSM:${SAMPLE_RG}\tPL:Illumina' \
    $REF_INDEX \
    $fq1 \
    $fq2 2>>${WORK_DIR}/${OUT_PREFIX}/log/bwamem.${SAMPLE_RG}.log | \
    samtools sort -@10 -m 5G -o ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem.sort.bam - " >>${WORK_DIR}/${OUT_PREFIX}/run_bwamem.command.txt
  ###### 4.比对
  bwa mem \
    -t $THREAD \
    -R "@RG\tID:${SAMPLE_RG}\tSM:${SAMPLE_RG}\tPL:Illumina" \
    $REF_INDEX \
    $fq1 \
    $fq2 2>>${WORK_DIR}/${OUT_PREFIX}/log/bwamem.${SAMPLE_RG}.log | \
    samtools sort -@15 -m 5G -o ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem.sort.bam -
  ###### 5.建索引
  echo "samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem.sort.bam"
  samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem.sort.bam
  samtools flagstat ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bowtie2.sort.bam >${WORK_DIR}/${OUT_PREFIX}/stat/${SAMPLE_RG}.bwamem.sort.bam.flagstat.txt
}

########################################################
# run bwamem2
########################################################
function run_bwamem2() {
echo -e "\e[32m   运行bwa-mem2: $SAMPLE_RG  \e[0m"
  ###### 1.检查软件
  check_soft bwa-mem2    # 检查 hisat2  软件
  check_soft samtools  # 检查 samtools软件
  ###### 2.输出配置
  echo "INDEX	$REF_INDEX
WORK_DIR	${WORK_DIR}
THREAD	$THREAD
log_DIR	${WORK_DIR}/${OUT_PREFIX}/log
stat_DIR	${WORK_DIR}/${OUT_PREFIX}/stat" >${WORK_DIR}/${OUT_PREFIX}/configure.STAR.info.txt
  ###### 3. 输出命令
  echo "bwa-mem2 mem \
    -t $THREAD \
    -R '@RG\tID:${SAMPLE_RG}\tSM:${SAMPLE_RG}\tPL:Illumina' \
    $REF_INDEX $fq1 $fq2 2>>${WORK_DIR}/${OUT_PREFIX}/log/bwamem2.${SAMPLE_RG}.log | \
    samtools sort -@10 -m 5G -o ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem2.sort.bam -" >>${WORK_DIR}/${OUT_PREFIX}/run_bwamem2.command.txt
  ###### 4.比对
  bwa-mem2 mem \
    -t $THREAD \
    -R "@RG\tID:${SAMPLE_RG}\tSM:${SAMPLE_RG}\tPL:Illumina" \
    $REF_INDEX $fq1 $fq2 2>>${WORK_DIR}/${OUT_PREFIX}/log/bwamem2.${SAMPLE_RG}.log | \
    samtools sort -@15 -m 5G -o ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem2.sort.bam -
  ###### 5.建索引
  echo "samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem2.sort.bam"
  samtools index ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem2.sort.bam
  samtools flagstat ${WORK_DIR}/${OUT_PREFIX}/mapping/${SAMPLE_RG}.bwamem2.sort.bam >${WORK_DIR}/${OUT_PREFIX}/stat/${SAMPLE_RG}.bwamem2.sort.bam.flagstat.txt
}


########################################################
# 用于根据 MAP_TYPE 调用对应的函数
########################################################
run_mapping() {
  case "$MAP_TYPE" in
    hisat2)   run_hisat2   ;;
    star)     run_star     ;;
    bowtie2)  run_bowtie2  ;;
    bwamem)   run_bwamem   ;;
    bwamem2)  run_bwamem2  ;;
    *)  echo "Unsupported mapping type: $MAP_TYPE" ;;
  esac
}
###############################################################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################################################
###  Running  
###############################################################################################################################################################################################################################################################################################################
###############################################################################################################################################################################################################################################################################################################

########################################################
### 传参                                             ###
########################################################
while getopts ":o:r:m:c:t:p:" opt; do
  case $opt in
    o) WORK_DIR="$OPTARG";;
    r) REF_INDEX="$OPTARG";;
    m) MAP_TYPE="$OPTARG";;
    c) MAP_CONF="$OPTARG";;
    t) THREAD="$OPTARG";;
    p) OUT_PREFIX="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2;;
  esac
done


########################################################
### 变量                                             ###
########################################################
pwd_path=$(pwd)                        # 当前路径
check_variable WORK_DIR $pwd_path      # 默认工作目录为当前目录
check_and_create_directory  $WORK_DIR  # 检测工作目录
check_directory_permissions $WORK_DIR  # 检测权限
check_variable THREAD 60               # 默认线程60
: ${OUT_PREFIX:=$MAP_TYPE}             # 若OUT_PREFIX为空，则等于MAP_TYPE值

########################################################
### 判断目录，若不存在，则创建                       ###
########################################################
check_and_create_directory ${WORK_DIR}/${OUT_PREFIX}/mapping
check_and_create_directory ${WORK_DIR}/${OUT_PREFIX}/log
check_and_create_directory ${WORK_DIR}/${OUT_PREFIX}/stat

########################################################
### 配置文件                                         ###
########################################################
# 检查配置文件
[ ! -f "$MAP_CONF" ] && echo "Error: 配置文件 $MAP_CONF 不存在" >&2 && exit 1

########################################################
### 从配置文件读取样品信息并运行相应的比对软件       ###
########################################################
while read -r line || [[ -n "$line" ]]; do
  if [[ $line = \#* ]] || [[ -z $line ]]; then
    continue
  fi
  SAMPLE_RG=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $1}')
  fq1=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $2}')
  fq2=$(echo "$line" | awk 'BEGIN{FS="\t";OFS="\t"}{print $3}')
  
  # 检查文件存在性和可读性，如果有问题就跳过此次循环
  if ! check_file_readable "$fq1" || ! check_file_readable "$fq2"; then
    echo "文件问题，跳过样品 $SAMPLE_RG" | tee -a run.log
    continue
  fi

  echo "  样品: $SAMPLE_RG" | tee -a run.log
  echo "  R1: $fq1" | tee -a run.log
  echo "  R2: $fq2" | tee -a run.log
  echo "  工作目录: $WORK_DIR" | tee -a run.log
  echo "  索引文件: $REF_INDEX" | tee -a run.log
  echo "  线程数量: $THREAD" | tee -a run.log
  echo "  样品配置: $MAP_CONF" | tee -a run.log
  echo "  比对软件: $MAP_TYPE" | tee -a run.log
  run_mapping
done < "$MAP_CONF"

