#!/bin/bash

# to do all except core/raxml:
# ./dev/eastconst.sh test/src/* src/* src/core/* src/core/pll/* src/io/* src/net/* src/pipeline/* src/sample/* src/seq/* src/tree/* src/util/*

declare -a arr=(
	"auto" "double" "int" "unsigned int" "size_t" "bool" "uint64_t"
	"std::string" "std::exception" "std::basic_string< char >"
	"std::vector< size_t >" "std::vector< T >" "std::vector< std::string >" "std::vector< int >" "std::vector< double >"
	"Timer<>" "T" # timer has to come before T
	"PQuery< Slim_Placement >" "PQuery" # specific before general
	"Sample< T >" "Sample" # specific before general
	"Options" "Work" "Placement" "Slim_Placement" "Sequence"
	"uchar" "mask_type" "MPI_Comm" "seqid_type" "value_type"
	"stack_type" "hook_type" "Function" "stage_f" "token_status"
	"genesis::sequence::Sequence" "Matrix< T >"
	"MSA_Info" # broken but fixable
)
# these don't work:
#

for i in "${arr[@]}"
do
	rpl -x .cpp -x .hpp "const ${i}" "${i} const" $@
done

# only those that match with &
for i in "MSA"
do
	rpl -x .cpp -x .hpp "const ${i}&" "${i} const&" $@
done

# fix incorrects
rpl -x .cpp -x .hpp --quiet "MSA_Info const::mask_type" "MSA_Info::mask_type const" $@