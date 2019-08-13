
{if (ANN!="") {

	{n=split(ANN,AR,",")}
	{nb_ANN=0}
	{sep_ANN=""}
	{RA=""}
	{for (K in AR) {
		{nb_ANN++}
		{A=AR[nb_ANN]}
		{if (A!="") {
			{R="~"}
			{LEN=length(A)+2}
			{if (match($8,"[;\t]"A"=[^;\t]*")){R=substr($8,RSTART+LEN,RLENGTH-LEN)}}
			{RA=RA sep_ANN R}
			{sep_ANN="\t"}
		}}
	}}
	{RA=RA "\t"}

}}

{if (FIELDS!="") {

	{if (FIELDS=="INFO") {

		FIELDS_OUTPUT=$8

	} else {

		{o=split(FIELDS,FIELDS_ARRAY,",")}
		{nb_FIELDS=0}
		{sep_FIELDS=""}
		{FIELDS_OUTPUT=""}
		{for (L in FIELDS_ARRAY) {
			{nb_FIELDS++}
			{B=FIELDS_ARRAY[nb_FIELDS]}
			{if (B!="") {
				{LEN=length(B)+2}
				{F=""}
				{F_TRUE=0}
				{if (match(";"$8";","[;\t]"B"=[^;\t]*")){F=substr(";"$8";",RSTART+LEN,RLENGTH-LEN)}}
				{if (match(";"$8";","[;\t]"B"[;\t]")){F="TRUE"; F_TRUE=1}}
				{if (FORMAT=="TSV") {
					{FIELDS_OUTPUT=FIELDS_OUTPUT sep_FIELDS F}
					{sep_FIELDS="\t"}
				} else {
					{if (F!="") {
						{if (F_TRUE) {
							{FIELDS_OUTPUT=FIELDS_OUTPUT sep_FIELDS B}
							{sep_FIELDS=";"}
						} else {
							{FIELDS_OUTPUT=FIELDS_OUTPUT sep_FIELDS B "=" F}
							{sep_FIELDS=";"}
						}}
					}}

				}}
			}}
		}}
	}}

	{FIELDS_OUTPUT=FIELDS_OUTPUT "\t"}

}}

{l=split($0,LINE_ARRAY,"\t")}

{S=""}
{sep_SAMPLES=""}
{ for (i=9; i<=l; i++) {
	{S=S sep_SAMPLES $i}
	{sep_SAMPLES="\t"}
}}

{print RA $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t" FIELDS_OUTPUT S}
