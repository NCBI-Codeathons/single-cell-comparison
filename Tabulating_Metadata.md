## Projects analyzed
1. ERP003613.tsv
2. SRP011546.tsv
3. SRP050499.tsv
4. SRP057196.tsv
5. SRP061549.tsv
6. SRP066632.tsv

Metadata for each of these projects were downloaded from NCBI using [Entrez Direct](http://bit.ly/entrez-direct) and processed with `xsltproc` using a custom `xsl` file as shown below:
```bash
$ cat studies_list.txt
ERP003613
SRP011546
SRP050499
SRP057196
SRP061549
SRP066632

## fetch metadata in XML format for each project separately
$ for acc in `cat studies_list.txt ` ; do
    esearch -db sra -query ${acc} \
    | efetch -format native -mode xml \
    | xtract -format > ${acc}.xml ;
done

## wrap xml in <root>
$ for f in *.xml; do
    sed -i '/<EXPERIMENT_PACKAGE_SET>/i <root>' $f ;
    echo -e "</root>" >> $f ;
done

## generate metadata table
$ for f in *.xml ; do
    of=$(echo $f | sed 's/xml/tsv/g');
    xsltproc convert_title_desc.as.xsl $f > $of ;
done

## add header to the table
$ for f in *.tsv ; do
    sed -i '1i#SRA_Study\tSRA_Experiment\tSRA_Run\tTotal_Bases\tTotal_Spots\tSRA_Sample\tBiosample_Acc\tdbGaP_Acc\tLibrary_Layout\tPlatform\tLibrary_Source\tSample_Attribs\tTissue_Source\tSource_CellLine\tDevelopmental_Stage\tAnonymized_Name\tCommon_Name\tTitle\tLibrary_Construction_Protocol\tStudy_Title\tStrain\tDevelopmental_Stage\tDisease_Status\tCell_Type\tGEO_Identifer\tSomething_else' $f ;
done
```
Metadata for each project were manually inspected to determine which fields contain the relevant information to produce a list of runs originating from the same type of cells.
### SRP011546
```bash
cut -f1,3,18 SRP011546.tsv \
  | perl -pe 's/passage#(\d+)/passage_\1/g; s/ -Cell.*//g; s/[ ]*#\d+//g' \
  | datamash -H --sort --group 1,3 count 2 collapse 2
  > SRP011546.tbl
```
SRA_Study|Title |SRA_Run_Count|SRA_Run_Accessions
:------------------|:--------------|:-------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SRP011546          |2-cell embryo  |6             |SRR445729, SRR445728, SRR445727, SRR445726, SRR445725, SRR445724
SRP011546          |4-cell embryo  |12            |SRR490972, SRR490971, SRR490970, SRR490969, SRR490968, SRR490967, SRR490966, SRR490965, SRR490964, SRR490963, SRR490962, SRR490961
SRP011546          |8-cell embryo  |20            |SRR490992, SRR490991, SRR490990, SRR490989, SRR490988, SRR490987, SRR490986, SRR490985, SRR490984, SRR490983, SRR490982, SRR490981, SRR490980, SRR490979, SRR490978, SRR490977, SRR490976, SRR490975, SRR490974, SRR490973
SRP011546          |Late blastocyst|30            |SRR491038, SRR491037, SRR491036, SRR491035, SRR491034, SRR491033, SRR491032, SRR491031, SRR491030, SRR491029, SRR491028, SRR491027, SRR491026, SRR491025, SRR491024, SRR491023, SRR491022, SRR491021, SRR491020, SRR491019, SRR491018, SRR491017, SRR491016, SRR491015, SRR491014, SRR491013, SRR491012, SRR491011, SRR491010, SRR491009
SRP011546          |Morulae        |16            |SRR491008, SRR491007, SRR491006, SRR491005, SRR491004, SRR491003, SRR491002, SRR491001, SRR491000, SRR490999, SRR490998, SRR490997, SRR490996, SRR490995, SRR490994, SRR490993
SRP011546          |Oocyte         |3             |SRR445720, SRR445719, SRR445718
SRP011546          |Zygote         |3             |SRR445723, SRR445722, SRR445721
SRP011546          |hESC passage_0 |8             |SRR491045, SRR491043, SRR491042, SRR491041, SRR491040, SRR491039, SRR445731, SRR445730
SRP011546          |hESC passage_10|26            |SRR491090, SRR491089, SRR491088, SRR491087, SRR491086, SRR491085, SRR491084, SRR491083, SRR491082, SRR491081, SRR491080, SRR491079, SRR491078, SRR491077, SRR491076, SRR491075, SRR491074, SRR491073, SRR491072, SRR491071, SRR491070, SRR491069, SRR491068, SRR491067, SRR491066, SRR491065
