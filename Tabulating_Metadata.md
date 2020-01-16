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
cut -f1,3,18 SRP011546.tsv | perl -pe 's/passage#(\d+)/passage_\1/g; s/ -Cell.*//g; s/[ ]*#\d+//g' | datamash -H --sort --group 1,3 count 2 collapse 2
```
