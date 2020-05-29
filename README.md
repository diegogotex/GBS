# Protocolo transcriptoma GBS



Autores:<br>
Diego G. Teixeira<br><sup>diego.go.tex@gmail.com</sup>

Viviane B. Nogueira<br><sup>vivianebritonogueiraa@gmail.com</sup>

Roberto Teodoro<br><sup>robteodoro@gmail.com</sup>


## Amostras

Este arquivo é um protocolo para a análise de expressão gênica em pessoas com síndrome de Guillain-Barré. O estudo é uma abordagem transcriptômica de amostras de sangue total.

Ao todo, são 20 amostras: 5 GBS; 9 GBS Recuperado (GBS_rec); 6 Controle (CTL).

| Amostra | Diag. Eletroneuro | Status |coleta | Sexo | Idade |Cod. interno|
| :------: | :-----: |  :-----: | :--------: | :-------: |:-------: |:-------: |
|1  | DESMIELINIZANTE  | GBS | 15/4/2015 |  M | 39  |
|2  | DESMIELINIZANTE  | GBS | 18/5/2015 |  M | 57  |
|13 | DESMIELINIZANTE  | GBS | 23/5/2015 |  H | 46  |
 |25 | INCONCLUSIVO  | GBS | 23/5/2015 |  M | 53  |
|26 | INCONCLUSIVO  | GBS | 12/2/2016 |  H |  81  |
|4  | DESMIELINIZANTE  | GBS | 24/2/2016 |  M | 45  |
|27 | DESMIELINIZANTE  | REC | 7/3/2016 |  M |  31  |
|28 | DESMIELINIZANTE  | REC | 7/3/2016 |  H |  31  |
|7  | DESMIELINIZANTE  | REC | 10/3/2016 |  M | 57  |
|16 | DESMIELINIZANTE  | REC | 11/3/2017 |  H | 47  |
|8  | DESMIELINIZANTE  | REC | 18/3/2016 |  H | 41  |
|30 |INCONCLUSIVO  | REC | 15/3/2016 |  M | 53  |
|9  | DESMIELINIZANTE  | REC | 16/3/2016 |  M | 39  |
|38 | DESMIELINIZANTE  | REC | 19/8/2016 |  H | 15  |
|39 | AXONAL  | REC | 19/8/2016 |  M |  42  |


<center>

| Amostra | Grupo |Amostra | Grupo |
| :------: | :-----: | :--------: | :-------: |
|**1**   |GBS          |**9***     | GBS_rec     |
|**2**   |GBS          |**7***       | GBS_rec     |
|**13**  |GBS          |**16***       | GBS_rec     |
|**25**  |GBS          |**30***      | GBS_rec     |
|**4**   |GBS          |             |             |
|        |             |**27**       |GBS_rec      |
|        |             |**28**       |GBS_rec      |
|        |             |**8**        |GBS_rec      |
|        |             |**38**       |GBS_rec      |
|        |             |**39**       |GBS_rec      |
|**10**  |Controle     |             |             |
|**11**  |Controle     |             |             |
|**19**  |Controle     |             |             |
|**20**  |Controle     |             |             |
|**33**  |Controle     |             |             |
|**34**  |Controle     |             |             |

</center>





<<<<<<< HEAD
# Expressão Diferencial
=======
## Análises
>>>>>>> d07d46a53e2f9a9dc6db1e5ef90b18930370b0b5
### 1 - Peparando as bibiotecas

As bibliotecas estão no formato **fastq** e comprimidas pelo **GZIP**, cada uma delas tem o final **.fastq.gz**. Cada amostra coletada foi sequenciada em 8 lanes distintas, então teremos para cada uma 16 arquivos enumerados de Lane1, Lane2, .. Lane8, com os pares (paired-end) do sequenciamento R1 e R2. Desse modo, teremos 320 arquivos.
Para reduzir a quantidade de arquivos para mapear ao genoma de referência, nós iremos concatenar todas as repetições da mesma amostra em um arquivo só. deixaremos separados somente os pares, ex:

1_lane1_R1.fastq.gz + 1_lane2_R1.fastq.gz + 1_lane3_R1.fastq.gz + 1_lane4_R1.fastq.gz + 1_lane5_R1.fastq.gz + 1_lane6_R1.fastq.gz + 1_lane7_R1.fastq.gz + 1_lane8_R1.fastq.gz = 1.R1.fastq.gz.



```shell
#o loop apresenta outras amostras porque outras amostras foram sequenciadas juto às amostras de GBS
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48;
do
        cat ../Analysis_DualSeq/samples/${K}_lane*R1*.fastq.gz > $K.R1.fastq.gz
        cat ../Analysis_DualSeq/samples/${K}_lane*R2*.fastq.gz > $K.R2.fastq.gz
done

```
### 2 - Mapeamento

Agora que as amostras estão concatenadas, nós iremos mapear as _reads_ ao genoma humano. O genoma humano será obtido no site [ENSEMBL](https://www.ensembl.org/info/data/ftp/index.html), iremos baixar os arquivos DNA(fasta) e GFF3. No ftp do ENSEMBL do genoma, iremos baixar o arquivo **Homo_sapiens.GRCh38.dna.toplevel.fa.gz** o qual contém as sequências dos cromossomos humano (GRCh38.p13).

Em seguida, iremos indexar o genoma Humano, possibilitando utilizá-lo como referência para mapear as _reads_. O mapeapento será realizado com a ferramenta [HISAT2](http://daehwankimlab.github.io/hisat2/download/).

```shell
#descomprimir o genoma
gzip -d Homo_sapiens.GRCh38.dna.toplevel.fa.gz

#formatar a referência para o HISAT2
hisat2-build Homo_sapiens.GRCh38.dna.toplevel.fa Homo_sapiens
````

junto ao genoma de referência, também vamos utilizar o arquivo GFF para criar um pequeno dicionário para traduzir o ENSEMBL gene_ID em SYMBOL. existem infinitas formas de fazer esse processo, mas aqui nos vamos utilizar o próprio GFF.

```shell
 cut -f9 Homo_sapiens.GRCh38.99.GENE.gff3 |\ #selecionando a nona coluna do arquivo GFF
  sed -r 's/;/\t/g' |\ #substituindo o caractere ; por um tab
  cut -f1,2 |\ #selecionando a primeira e segunda colunas com o gene_ID e o respectivo SYMBOL
  sed -r 's/ID=gene://g' |\ #removendo o cabeçalho do identificador
  sed -r 's/Name=//g' \
  >  geneID_symbol_dict.txt
```



Uma vez que o processo de formatação da referência está terminado, iremos mapear as _reads_ ao genoma de referência.

```shell
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48;

do

#mapeando as reads
hisat2 -p 32 --dta-cufflinks -x ref/Homo_sapiens \
           -1 ../amostras_concat/$K.R1.fastq.gz \
           -2 ../amostras_concat/$K.R2.fastq.gz \
           -S sam/$K.sam

#transformando em arquivo bam
samtools sort -@ 32 -o bam/$K.bam sam/$K.sam

done
```


### 3 - Importando a tabela de expressão no R

```R
library(DESeq2)
library(dplyr)
library(plyr)

#reading counts table
countdata <- read.table("~/Dropbox/DualSeq_IMT/GBS_NEW/Counts_GENE.txt",
                        header=TRUE,
                        row.names=1,
                        check.names = FALSE)


countdata <- read.table("/Users/diego/Dropbox/DualSeq_IMT/GBS_NEW/Counts_GENE.txt",
                        header = T,
                        row.names = 1,
                        check.names = F)


#tira os .bam do nome das amostras
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("bam/", "", colnames(countdata))

#filtro pros dados a serem analisados
countdata <- countdata %>% dplyr:: select(grep("01", names(countdata)), #GBS_dis_pair
                                          grep("02", names(countdata)), #GBS_dis_pair
                                          grep("13", names(countdata)), #GBS_dis_pair
                                          grep("25", names(countdata)), #GBS_dis_pair
                                          grep("26", names(countdata)), #GBS_dis_pair
                                          grep("04", names(countdata)), #GBS_dis

                                          grep("09", names(countdata)), #GBS_rec_pair
                                          grep("07", names(countdata)), #GBS_rec_pair
                                          grep("16", names(countdata)), #GBS_rec_pair
                                          grep("30", names(countdata)), #GBS_rec_pair
                                          grep("29", names(countdata)), #GBS_rec_pair
                                          grep("27", names(countdata)), #GBS_rec
                                          grep("28", names(countdata)), #GBS_rec
                                          grep("08", names(countdata)), #GBS_rec
                                          grep("38", names(countdata)), #GBS_rec
                                          grep("39", names(countdata)), #GBS_rec

                                          grep("10", names(countdata)), #control
                                          grep("11", names(countdata)), #control
                                          grep("19", names(countdata)), #control
                                          grep("20", names(countdata)), #control
                                          grep("33", names(countdata)), #control
                                          grep("34", names(countdata)) #control
)


countdata <- as.matrix(countdata)

(condition <- factor(c(rep("GBS",6), #1
                       rep("REC",10), #9
                       rep("CTL",6) #10
)
)
)

(coldata <- data.frame(row.names=colnames(countdata), condition))

#gerar o input pro DESeq
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)


```
<<<<<<< HEAD
### 4 - PCA

Calculando a PCA e pegando as informações para gerar um plit independente do pacote (DESes2):
```R

### PCA

#cauculando a PCA com a função do DESEQ2
rld <- rlogTransformation(dds)
pca <- plotPCA(rld, ntop = nrow(dds),returnData=F)
pcadata <- pca$data

#colocando as cores para cada grupo de amostras de acordo com o PCA$data
col <- c(rep("#C35B52", 6),
         rep("#69B1B7", 10),
         rep("#303030", 6)
         )


condition <- as.data.frame(col, row.names = rownames(pca$data))

condition$cond <- c(rep("GBS", 6), #DIS
                    rep("REC", 10), #REC
                    rep("CTL", 6) #CTL

)
#plotting the PCA with the pre-set colors
plot(pca$data$PC1,
     pca$data$PC2,
     pch=21,
     col="black",
     bg=col,
     cex=2,
     axes=F,
     xlab = ("PC1 (22%)"),
     ylab = ("PC2 (15%)")
)
#box()
abline(h=0,v=0, lty=2)
```
Utilizando as coordenadas da PCA para realizar uma ANOVA, associando o fenótipo das pessoas com a posição na coordenada da PCA:

```R
#Fazendo uma ANOVA nas amostras utilizando a PC1
fitPC1 <- aov(PC1 ~ condition, data=pcadata)
#test tukey for the previous anova
TukeyHSD(fitPC1)
#performing the anova for PC2 coordinates
fitPC2 <- aov(PC2 ~ group, data=pcadata)
TukeyHSD(fitPC2)

#plotting the boxplot for PC1 data
ggplot2::ggplot(pcadata, aes(x=group, y=PC1, fill=group)) +
  scale_fill_manual(values = c("#303030","#C35B52","#69B1B7"))+
  geom_boxplot(fill="white")+
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=1.2)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(face = "bold",size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"))


```

### 5 - Análise de Expressão Diferencial
=======

### 4 - Análise de Expressão Diferencial
>>>>>>> d07d46a53e2f9a9dc6db1e5ef90b18930370b0b5

gerando as tabelas com a análise de expressão diferencial entre os grupos GBSxCTL; GBSxREC; RECxCTL.

```R
#expressao diferencial
dds <- DESeq(dds)

#### DE contraste ####

dict <- read.table("/Users/diego/Dropbox/DualSeq_IMT/GBS_NEW/ref/geneID_symbol_dict.txt",
                   stringsAsFactors = F,
                   header = F)
colnames(dict) <- c("gene_id","SYMBOL")


#extracting the results from Differential expression analysis
contr_GBSctl <- as.data.frame(results(dds, contrast=c('condition','GBS','CTL')))
#creating a new column with the gene SYMBOL name
contr_GBSctl$gene_id <- rownames(contr_GBSctl)
#adicionando os SYMBOL aos respectivos gene_ID
contr_GBSctl$SYMBOL <- subset(dict$SYMBOL, dict$gene_id %in% contr_GBSctl$gene_id)
#removing the NAs from padj column
contr_GBSctl$padj[is.na(contr_GBSctl$padj)] <- 1


contr_GBSrec <- as.data.frame(results(dds, contrast=c('condition','GBS','REC')))
contr_GBSrec$gene_id <- rownames(contr_GBSrec)
contr_GBSrec$SYMBOL <- subset(dict$SYMBOL, dict$gene_id %in% contr_GBSrec$gene_id)
contr_GBSrec$padj[is.na(contr_GBSrec$padj)] <- 1

contr_RECctl <- as.data.frame(results(dds, contrast=c('condition','REC','CTL')))
contr_RECctl$gene_id <- rownames(contr_RECctl)
contr_RECctl$SYMBOL <- subset(dict$SYMBOL, dict$gene_id %in% contr_RECctl$gene_id)
contr_RECctl$padj[is.na(contr_RECctl$padj)] <- 1


#selecting only DEGs from each dataframe
DEGs_GBSctl <- subset(contr_GBSctl, padj <= 0.05 & abs(log2FoldChange) > 1)
DEGs_GBSrec <- subset(contr_GBSrec, padj <= 0.05 & abs(log2FoldChange) > 1)
DEGs_RECctl <- subset(contr_RECctl, padj <= 0.05 & abs(log2FoldChange) > 1)
```

Preparando uma única tabela com os valores de expressão para todas as 3 comparações.

```R
#creating a new column called exp with only 0
contr_GBSctl$exp <- "0"
contr_GBSrec$exp <- "0"
contr_RECctl$exp <- "0"


#setting values of 1 for up regulated genes and 0 fow down into the column exp
for (i in 1:nrow(contr_GBSctl)){
  if ( contr_GBSctl[i,"log2FoldChange"] > 1 &&  contr_GBSctl[i,"padj"] < 0.05){
    contr_GBSctl[i, "exp"] <- "1"
  }
  else if ( contr_GBSctl[i,"log2FoldChange"] < -1 &&  contr_GBSctl[i,"padj"] <= 0.05){
    contr_GBSctl[i, "exp"] <- "-1"
  }
}
#creating a new column with logFC where non DEGs recieve 0 and for DEGs the value remains the same
for (i in 1:nrow(contr_GBSctl)){
  if ( contr_GBSctl[i,"exp"] != 0 ){
    contr_GBSctl[i, "log2FoldChange_mod"] <- contr_GBSctl[i, "log2FoldChange"]
  }
  else {
    contr_GBSctl[i, "log2FoldChange_mod"] <- 0
  }
}


for (i in 1:nrow(contr_GBSrec)){
  if ( contr_GBSrec[i,"log2FoldChange"] > 1 &&  contr_GBSrec[i,"padj"] <= 0.05){
    contr_GBSrec[i, "exp"] <- "1"
  }
  else if ( contr_GBSrec[i,"log2FoldChange"] < -1 &&  contr_GBSrec[i,"padj"] <= 0.05){
    contr_GBSrec[i, "exp"] <- "-1"
  }
}
for (i in 1:nrow(contr_GBSrec)){
  if ( contr_GBSrec[i,"exp"] != 0 ){
    contr_GBSrec[i, "log2FoldChange_mod"] <- contr_GBSrec[i, "log2FoldChange"]
  }
  else {
    contr_GBSrec[i, "log2FoldChange_mod"] <- 0
  }
}


for (i in 1:nrow(contr_RECctl)){
  if ( contr_RECctl[i,"log2FoldChange"] > 1 &&  contr_RECctl[i,"padj"] <= 0.05){
    contr_RECctl[i, "exp"] <- "1"
  }
  else if ( contr_RECctl[i,"log2FoldChange"] < -1 &&  contr_RECctl[i,"padj"] <= 0.05){
    contr_RECctl[i, "exp"] <- "-1"
  }
}
for (i in 1:nrow(contr_RECctl)){
  if ( contr_RECctl[i,"exp"] != 0 ){
    contr_RECctl[i, "log2FoldChange_mod"] <- contr_RECctl[i, "log2FoldChange"]
  }
  else {
    contr_RECctl[i, "log2FoldChange_mod"] <- 0
  }
}

DE_all <- join_all(list(contr_GBSctl,
                        contr_GBSrec,
                        contr_RECctl),
                   by = "gene_id")

#renaming the columns
colnames(DE_all) <- c("baseMean_GBSctl", "log2FoldChange_GBSctl", "lfcSE_GBSctl", "stat_GBSctl", "pvalue_GBSctl", "padj_GBSctl", "gene_id", "SYMBOL", "exp_GBSctl", "log2FoldChange_mod_GBSctl",
                      "baseMean_GBSrec", "log2FoldChange_GBSrec", "lfcSE_GBSrec", "stat_GBSrec", "pvalue_GBSrec", "padj_GBSrec", "SYMBOL1","exp_GBSrec", "log2FoldChange_mod_GBSrec",
                      "baseMean_RECctl", "log2FoldChange_RECctl", "lfcSE_RECctl", "stat_RECctl", "pvalue_RECctl", "padj_RECctl","SYMBOL2", "exp_RECctl", "log2FoldChange_mod_RECctl")
#renaming the rows
rownames(DE_all) <- DE_all$gene_id

#creating a cleaner DF with logFC, padj and exp values
DE_mod_all <- DE_all[,c(9,18,27,6,16,25,10,19,28,7,8)]

```
