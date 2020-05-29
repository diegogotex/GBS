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



# Expressão Diferencial

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
![resultado da PCA](https://github.com/diegogotex/GBS/blob/master/Figs/PCA.png)

<br/>
<br/>
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
![boxplot](https://github.com/diegogotex/GBS/blob/master/Figs/boxplot.png)
<br/>

### 5 - Análise de Expressão Diferencial

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

Agora nós vamos calcular o z-value para os valores de expressão de cada um dos genes para isso, primeiro vamos pegar a normalização da expressão de cada gene com valores em FPKM:

```R
mcols(dds)$basepairs <- subset(countdata_raw$Length, rownames(countdata_raw) %in% rownames(counts(dds)))
fpkm_all  <- fpkm(dds, robust = T)
colnames(fpkm_all) <- c("01.GBS", "02.GBS", "13.GBS", "25.GBS", "26.GBS", "04.GBS",
                        "09.REC", "07.REC", "16.REC", "30.REC", "29.REC", "27.REC", "28.REC", "08.REC", "38.REC", "39.REC",
                        "10.CTL", "11.CTL", "19.CTL", "20.CTL", "33.CTL", "34.CTL")

```

agora vamos calcular o z-value:

```R

log_fpkm <- log2(fpkm_all+0.01)
#calculate the mean row by bow
sum_fpkm <- rowMeans(log_fpkm)
#calculate the standard deviation row by bow
sd_fpkm <- rowSds(log_fpkm)

#creating a dataframe with the same number of rows and columns as log_fpkm df without data
zscore_fpkm <- matrix(data=NA, nrow = nrow(log_fpkm), ncol = ncol(log_fpkm))
colnames(zscore_fpkm) <- colnames(log_fpkm)
rownames(zscore_fpkm) <- rownames(log_fpkm)

#populating the dataframe with zscore values
for (i in 1:nrow(log_fpkm)){
  zscore_fpkm[i,] = (log_fpkm[i,]-sum_fpkm[i])/sd_fpkm[i]
}

```


Plotando o volcano plot para os três contrastes:

```R

par(mfrow=c(1,3))

#GBSvsCTL
with(DE_all, plot(log2FoldChange_mod_GBSctl, -log10(padj_GBSctl), pch=16, axes=T,
                      xlab = "Log2(FC)", ylab = "-Log10(Pvalue-Adjusted)", main = "GBS vs CTL" ,
                      xlim = c(-5,8), ylim = c(0,7), col=NA))

with(subset(DE_all, exp_GBSctl == 0 & abs(padj_GBSctl) != 1), points(log2FoldChange_GBSctl,-log10(padj_GBSctl), pch=16, col=alpha("black", 0.3)))
with(subset(DE_all, exp_GBSctl == 1), points(log2FoldChange_GBSctl,-log10(padj_GBSctl), pch=16, col=alpha("#C35B52", 0.3)))
with(subset(DE_all, exp_GBSctl == -1), points(log2FoldChange_GBSctl,-log10(padj_GBSctl), pch=16, col=alpha("#69B1B7", 0.3)))
abline(h=1.3,col="green", lty = 2, cex= 3)
abline(v=1,col="green", lty = 2, cex= 3)
abline(v=-1,col="green", lty = 2, cex= 3)

#GBSvsREC
with(DE_all, plot(log2FoldChange_mod_GBSrec, -log10(padj_GBSrec), pch=16, axes=T,
                  xlab = "Log2(FC)", ylab = "-Log10(Pvalue-Adjusted)", main = "GBS vs REC" ,
                  xlim = c(-5,8), ylim = c(0,7), col=NA))

with(subset(DE_all, exp_GBSrec == 0 & abs(padj_GBSrec) != 1), points(log2FoldChange_GBSrec,-log10(padj_GBSrec), pch=16, col=alpha("black", 0.3)))
with(subset(DE_all, exp_GBSrec == 1), points(log2FoldChange_GBSrec,-log10(padj_GBSrec), pch=16, col=alpha("#C35B52", 0.3)))
with(subset(DE_all, exp_GBSrec == -1), points(log2FoldChange_GBSrec,-log10(padj_GBSrec), pch=16, col=alpha("#69B1B7", 0.3)))
abline(h=1.3,col="green", lty = 2, cex= 3)
abline(v=1,col="green", lty = 2, cex= 3)
abline(v=-1,col="green", lty = 2, cex= 3)


#GBSvsREC
with(DE_all, plot(log2FoldChange_mod_GBSrec, -log10(padj_RECctl), pch=16, axes=T,
                  xlab = "Log2(FC)", ylab = "-Log10(Pvalue-Adjusted)", main = "REC vs CTL" ,
                  xlim = c(-4,4), ylim = c(0,4), col=NA))

with(subset(DE_all, exp_RECctl == 0 & abs(padj_RECctl) != 1), points(log2FoldChange_RECctl,-log10(padj_RECctl), pch=16, col=alpha("black", 0.3)))
with(subset(DE_all, exp_RECctl == 1), points(log2FoldChange_RECctl,-log10(padj_RECctl), pch=16, col=alpha("#C35B52", 0.3)))
with(subset(DE_all, exp_RECctl == -1), points(log2FoldChange_RECctl,-log10(padj_RECctl), pch=16, col=alpha("#69B1B7", 0.3)))
abline(h=1.3,col="green", lty = 2, cex= 3)
abline(v=1,col="green", lty = 2, cex= 3)
abline(v=-1,col="green", lty = 2, cex= 3)

```
![Volcanoplot](https://github.com/diegogotex/GBS/blob/master/Figs/volcanoplot.png)
<br/>


Plotando o scatterplot comparando **GBS vs REC** com **GBS vs CTL**:

```R
DE_mod_all_sub <- subset(DE_mod_all, DE_mod_all$exp_GBSctl == 1 & DE_mod_all$exp_GBSrec == 1 | DE_mod_all$exp_GBSctl == -1 & DE_mod_all$exp_GBSrec == -1)

with(DE_mod_all, plot(log2FoldChange_mod_GBSctl, log2FoldChange_mod_GBSrec, pch=16, xlab = "GBSctl Genes Log2(FC)", ylab= "GBSrec Genes Log2(FC)", col = "white", xlim = c(-8,8), ylim = c(-8,8), main = NA))
with(subset(DE_mod_all, log2FoldChange_mod_GBSctl > 0 & log2FoldChange_mod_GBSrec > 0), points(log2FoldChange_mod_GBSctl, log2FoldChange_mod_GBSrec, pch=16, col = alpha("#C35B52", 0.3) , cex = 1))
with(subset(DE_mod_all, log2FoldChange_mod_GBSctl > 0 & log2FoldChange_mod_GBSrec < 0), points(log2FoldChange_mod_GBSctl, log2FoldChange_mod_GBSrec, pch=16, col = alpha("gray44", 0.5), cex = 1))
with(subset(DE_mod_all, log2FoldChange_mod_GBSctl < 0 & log2FoldChange_mod_GBSrec < 0), points(log2FoldChange_mod_GBSctl, log2FoldChange_mod_GBSrec, pch=16, col = alpha("#69B1B7", 0.3), cex = 1))
with(subset(DE_mod_all, log2FoldChange_mod_GBSctl < 0 & log2FoldChange_mod_GBSrec > 0), points(log2FoldChange_mod_GBSctl, log2FoldChange_mod_GBSrec, pch=16, col = alpha("gray44", 0.5), cex = 1))
abline(h=0, col="gray44", lty = 3)
abline(v=0, col="gray44", lty = 3)
abline(lm(DE_mod_all_sub$log2FoldChange_mod_GBSrec~DE_mod_all_sub$log2FoldChange_mod_GBSctl, data = DE_mod_all_sub),
       col="green")
text(c(0.5,0.5),c(-1,-2), labels = c(expression('R'^"2"*'=0.75'), expression('p-value=2.2'^"-16")), font = 2, col = alpha("green",0.7), pos = 4)


text(c(2,-2),
     c(7,-7),
     labels = c(nrow(subset(DE_mod_all, log2FoldChange_mod_GBSctl > 0 & log2FoldChange_mod_GBSrec > 0)),
                #nrow(subset(DE_mod_all, log2FoldChange_mod_GBSctl > 0 & log2FoldChange_mod_GBSrec < 0)),
                nrow(subset(DE_mod_all, log2FoldChange_mod_GBSctl < 0 & log2FoldChange_mod_GBSrec < 0))
                #nrow(subset(DE_mod_all, log2FoldChange_mod_GBSctl < 0 & log2FoldChange_mod_GBSrec > 0))
     ),
     #pos = 1,
     cex = 1.5,
     font = 2,
     col = alpha("gray44", 0.7))

```
![Scatterplot](https://github.com/diegogotex/GBS/blob/master/Figs/xplot.png)



### 6 - Análise de enriquecimento

Para fazer a análise de enriquecimento, nós vamos utilizar o pacote **ClusterProfiler** junto ao pacote **org.Hs.eg.db**:

```R
CP_comb_GO <-enrichGO(gene = c(DEGs_GBSctl$SYMBOL, DEGs_GBSrec$SYMBOL),
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      qvalueCutoff = 0.05)

CP_comb_GO_tab <- subset(CP_comb_GO@result, CP_comb_GO@result$qvalue <= 0.05)

emapplot(CP_comb_GO, layout="kk")
```
![Emaplot](https://github.com/diegogotex/GBS/blob/master/Figs/emmaplot.png)
<br/>

Para reduzir a quantidade de Termos enriquecidos, vamos utilizar a ferramenta [**REVIGO**](http://revigo.irb.hr/), com a opção Tiny, para pegarmos somente os termos com menos redundância e mais significativos. Em seguida importamos no R para reduzir na tabela interna:

```R
CP_comb_GO_tab_REVIGO <- read.csv("CP_comb_GO_tab_REVIGO.csv", stringsAsFactors = F)
CP_comb_GO_tab_REVIGO <- CP_comb_GO_tab_REVIGO[CP_comb_GO_tab_REVIGO$eliminated == 0,]
CP_comb_GO_tab_kept <- subset(CP_comb_GO_tab, CP_comb_GO_tab$ID %in% CP_comb_GO_tab_REVIGO$term_ID)

```
<br/>

| ID         | Description                                  | GeneRatio | BgRatio | p.adjust |
|------------|----------------------------------------------|:---------:|:-------:|:--------:|
| GO:0050808 | synapse organization                         |     65    |   394   |  5.3e-13 |
| GO:0042391 | regulation of membrane potential             |     58    |   374   |  2.0e-10 |
| GO:0007416 | synapse assembly                             |     37    |   172   |  2.0e-10 |
| GO:0006836 | neurotransmitter transport                   |     44    |   241   |  3.5e-10 |
| GO:0030198 | extracellular matrix organization            |     51    |   334   |  3.9e-09 |
| GO:0050804 | modulation of chemical synaptic transmission |     56    |   393   |  5.1e-09 |
| GO:0034765 | regulation of ion transmembrane transport    |     61    |   453   |  5.4e-09 |
| GO:0043062 | extracellular structure organization         |     54    |   387   |  2.0e-08 |
| GO:0007215 | glutamate receptor signaling pathway         |     20    |    87   |  5.5e-06 |
| GO:0034776 | response to histamine                        |     7     |    11   |  4.4e-05 |
| GO:0099504 | synaptic vesicle cycle                       |     24    |   162   |  3.7e-04 |
| GO:0050919 | negative chemotaxis                          |     9     |    26   |  3.8e-04 |
| GO:0008037 | cell recognition                             |     20    |   129   |  9.1e-04 |
| GO:0017156 | calcium ion regulated exocytosis             |     19    |   127   |  2.0e-03 |
| GO:0007158 | neuron cell-cell adhesion                    |     6     |    15   |  3.4e-03 |
| GO:0034330 | cell junction organization                   |     30    |   277   |  5.5e-03 |
| GO:0007218 | neuropeptide signaling pathway               |     12    |    66   |  5.6e-03 |
| GO:0034332 | adherens junction organization               |     18    |   134   |  8.1e-03 |
| GO:0048511 | rhythmic process                             |     30    |   287   |  9.0e-03 |
| GO:0070252 | actin-mediated cell contraction              |     15    |   115   |  2.4e-02 |
| GO:0019932 | second-messenger-mediated signaling          |     33    |   360   |  3.0e-02 |
| GO:0071229 | cellular response to acid chemical           |     21    |   199   |  3.7e-02 |
| GO:0043030 | regulation of macrophage activation          |     10    |    66   |  4.1e-02 |
| GO:0089718 | amino acid import across plasma membrane     |     5     |    19   |  4.4e-02 |
| GO:0050678 | regulation of epithelial cell proliferation  |     32    |   369   |  5.8e-02 |
| GO:0043252 | sodium-independent organic anion transport   |     4     |    13   |  5.8e-02 |

<br/>

----
# **Co-expressão**

Para as análises de co-expressão, eu vou carregar novamente todas as bibliotecas e arquivos, já que eu fiz separado das análises de expressão diferencial.

### 1 - Carregando os dados

```R
library(dplyr)

countdata <- read.table("~/Dropbox/DualSeq_IMT/GBS_NEW/Counts_SYMBOL.txt",
                        header = T,
                        row.names = 1,
                        check.names = F)

#tira os .bam do nome das amostras
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("bam/", "", colnames(countdata))

#filtro pros dados a serem analisados
countdata <- countdata %>% dplyr:: select(grep("Length", names(countdata)),
                                          grep("01", names(countdata)), #GBS_dis_pair
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

```

Nessa parte eu irei trabalhar com a expressão normalizada em CPM (Counts Per Million). Utilizarei a função cpm() do pacote [edgeR](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html)

```R
library(edgeR)

#trasnformando a express?o dos genes em CPM
gbs_CPM <- as.data.frame(cpm(countdata[,2:ncol(countdata)], log=T, gene.length=countdata$Length))

```

Para identificar os módulos de genes co-expressos, eu vou utilizar o pacote [CEMiTool](https://www.bioconductor.org/packages/release/bioc/html/CEMiTool.html). Primeiro eu vou fazer a análise para o grupo de amostras coletadas de pessoas com GBS:

### 2 - Identificando módulos de co-expressão

```R
#selecionando as amostras de GBS
gbs_CPM_dis <- gbs_CPM[,c("01","02","13","25","26","04")]

#rodando a análise de identificacao dos m?dulos de co-express?o
cem <- cemitool(expr= gbs_CPM_dis,
                force_beta = T)

#gerando uma tabela mostrando os genes coexpressos e o respectivo modulo
Modules_genes <- module_genes(cem)
```

Vamos calcular o Z-value e plotar o heatmap mostrando o genes nos módulos.

```R
library(pheatmap)
library(genefilter)

#calculando a soma linha por linha
sum_CPM <- rowMeans(gbs_CPM)
#calculate the standard deviation row by bow
sd_CPM <- rowSds(gbs_CPM)

#creating a dataframe with the same number of rows and columns as log_fpkm df without data
zscore_CPM <- matrix(data=NA, nrow = nrow(gbs_CPM), ncol = ncol(gbs_CPM))
colnames(zscore_CPM) <- colnames(gbs_CPM)
rownames(zscore_CPM) <- rownames(gbs_CPM)

zscore_CPM <- as.data.frame(zscore_CPM)

#populating the dataframe with zscore values
for (i in 1:nrow(gbs_CPM)){
  zscore_CPM[i,] = as.vector((gbs_CPM[i,]-sum_CPM[i])/sd_CPM[i])
  #print((gbs_CPM[i,]-sum_CPM[i])/sd_CPM[i])
}

##########################
### plotando o heatmap ###
##########################

mat_col = data.frame(Condition=factor(c(rep("GBS",6), #01,02,13,25,26,04
                                        rep("REC",10), #09.07,16,30,29,27,28,08,38,39
                                        rep("CTL",6) #10,11,19,20,33,34
                                        )
                                      )
                     )

row.names(mat_col) = colnames(zscore_CPM)

mat_row <- Modules_genes
row.names(mat_row) <- mat_row$genes
mat_row <- as.matrix(mat_row)
mat_row <- as.data.frame(mat_row[,2])
colnames(mat_row) <- "Modules"
mat_row$Modules <- as.factor(mat_row$Modules)


mat_colors <- list(Condition = c('GBS' = "#C35B52",
                                 'REC' = "#69B1B7",
                                 'CTL' = "#303030"),
                   Modules = c('M1'= "#b40000ff",'M2'= "#b34200ff",'M3'= "#b28400ff",
                               'M4'= "#9bb200ff",'M5'= "#57b300ff",'M6'= "#14b400ff",
                               'M7'= "#00b42bff",'M8'= "#00b46fff",'M9'= "#00b4b4ff",
                               'M10'= "#006fb4ff",'M11'= "#002bb3ff",'M12'= "#1400b4ff")
)


pheatmap(t(subset(zscore_CPM, rownames(zscore_CPM) %in% Modules_genes$genes)),
         color = colorRampPalette(c("blue","gray20","red"))(n = 500),   #colors from RColorBrewer,
         breaks =  c(seq(-3.8,-0.97,length=166),                            # for blue
                     seq(-0.971,0.99,length=166),                         # for black
                     seq(0.991,4.4,length=166)
         ),
         show_rownames = T,
         show_colnames = F,
         treeheight_col = 0,
         annotation_col = mat_row,
         annotation_row = mat_col,
         annotation_colors = mat_colors,
         border_color = NA,
         #clustering_method = "ward.D2",
         legend = T,
         fontsize_col = 10)


```
![](https://github.com/diegogotex/GBS/blob/master/Figs/heatmap_ModGBS.png)
<br/>

### 3 - Enriquecimento com termos de GO

Agora que eu tenho os módulos de genes co-expressos, eu vou realizar uma análise de enriquecimento com termos de GO utilizando o pacote [enrichR](https://cran.r-project.org/package=enrichR)

```R
library(enrichR)

#enriquecimento com termos de GO para processos biologicos
enrichr_M1 <- enrichr(Modules_genes[Modules_genes$modules %in% "M1",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
#selecionando termos com pvalor ajustado unferior a 0.05
enrichr_M1 <- subset(enrichr_M1, enrichr_M1$Adjusted.P.value < 0.05)
#adicionando uma coluna identificando o m?dulo
enrichr_M1$Module <- rep("M1", nrow(enrichr_M1))

enrichr_M2 <- enrichr(Modules_genes[Modules_genes$modules %in% "M2",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M2$Module <- rep("M2", nrow(enrichr_M2))

enrichr_M3 <- enrichr(Modules_genes[Modules_genes$modules %in% "M3",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M3$Module <- rep("M3", nrow(enrichr_M3))

enrichr_M4 <- enrichr(Modules_genes[Modules_genes$modules %in% "M4",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M4$Module <- rep("M4", nrow(enrichr_M4))

enrichr_M5 <- enrichr(Modules_genes[Modules_genes$modules %in% "M5",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M5$Module <- rep("M5", nrow(enrichr_M5))

enrichr_M6 <- enrichr(Modules_genes[Modules_genes$modules %in% "M6",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M6$Module <- rep("M6", nrow(enrichr_M6))

enrichr_M7 <- enrichr(Modules_genes[Modules_genes$modules %in% "M7",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M7 <- subset(enrichr_M7, enrichr_M7$Adjusted.P.value < 0.05)
enrichr_M7$Module <- rep("M7", nrow(enrichr_M7))

enrichr_M8 <- enrichr(Modules_genes[Modules_genes$modules %in% "M8",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M8$Module <- rep("M8", nrow(enrichr_M8))

enrichr_M9 <- enrichr(Modules_genes[Modules_genes$modules %in% "M9",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M9 <- subset(enrichr_M9, enrichr_M9$Adjusted.P.value < 0.05)
enrichr_M9$Module <- rep("M9", nrow(enrichr_M9))

enrichr_M10 <- enrichr(Modules_genes[Modules_genes$modules %in% "M10",]$genes,
                      "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M10$Module <- rep("M10", nrow(enrichr_M10))

enrichr_M11 <- enrichr(Modules_genes[Modules_genes$modules %in% "M11",]$genes,
                       "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M11$Module <- rep("M11", nrow(enrichr_M11))

enrichr_M12 <- enrichr(Modules_genes[Modules_genes$modules %in% "M12",]$genes,
                       "GO_Biological_Process_2018")$GO_Biological_Process_2018
enrichr_M12$Module <- rep("M12", nrow(enrichr_M12))

#juntando todas as tabelas com os termos de GO dos m?dulos
GO_enrichment <- rbind(enrichr_M1, enrichr_M2, enrichr_M3, enrichr_M4, enrichr_M5, enrichr_M6,
                       enrichr_M7, enrichr_M8, enrichr_M9, enrichr_M10, enrichr_M11, enrichr_M12)

#removendo colunas in?teis
GO_enrichment <- GO_enrichment[,-c(5,6)]

#selecionando termos com pvalor ajustado inferior a 0.05
GO_enrichment_Pv <- subset(GO_enrichment, GO_enrichment$Adjusted.P.value < 0.05)
```

Dentro da tabela, só foram identificados termos significantemente enriquecidos para os módulos M1, M7 e M9. Os genes do módulo M9 relacionam-se com os termos _peptidyl-tyrosine modification (GO:0018212)_ e _embryonic digestive tract development (GO:0048566)_ já os módulos M1 e M7 apresentam termos que parecem relevantes para as análises de GBS. No entando, o M1 tem 91 termos e para reduzir eu vou utilizar novamente o REVIGO e uma função que eu mesmo criei, chamada de [**cool_table**](https://github.com/diegogotex/Rcodes). Essa função serve para tornar a tabela do enrichR com uma visualização melhor.

Não achei necessário reduzir o número de termos para o módulo M7 já que ele tem somente 10 termos.

```R
source("https://github.com/diegogotex/Rcodes/blob/master/cooltable.R")
library(ggplot2)
library(tidyr)


#preparando uma tabela com termos de GO para o M1
#a funcaoo cool_table tem que ser carregada
enrichr_M1 <- cool_table(enrichment.obj = enrichr_M1, enrich.type = 2, up = Modules_genes$genes, down = NA)
#pegando os IDs dos termos de GO
enrichr_M1$GO_id <- separate(data = enrichr_M1, col = Term, into = c("left", "right"), sep = "\\(")$right
enrichr_M1$GO_id <- gsub(")", "", enrichr_M1$GO_id)

#selecionando os termos de GO com menos redundancia
enrichr_M1_REVIGO <- read.csv("enrich_M1_REVIGO.csv", stringsAsFactors = F)

enrich_M1_kept <- subset(enrichr_M1, enrichr_M1$GO_id %in% enrichr_M1_REVIGO[enrichr_M1_REVIGO$eliminated == 0,]$term_ID)

ggplot(data = enrich_M1_kept, aes(x=reorder(Term,-Adjusted.P.value,), y = -log10(Adjusted.P.value)))+
  geom_col(aes(fill = -log10(Adjusted.P.value))) +
  scale_fill_gradient2(low = "white",
                       high = "#b40000ff")+
  coord_flip() +
  labs(x = "",
       y = "",
       title = "GO (BP) Terms for M1")+
  theme_bw()

ggplot(data = enrichr_M7, aes(x=reorder(Term,-Adjusted.P.value,), y = -log10(Adjusted.P.value)))+
  geom_col(aes(fill = -log10(Adjusted.P.value))) +
  scale_fill_gradient2(low = "white",
                       high = "#00b42bff")+
  coord_flip() +
  labs(x = "",
       y = "",
       title = "GO (BP) Terms for M7")+
  theme_bw()

```

### Rede PPI com dados de Co-expressão

Essa etapa é a mais demorada pois depende do download de alguns dados do ENSEMBL e do STRING.

**1.** Na primeira parte eu vou importar um dataframe que tem a associação entre o gene_id e o SYMBOL, eu gerei isso a partir do arquivo .GFF.
```R
ref <- read.table("/Users/diego/Dropbox/DualSeq_IMT/GBS_NEW/ref/geneID_symbol_dict.txt",
                  stringsAsFactors = F,
                  header = F)
colnames(ref) <- c("gene_id","SYMBOL")

```

**2.** Em seguida eu vou puxar as informações do ensembl

```R
library(biomaRt)
# Conectar a um database
ensembl <- useMart("ensembl")

# selecionando o dataset
ensembl <- useDataset("hsapiens_gene_ensembl", ensembl)


# retirar outras informacoes sobre nossa lista de proteinas. Rode o comando abaixo:
query <- getBM(attributes=c("ensembl_gene_id","ensembl_peptide_id","entrezgene_id","hgnc_symbol"),
               mart=ensembl, filters="ensembl_gene_id", values=ref$gene_id)

```

**3.** Pegando as informações de interação PPI no STRING

```R
library(STRINGdb)

#gerando um objeto do string associado a versao e ao organismo
string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory="" )

#mapeando os IDs no string
onto_query_mapped <- string_db$map(query, colnames(query)[2], removeUnmappedRows = TRUE )

# Recupere as interacoes do string online
string_query_dat <- string_db$get_interactions(onto_query_mapped$STRING_id)

```

**4.** Vou remover as interações com os scores baixos (inferiores a 0.9), com base nas informações de **coexpression**, **experiments** e **database** do STRING 10.

O script para fazer essa filtragem foi escrito por Iara, chama-se [source_string.R](https://github.com/diegogotex/Rcodes/blob/master/source_stringdb.R).

```R

# Carregue a funcao 'combinescores' do arquivo 'source_stringdb.R', filtre as interacoes segundo as evidencias que voce julgue relevantes
# e combine os scores associados.
source("source_stringdb.R")
res_query <- combinescores(string_query_dat, evidences = c("coexpression","experiments","database"), confLevel=0.9)
# Remova as  informacoes sem utilidade e remova os id/links nao anotados.
res_query$from <- substr(res_query$from,6,1000)
res_query$to <- substr(res_query$to,6,1000)
query_inter <- res_query

#gerando uma tabela de interacao por SYMBOL
query_inter <- merge(query_inter, query, by.x = "from", by.y = "ensembl_peptide_id", all.x = T)
query_inter <- merge(query_inter, query, by.x = "to", by.y = "ensembl_peptide_id", all.x = T)
query_inter <- query_inter[,c(6,9,3)]
colnames(query_inter) <- c("from", "to","combined_score")



#selecionando somente as interacoes no M1
idx1 <- query_inter$from %in% Modules_genes$genes
idx2 <- query_inter$to %in% Modules_genes$genes
inter <- query_inter[idx1 & idx2,]
inter <- inter[,c(1,2)]

#nomeando as colunas de moco que o CEMiTools reconheca
colnames(inter) <- c("gene1symbol", "gene2symbol")

```

**5.** Associando as informações de PPI ao dado de co-expressão do CEMiTools
```R
interactions_data(cem) <- inter # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module

#plotando o módulo 1
plots[1]
```
![](https://github.com/diegogotex/GBS/blob/master/Figs/M_ALL.png)

### 4 - Identificando os tipos celulares

Para essa etapa das análises, irei utilizar o (CTen)[http://www.influenza-x.org/~jshoemaker/cten/upload.php]

O dado de entrada é a tabela de Gene e o respectivo módulo.

```R
CTen.GBS <- read.csv("~/Dropbox/DualSeq_IMT/GBS_NEW/Coexp/CellType/CTen/Enrichment_GBS.csv",
                     stringsAsFactors = F,
                     row.names = 1)

CTen.GBS <- CTen.GBS[,c("M1","M2","M3","M4","M5","M6",
                        "M7","M8","M9","M10","M11","M12")]

#colocando os valores abaixo de 1 como 0
CTen.GBS[CTen.GBS < 1] <- 0

pheatmap(#mat = CTen.GBS,
  mat = CTen.GBS[!rowSums(CTen.GBS) <  1.5,],
  cluster_rows = T,
  cluster_cols = F,
  fontsize_row = 8,
  border_color = "white")


CTen.REC <- read.csv("~/Dropbox/DualSeq_IMT/GBS_NEW/Coexp/CellType/CTen/Enrichment_REC.csv",
                     stringsAsFactors = F,
                     row.names = 1)

CTen.REC <- CTen.REC[,c("M1","M2","M3","M4","M5","M6","M7","M8"
                        ,"M9","M10","M11","M12","M13","M14","M15","M16","Not.Correlated")]

CTen.REC <- CTen.REC[,-17]

CTen.REC[CTen.REC < 1] <- 0

pheatmap(mat = CTen.REC[!rowSums(CTen.REC) <  1.5,],
         cluster_rows = T,
         cluster_cols = F,
         fontsize_row = 8,
         border_color = "white")

```
