Figures of the Ratio manuscript
================
<gaospecial@gmail.com>
November 11, 2024

# load packages

``` r
library(tidyverse)
library(cowplot)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(vegan)

theme_set(theme_bw())
```

# 读取数据

本研究所使用的数据主要来自于两个实验：首先是两个单独培养体系和三个共培养体系在BIOLOG板中含有的71种碳源中的碳源利用情况，其次是三个共培养体系中所含的两个物种的数量。碳源利用情况以A750的吸光值确定，物种数量则由qPCR实验确定。两部分的数据分别保存在“biolog”和“qPCR”数据文件中。

为了方便后续的分析过程,这两部分的数据都经过了必要的预处理，规范了数据格式。

我们首先载入对应数据，查看一下数据的格式。

``` r
biolog_24h <- read_csv("data/biolog.csv")
qPCR_data <- read_csv(file="data/qPCR.csv")
head(biolog_24h)
#> # A tibble: 6 × 5
#>   ratio0 plate carbon_id  A590  A750
#>   <chr>  <dbl>     <dbl> <dbl> <dbl>
#> 1 equal      1         1 0.191 0.138
#> 2 equal      1         2 0.344 0.181
#> 3 equal      1         3 1.37  0.561
#> 4 equal      1         4 1.25  0.778
#> 5 equal      1         5 0.191 0.137
#> 6 equal      1         6 0.259 0.142
head(qPCR_data)
#> # A tibble: 6 × 6
#>   ratio0 plate carbon_id         EC         PP   ratio1
#>   <chr>  <dbl>     <dbl>      <dbl>      <dbl>    <dbl>
#> 1 less       1        10    177416. 175245243. 0.00101 
#> 2 less       1        11    239521. 148368132. 0.00161 
#> 3 less       1        12  29837437. 142288704. 0.210   
#> 4 less       1        13    645563. 162324668. 0.00398 
#> 5 less       1        14     52481. 142197847. 0.000369
#> 6 less       1        15 164023932. 156667034. 1.05
```

每一行`biolog`数据中包含了初始接种比例（`ratio0`），实验编号（`plate`），碳源的编号，24h时的碳源利用情况（`A750`）等信息。

每一行的qPCR数据中包含了初始接种比例（`ratio0`），实验编号（`plate`），碳源编号（`carbon_id`），*E.
coli* 的定量结果（`EC`），*P. putida*
的定量结果（`PP`），以及最终的物种比例（`ratio1`）。

最后，读取碳源的名称。

``` r
carbon_name <- read_csv("data/carbon.csv")
head(carbon_name)
#> # A tibble: 6 × 2
#>   carbon_id carbon_source
#>       <dbl> <chr>        
#> 1         2 Dextrin      
#> 2         3 D-Maltose    
#> 3         4 D-Trehalose  
#> 4         5 D-Cellobiose 
#> 5         6 Gentiobiosse 
#> 6         7 Sucrose
```

## 数据标准化

``` r
qPCR_data <- qPCR_data %>%
   mutate(ratio0 = factor(ratio0, levels = c("less","equal","more")))

# Normalization
biolog_24h <- biolog_24h %>% 
  mutate(ratio0 = factor(ratio0, levels = c("none","less","equal","more","all"))) %>%
  group_by(plate,ratio0) %>% 
  mutate(A590=A590-A590[carbon_id==1],A750=A750-A750[carbon_id==1]) %>% 
  filter(carbon_id!=1) %>%
  ungroup()
biolog_mono_24h <- biolog_24h %>% 
  filter(ratio0 %in% c("none","all")) %>% 
  mutate(species=factor(ratio0,levels = c("all","none"),labels = c("E. coli","P. putida"))) %>% 
  select(-ratio0)
biolog_coculture_24h <- biolog_24h %>% 
  filter(ratio0 %in% c("less","equal","more")) %>%
  mutate(ratio0 = factor(ratio0, levels = c("less","equal","more")))
```

## 对碳源进行聚类

定义碳源利用的分组：U1，U2，U3。分组采用了层级聚类方法，将碳源分为三类。参见
@ref(table1)。

``` r
M_A750_24h <- biolog_24h %>% mutate(sample=paste(ratio0,plate,sep="-")) %>%
  select(sample,carbon_id,A750) %>%
  spread(key=sample,value=A750) %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var="carbon_id")
k3 <- cutree(hclust(dist(M_A750_24h)),k=3)
carbon_group <-  data.frame(usage=k3) %>%
  rownames_to_column(var="carbon_id") %>%
  mutate(carbon_id=as.numeric(carbon_id)) %>%
  mutate(usage=paste("U",usage,sep=""))
```

## 定义碳源偏好性

根据 *E. coli* 和 *P.putida*
在单独培养时利用碳源能力的差异，将碳源又分为 *E. coli* 偏好性碳源、 *P.
putida* 偏好性碳源及其它碳源。参见 @ref(table1)。

``` r
biolog_mono_A750_24h <- biolog_mono_24h %>% 
  select(plate,carbon_id,species,A750) %>% 
  spread(species,A750) 

PP_prefered <- biolog_mono_A750_24h %>% 
  group_by(carbon_id) %>%  
  summarise(p=t.test(`P. putida`,`E. coli`,alternative = "greater")$p.value) %>% 
  filter(p<0.05)
EC_prefered <- biolog_mono_A750_24h %>% 
  group_by(carbon_id) %>%  
  summarise(p=t.test(`P. putida`,`E. coli`,alternative = "less")$p.value) %>% 
  filter(p<0.05)

carbon_prefer <- data.frame("carbon_id"=carbon_name$carbon_id,
                            "prefer"="none",
                            stringsAsFactors = F)
carbon_prefer[carbon_prefer$carbon_id %in% EC_prefered$carbon_id,"prefer"] <- "EC"
carbon_prefer[carbon_prefer$carbon_id %in% PP_prefered$carbon_id,"prefer"] <- "PP"
```

# 根据碳源利用效率定义相互作用关系

相互作用通常被认为是物种间的固有属性。负的相互作用可能来自于底物竞争、抗生素的合成等，而正的相互作用可能来自于代谢的耦联等。

对于共培养体系而言，其总的碳源利用效率与其中包括的两个物种的碳源利用能力和它们在体系中所占的比例相关。在最简单的情况下，如果总的碳源利用效率高于单独培养时最大的碳源利用效率，那么显然是正的相互作用；如果总的碳源利用效率低于单独培养时最小的碳源利用效率，那么则显然是负的相互作用。但是，实际上共培养体系的碳源利用效率较少出现这两种情况，而多介于最大和最小之间，这时就需要综合考虑两个物种所占有的比例。基于这一思考，我们提出了下列理论模型。

<img src="figures/Figure-definition_of_social_type-1.png" width="70%" />

# Table 1 碳源及其分类

Table 1
包括了实验中的71中碳源的编号、名称，及按照前述分类指标的分类情况。

``` r
library(kableExtra)
left_join(carbon_name,carbon_group) %>% 
  left_join(carbon_prefer) %>% 
  kable()
```

<table>

<thead>

<tr>

<th style="text-align:right;">

carbon_id
</th>

<th style="text-align:left;">

carbon_source
</th>

<th style="text-align:left;">

usage
</th>

<th style="text-align:left;">

prefer
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

2
</td>

<td style="text-align:left;">

Dextrin
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

3
</td>

<td style="text-align:left;">

D-Maltose
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

4
</td>

<td style="text-align:left;">

D-Trehalose
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

5
</td>

<td style="text-align:left;">

D-Cellobiose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

6
</td>

<td style="text-align:left;">

Gentiobiosse
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

7
</td>

<td style="text-align:left;">

Sucrose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

8
</td>

<td style="text-align:left;">

D-Turanose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

9
</td>

<td style="text-align:left;">

Stachyose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

10
</td>

<td style="text-align:left;">

D-Raffinose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

11
</td>

<td style="text-align:left;">

α-D-Lactose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

12
</td>

<td style="text-align:left;">

D-Melibiose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

13
</td>

<td style="text-align:left;">

β-Methyl-D-Glucoside
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

14
</td>

<td style="text-align:left;">

D-Salicin
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

15
</td>

<td style="text-align:left;">

N-Acetyl-D-Glucosamine
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

16
</td>

<td style="text-align:left;">

N-Acetyl-β-Dmannosamine
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

17
</td>

<td style="text-align:left;">

N-Acetyl-D-Galactosamine
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

18
</td>

<td style="text-align:left;">

N-AcetylNeuraminic Acid
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

19
</td>

<td style="text-align:left;">

α-D-Glucose
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

20
</td>

<td style="text-align:left;">

D-Mannose
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

21
</td>

<td style="text-align:left;">

D-Fructose
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

22
</td>

<td style="text-align:left;">

D-Galactose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

23
</td>

<td style="text-align:left;">

3-MethylGlucose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

24
</td>

<td style="text-align:left;">

D-Fucose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

25
</td>

<td style="text-align:left;">

L-Fucose
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

26
</td>

<td style="text-align:left;">

L-Rhamnose
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

27
</td>

<td style="text-align:left;">

Inosine
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

28
</td>

<td style="text-align:left;">

D-Sorbitol
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

29
</td>

<td style="text-align:left;">

D-Mannitol
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

30
</td>

<td style="text-align:left;">

D-Arabitol
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

31
</td>

<td style="text-align:left;">

myo-Inositol
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

32
</td>

<td style="text-align:left;">

Glycerol
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

33
</td>

<td style="text-align:left;">

D-Glucose-6-PO4
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

34
</td>

<td style="text-align:left;">

D-Fructose-6-PO4
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

35
</td>

<td style="text-align:left;">

D-Aspartic Acid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

36
</td>

<td style="text-align:left;">

D-Serine
</td>

<td style="text-align:left;">

U2
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

37
</td>

<td style="text-align:left;">

Gelatin
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

38
</td>

<td style="text-align:left;">

Glycyl-L-Prolin
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

39
</td>

<td style="text-align:left;">

L-Alanine
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

40
</td>

<td style="text-align:left;">

L-Arginine
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

41
</td>

<td style="text-align:left;">

L-Aspartic Acid
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

42
</td>

<td style="text-align:left;">

L-Glutamic
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

43
</td>

<td style="text-align:left;">

L-Histidine
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

44
</td>

<td style="text-align:left;">

L-Pyroglutamic
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

45
</td>

<td style="text-align:left;">

L-Serine
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

46
</td>

<td style="text-align:left;">

Pectin
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

EC
</td>

</tr>

<tr>

<td style="text-align:right;">

47
</td>

<td style="text-align:left;">

D-Galacturonic Acid
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

48
</td>

<td style="text-align:left;">

L-Galactonic Acid Lactone
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

49
</td>

<td style="text-align:left;">

D-Gluconic
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

50
</td>

<td style="text-align:left;">

D-Glucuronic
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

51
</td>

<td style="text-align:left;">

Glucuronamid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

52
</td>

<td style="text-align:left;">

Mucic Acid
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

53
</td>

<td style="text-align:left;">

Quinic Acid
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

54
</td>

<td style="text-align:left;">

D-Saccharic
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

55
</td>

<td style="text-align:left;">

p-Hydroxy-Phenylacetic Acid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

56
</td>

<td style="text-align:left;">

Methyl Pyruvate
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

57
</td>

<td style="text-align:left;">

D-Lactic Acid Methyl Ester
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

58
</td>

<td style="text-align:left;">

L-Lactic Acid
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

59
</td>

<td style="text-align:left;">

Citric Acid
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

60
</td>

<td style="text-align:left;">

α-Keto-Glutaric Acid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

61
</td>

<td style="text-align:left;">

D-Malic Acid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

62
</td>

<td style="text-align:left;">

L-Malic Acid
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

63
</td>

<td style="text-align:left;">

Bromo-Succinic
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

64
</td>

<td style="text-align:left;">

Tween 40
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

65
</td>

<td style="text-align:left;">

y-Amino-ButryricAcid
</td>

<td style="text-align:left;">

U3
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

66
</td>

<td style="text-align:left;">

α-Hydroxy-ButyricAcid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

67
</td>

<td style="text-align:left;">

β-Hydroxy-D,L-ButyricAcid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

68
</td>

<td style="text-align:left;">

α-Keto-ButyricAcid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

69
</td>

<td style="text-align:left;">

Acetoacetic Acid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

70
</td>

<td style="text-align:left;">

Propionic Acid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

PP
</td>

</tr>

<tr>

<td style="text-align:right;">

71
</td>

<td style="text-align:left;">

Acetic Acid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

<tr>

<td style="text-align:right;">

72
</td>

<td style="text-align:left;">

Formic Acid
</td>

<td style="text-align:left;">

U1
</td>

<td style="text-align:left;">

none
</td>

</tr>

</tbody>

</table>

# Figure 2 共培养体系中最终比例的差异

Figure 2 主要展示了共培养体系中最终比例的差异。最终比例以 *E. coli* 与
*P. putida* 定量结果的比值为依据。采用 **ANOVA**
分析三个共培养体系的最终比例是否有差异。P值使用 “BH” 方法做矫正。

``` r
aov_p <- compare_means(ratio1~ratio0,
                       group.by = "carbon_id",
                       data=qPCR_data,
                       method = "anova",
                       p.adjust.method = "BH") %>% 
  arrange(p.adj)
```

为了展示结果，针对于存在显著差异和不存在显著差异的结果，分别选择了5个典型作为示例。

``` r

carbon_name_labeller <- function(x){
  name_of <- carbon_name$carbon_source
  names(name_of) <- carbon_name$carbon_id
  return(as.character(name_of[x]))
}
selected_significant_carbon_id <- c(29,32,36,39,46)
selected_nonsignificant_carbon_id <- c(3,4,7,19,45)
p1 <- ggplot(
  data=filter(qPCR_data,carbon_id %in% selected_significant_carbon_id) %>%
    left_join(aov_p), 
  mapping = aes(ratio0,ratio1)) 
p2 <- ggplot(
  data=filter(qPCR_data,carbon_id %in% selected_nonsignificant_carbon_id) %>%
    left_join(aov_p),
  mapping = aes(ratio0,ratio1)) 

plots <- lapply(list(p1,p2),function(x){
  x + geom_boxplot() + geom_jitter() +
    geom_text(aes(x="equal", y=0.55,label= paste("p.adj = ",p.adj)),check_overlap = T) +
    geom_text(aes(x="less",y=.65,label=carbon_id),color="grey",size=3) +
    facet_wrap(~carbon_id,
               ncol=5,
               labeller = labeller(carbon_id=carbon_name_labeller)) + 
    # stat_compare_means(method="aov") +
    labs(x="",y="final ratio (EC/PP)")
})


plot_grid(plotlist = plots,ncol=1,labels="AUTO")
#> Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : conversion failure on 'α-D-Glucose' in 'mbcsToSbcs': dot
#> substituted for <ce>
#> Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : conversion failure on 'α-D-Glucose' in 'mbcsToSbcs': dot
#> substituted for <b1>
#> Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : conversion failure on 'α-D-Glucose' in 'mbcsToSbcs': dot
#> substituted for <ce>
#> Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : conversion failure on 'α-D-Glucose' in 'mbcsToSbcs': dot
#> substituted for <b1>
#> Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : conversion failure on 'α-D-Glucose' in 'mbcsToSbcs': dot
#> substituted for <ce>
#> Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : conversion failure on 'α-D-Glucose' in 'mbcsToSbcs': dot
#> substituted for <b1>
#> Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : conversion failure on 'α-D-Glucose' in 'mbcsToSbcs': dot
#> substituted for <ce>
#> Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, : conversion failure on 'α-D-Glucose' in 'mbcsToSbcs': dot
#> substituted for <b1>
```

<img src="figures/Figure-unnamed-chunk-5-1.png" width="70%" />

同时，在子图中展示了P值的分布情况。

``` r
p.cutoff <- 0.05
p1 <- ggplot(aov_p,aes(p.adj)) + 
  # geom_histogram(bins=30) + 
  geom_line(stat = "density",lwd=1) +
  geom_density(lwd=0,color=NA,fill="lightblue") +
  geom_vline(xintercept = p.cutoff,lwd=1,lty="dashed",color="firebrick") +
  geom_text(x=0.06,y=0,label=p.cutoff,
            vjust="top",
            hjust="left",
            color="firebrick")

aov_p$p.adj.signif <- cut(aov_p$p.adj,breaks = c(0,0.01,0.05,1),labels = c("**","*","ns"))
freq <- as.data.frame(table(aov_p$p.adj.signif))
p2 <- ggplot(freq,aes(Var1,Freq)) + geom_col() + 
  labs(x="significance of p.adj",y="frequency") + 
  geom_text(aes(label=Freq),vjust=0,nudge_y = 1) +
  ylim(c(0,50))

library(grid)
vp <- viewport(width=0.4,height=0.6,x=0.7,y=0.6)
pushViewport(vp)
```

<img src="figures/Figure-unnamed-chunk-6-1.png" width="70%" />

``` r
p1
print(p2,vp=vp)
```

<img src="figures/Figure-unnamed-chunk-7-1.png" width="70%" />

# Figure S2. Final ratios of all cultures

与@ref(fig2) 对应，Figure S2展示了全部共培养体系的最终比例差异。

``` r
ggplot(qPCR_data,aes(ratio0,ratio1)) + geom_boxplot() +
  geom_text(aes(x="equal",y=.001,label=carbon_id),color="grey",vjust=1,size=3,show.legend = F) +
  # ggpubr::stat_compare_means(method="aov",label="p.signif") +
  facet_wrap(~carbon_id) +
  # geom_jitter() +
  geom_text(aes(x="equal", y=1,label= p.adj),check_overlap = T, data=aov_p) +
  # scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),labels=c("0.001","0.01","0.1","1","10")) +
  xlab("") + ylab("final ratio (EC/PP)") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
      legend.position = "top",
      legend.direction = "horizontal",
      strip.background = element_blank(),  # remove facet label - "strip"
      strip.text = element_blank())
```

<img src="figures/Figure-unnamed-chunk-8-1.png" width="70%" />

# Figure 3. 碳源偏好性对共培养体系最终比例的影响

Figure 3展示了碳源偏好性对共培养体系最终比例的影响。

这个影响可以从两个方面来看，一是对于一个确定初始比例的共培养体系，由于碳源偏好性的差异，可以得到：在*E.
coli*偏好性的碳源中生长有利于提高最终比例；在 *P. putida*
偏好性的碳源中生长则有利于降低最终比例。

另一方面，是在碳源确定的情况下，我们发现，在 *E. coli*
偏好性的碳源中，初始比例差异较大的三个共培养体系的最终比例趋于一致。其总体上差异已经不再显著。而在其它类型碳源中则没有观察到这一现象。

``` r
p <- left_join(biolog_mono_24h,carbon_prefer) %>% 
  filter(prefer!="none") %>%
  ggplot(aes(factor(carbon_id),A750,color=species)) + 
  geom_boxplot() + 
  facet_wrap(~prefer,scales="free_y",ncol=1) +
  labs(y="CUE",x="carbon id")+
  coord_flip() +
  theme(legend.position = c(0.8,0.2),
        legend.text = element_text(face = "italic"))

library(grid)
library(gtable)

# adjust panel height
gp <- ggplotGrob(p)
# gtable_show_layout(gp)
facet.rows <- gp$layout$t[grepl("panel",gp$layout$name)]
y.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                function(l) length(l$range$range))
gp$heights[facet.rows] <- gp$heights[facet.rows] * y.var
```

``` r
p1 <- left_join(qPCR_data,carbon_prefer) %>% ggplot(aes(x=prefer,y=ratio1)) + 
  geom_boxplot() + 
  facet_wrap(~ratio0) + 
  stat_compare_means(method="wilcox.test",comparisons = list(c("EC","none"),c("none","PP"),c("EC","PP"))) +
  labs(x="type of carbon perference", y="final ratio (EC/PP)") +
  ylim(c(NA,1.5))

p2 <- left_join(qPCR_data,carbon_prefer) %>% ggplot(aes(x=ratio0,y=ratio1)) + 
  geom_boxplot() + 
  facet_wrap(~prefer) + 
  stat_compare_means(method="wilcox.test",comparisons = list(c("less","equal"),c("equal","more"),c("less","more"))) +
  labs(x="inital ratio (EC/PP)", y="final ratio (EC/PP)") +
  ylim(c(NA,1.5))
```

``` r
plot_grid(gp,plot_grid(p1,p2,labels = c("B","C"),ncol = 2),labels = c("A",""),ncol = 1,rel_heights = c(1.5,1))
```

<img src="figures/Figure-unnamed-chunk-11-1.png" width="70%" />

# Figure 4. 初始比例对共培养体系代谢能力的调控

在接下来的分析中，我们把三个共培养体系视为一个简单的合成群落，通过分析群落对71种碳源的利用能力比较群落功能的差异。

首先，我们分析3个共培养体系CUE的数值分布是不同的。中位数的比较可以窥其一斑。PCA的结果则明确指出了“equal”和“more”的共培养体系具有与单独培养体系和“less”共培养体系完全不同的CUE谱。这一结果在热图中进一步得到体现。

由此，我们可以发现，U2类的碳源在“equal”和“more”中的利用效率显著提高。而在“less”中没有这种现象。说明初始比例的差异最终能够导致群落功能上走向不同的终点。

``` r
p_cue <- ggplot(biolog_24h,aes(ratio0,A750)) + 
  geom_boxplot() +
  # stat_compare_means(ref.group = ".all.",
  #                    label = "p.format",
  #                    method="wilcox.test") +
  geom_hline(aes(yintercept = median(A750)),lty=2,color="firebrick") +
  geom_text(aes(ratio0,y,label=y),
            inherit.aes = F, 
            # color="firebrick",
            data = biolog_24h %>% group_by(ratio0) %>% 
              summarise(y=median(A750)),
            vjust=-0.3) +
  labs(x="initial ratio (EC/PP)",y="CUE")
```

``` r
# 定义一个在PCA图中绘制置信区间的函数
add_ellipase <- function(p, x="PC1", y="PC2", group="group",
                         ellipase_pro = 0.95,
                         linetype="dashed",
                         colour = "black",
                         lwd = 1,...){
  obs <- p$data[,c(x, y, group)]
  colnames(obs) <- c("x", "y", "group")
  ellipse_pro <- ellipase_pro
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- plyr::ddply(obs, 'group', function(x) {
    if(nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$x, x$y))
    mu <- c(mean(x$x), mean(x$y))
    ed <- sqrt(qchisq(ellipse_pro, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
    })
  names(ell)[2:3] <- c('x', 'y')
  
  ell <- plyr::ddply(ell, plyr::.(group) , function(x) x[chull(x$x, x$y), ])
  p <- p + geom_polygon(data = ell, 
                        aes(x=x,y=y,group = group), 
                        colour = colour,
                        linetype=linetype,
                        lwd =lwd,
                        inherit.aes = F,
                        ...)
  return(p)
}
```

``` r
library(vegan)
pca <- rda(t(M_A750_24h))
percent_var <- pca$CA$eig/pca$tot.chi
df <- scores(pca)$sites  %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="sample") %>%
  separate(sample,c("ratio0","rep"),sep="-",remove = F)
df$ratio0 <- factor(df$ratio0, levels = c("none","less","equal","more","all"))
group <- cutree(hclust(dist(t(M_A750_24h))),k=3)
clustered_group <- as.data.frame(group) %>% tibble::rownames_to_column(var = "sample")
df %<>% left_join(clustered_group) 
p <- ggplot(df, aes(PC1,PC2,label=ratio0,color=ratio0))+
  geom_point(size=3,show.legend = F) +
  scale_color_manual(values = brewer.pal(9,"YlOrRd")[5:9],name="Initial Ratio")+
  xlab(paste0("PC1: ", round(percent_var[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percent_var[2] * 100), "% variance")) +
  directlabels::geom_dl(method="smart.grid",size=3)

p_pca <- add_ellipase(p,alpha=0.1,show.legend = F,lwd=1)
```

``` r
anno_carbon_group <- carbon_group %>% 
  left_join(carbon_prefer) %>%
  column_to_rownames(var="carbon_id")
p_heatmap <- pheatmap(t(M_A750_24h),
         annotation_col = anno_carbon_group[c(2,1)],
         cutree_cols = 3,
         # cutree_rows = 3,
         fontsize_col = 6,
         silent = T)
```

``` r
plot_grid(plot_grid(p_cue,p_pca,ncol = 2,labels = c("A","B")),
          ggplotify::as.ggplot(p_heatmap),ncol = 1,labels = c("","C"),
          rel_heights = c(0.8,1))
```

<img src="figures/Figure-unnamed-chunk-15-1.png" width="70%" />

# Figure 5. Final ratio and CUE in U2 carbon sources

由于U2类碳源表现了独特现象。我们列出U2类碳源中的CUE和最终比例。

``` r
biolog_24h_U2 <- left_join(biolog_24h,carbon_group) %>% filter(usage=="U2") 
hsd_group <- lapply(unique(biolog_24h_U2$carbon_id), function(x){
  m <- aov(A750~ratio0,data=filter(biolog_24h_U2,carbon_id==x))
  library(agricolae)
  g <- HSD.test(m,"ratio0",group=TRUE)$groups
  d <- rownames_to_column(g,var="ratio0")
  d$carbon_id <- x
  return(d[-2])
})
hsd_group <- do.call("rbind",hsd_group)
hsd_group$ratio0 <- factor(hsd_group$ratio0, 
                           levels = c("none","less","equal","more","all"))
# add group on top of boxplot
hsd_group <- biolog_24h_U2 %>% group_by(ratio0,carbon_id) %>% summarize(q3=quantile(A750)[3]) %>% left_join(hsd_group)
```

``` r
u2_p1 <- ggplot(biolog_24h_U2, aes(ratio0,A750)) + 
  geom_boxplot() + 
  geom_text(aes(x="none",y=max(A750)*1.1,label=carbon_id),color="grey",vjust=1,size=3,show.legend = F) +
  geom_text(aes(x=ratio0,y=q3,label=groups),show.legend = F,
            data = hsd_group,inherit.aes = F,
            vjust=0,nudge_y = .2,hjust=0) +
  facet_wrap(~carbon_id,ncol=3,
             labeller = labeller(carbon_id=carbon_name_labeller)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  labs(x="",y="CUE") + 
  # ggpubr::stat_compare_means(method="aov",label="p.format") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = "top",
        legend.direction = "horizontal"
  )
```

``` r
u2_p2 <- left_join(qPCR_data,carbon_group) %>%
  left_join(aov_p) %>%
  filter(usage=="U2") %>%
  ggplot(aes(ratio0,ratio1)) +
  geom_boxplot() +
  facet_wrap(~carbon_id) +
  # geom_jitter() +
  geom_text(aes(x="equal", y=1,label= paste("p.adj = ",p.adj)),check_overlap = T) +
  geom_text(aes(x="less",y=1,label=carbon_id),color="grey",size=3,hjust=1,vjust=1,nudge_x = -0.1,nudge_y = 0.1) +
  facet_wrap(~carbon_id,
             ncol=3,
             labeller = labeller(carbon_id=carbon_name_labeller)) + 
  # stat_compare_means(method="aov") +
  labs(x="",y="final ratio (EC/PP)") +
  # scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),labels=c("0.001","0.01","0.1","1","10")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)  )
```

``` r
plot_grid(u2_p1,u2_p2,labels = "AUTO",ncol=1)
```

<img src="figures/Figure-fig5-1.png" width="70%" />

# Figure S3 所有体系的 CUE

与@ref(fig:fig5) 对应，Figure S3展示了所有的CUE利用情况及差异。

``` r
ggplot(biolog_24h, aes(ratio0,A750)) + 
  geom_boxplot() + 
  geom_text(aes(x="less",y=max(A750)*1.1,label=carbon_id),
            color="grey",
            vjust=1,size=3,show.legend = F) +
  facet_wrap(~carbon_id,ncol=9) +
  scale_y_continuous(breaks = c(0,0.5,1),name = "CUE") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        legend.position = "top",
        legend.direction = "horizontal",
        strip.background = element_blank(),  # remove facet label - "strip"
        strip.text = element_blank())
```

<img src="figures/Figure-unnamed-chunk-19-1.png" width="70%" />

# Figure 6. 初始比例对相互作用的影响

## 计算相互作用类型

定义相互作用的类型。首先，根据 monoculture 的 CUE 和 coculture
中的物种比例，计算理论值。然后，拿理论CUE值与实际测得的值进行比较，从而确定相互作用类型。

``` r
ratio1 <- qPCR_data %>% filter(ratio0 %in% c("less","equal","more")) %>%
  complete(ratio0,carbon_id,plate) %>% 
  group_by(ratio0,carbon_id) %>% 
  select(ratio0,plate,carbon_id,ratio1) %>% 
  mutate(ratio1_mean=mean(ratio1,na.rm = T)) %>% 
  mutate(ratio1=ifelse(is.na(ratio1),ratio1_mean,ratio1)) %>% 
  select(-ratio1_mean)

mono_A750 <- biolog_mono_24h %>% 
  group_by(carbon_id,species) %>% 
  summarise(A750=mean(A750)) %>% 
  spread(key="species",value="A750") 

A750_caculated <- left_join(ratio1,mono_A750) %>% mutate(A750_cac=(`P. putida`+ratio1*`E. coli`)/(1+ratio1))

social <- biolog_coculture_24h %>% select(plate,carbon_id,ratio0,A750) %>%
  left_join(A750_caculated) %>%
  group_by(carbon_id,ratio0) %>% 
  summarise(p_pos=t.test(x=A750,y=A750_cac,alternative = "greater")$p.value,
            p_neg=t.test(x=A750,y=A750_cac,alternative = "less")$p.value)

# social$p_pos.adj <- p.adjust(social$p_pos,method = "BH")
# social$p_neg.adj <- p.adjust(social$p_neg,method = "BH")

# 没有同时显著的情况
# any(social$p_neg.adj < 0.05 & social$p_pos.adj < 0.05)

social <- social %>%
  mutate(social_type=ifelse(
    p_pos<0.05,"+",
    ifelse(p_neg<0.05,"-","unresolved"))
    ) %>% 
  ungroup() %>%
  mutate(ratio0=factor(ratio0,levels = c("less","equal","more")))

biolog_with_social <- biolog_coculture_24h  %>% left_join(social)
```

``` r
carbon_group <- left_join(carbon_group,carbon_prefer)
```

``` r
# social vs ratio0
plots <- lapply(list(c("ratio0","social_type"),c("ratio0","usage","social_type"),c("ratio0","prefer","social_type")), function(x){
  df <- social %>% left_join(carbon_group,by="carbon_id") %>%
    group_by(.dots=x) %>%
    summarise(count=n()) %>%
    mutate(proportion=count/sum(count)) %>%
    mutate(label=paste(round(proportion*100),"%",sep=""))
  ggplot(df,aes_string("ratio0","proportion",fill="social_type")) +
    geom_col() +
    geom_text(aes(label=label),color="white",
              position = position_stack(vjust=0.5),
              size=3) +
    scale_fill_manual(values = c("+"="firebrick","-"="royalblue","unresolved"="grey"),
                     name="effect of mixing") +
   theme(legend.position = "none",
          legend.title = element_text(face="bold")) +
    xlab("") 
})
#> Warning: The `.dots` argument of `group_by()` is deprecated as of dplyr 1.0.0.
#> ℹ The deprecated feature was likely used in the dplyr package.
#>   Please report the issue at <https://github.com/tidyverse/dplyr/issues>.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

legend <- get_legend(plots[[1]] + theme(legend.position = c(0.5,0.7)))
```

``` r
blank <- ggplot()+theme_void()
plot_grid(plot_grid(blank, plots[[1]], legend, rel_widths = c(0.3,0.8,0.7),ncol=3),
          plots[[2]] + facet_wrap(~usage ),
          plots[[3]] + facet_wrap(~prefer),
          ncol=1,
          labels = c("AUTO"))
```

<img src="figures/Figure-social_vs_groups_barplot-1.png" width="70%" />

# Figure S4 所有体系中的相互作用类型

Figure S4 与 Figure 6 对应，展示所有组合中的相互作用类型。

``` r
ggplot(biolog_with_social,aes(ratio0,A750,color=social_type)) + 
    geom_boxplot() + 
    geom_text(aes(x="less",y=max(A750)*1.1,label=carbon_id),vjust=1,color="grey",size=3) +
    facet_wrap(~carbon_id,ncol=9) + 
    scale_color_manual(values=c("+"="firebrick","-"="royalblue","unresolved"="grey"),
                       name="effect of mixing: ") +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    labs(x="",y="CUE") +
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
          legend.position = "top",
          strip.background = element_blank(),  # remove facet label - "strip"
          strip.text = element_blank())
```

<img src="figures/Figure-unnamed-chunk-23-1.png" width="70%" />
