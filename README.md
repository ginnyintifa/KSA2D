KSA2D

## 1 Introduction 

This package aims to help users identify differentially expressed proteins and phosphorylation sites with significantly increased/reduced levels. 
The main algorithm for 1-dimensional discovery is based on
>"Efron, Bradley, et al. "Empirical Bayes analysis of a microarray experiment." Journal of the American statistical association 96.456 (2001): 1151-1160.
Efron, B., Tibshirani, R., Storey, J. "

We extended the algorithm to 2-dimensional setting to discover joint variation of kinase and substrates.  (cite our paper)

This package can achieve the following functions. 
## 2 Installation and preparation 
KSA2D can be downloaded and installed in R. Installation of GPD requires devtools as a prerequisite:

```{r}
install.packages("devtools")
```
Next, install KSA2D by:

```{r}
library("devtools")
devtools::install_github("ginnyintifa/KSA2D")
library(KSA2D)
```

Suppose we wish to identify differentially expressed proteins and phosphorylation sites with different abundance levels between 2 time points(for example 5 min vs 0 min, or 30 min vs 0 min) upon drug treatment.
Input files for KSA2D include 

* protData.tsv: Data frame of proteome data, each row is a protein, columns are samples' proteome data at two time points. 
* psiteData.tsv: Data frame of phosphoproteome data, each row is a protein, columns are samples' proteome data at two time points. 
* Files for annotating proteins and mapping kinase-substrate relationships can be downloaded from our [WEBSITE](https://drive.google.com/drive/folders/1k3YjPWVav5z9zbuuiP14kiiB5Fu856DX).

Please notice that there are specific format requirements for both protData and psiteData. 

protData.tsv:

```
accession	proteinName	geneName	p_0min_1	p_5min_1	p_30min_1	p_0min_6	p_5min_6	p_30min_6	p_0min_10	p_5min_10	p_30min_10
A0AVF1	Intraflagellar transport protein 56	TTC26	46017.26824	51956.02091	58420.38931	56927.73264	51625.471	49119.23734	45616.15195	52712.8303	63683.23818
A0AVT1	Ubiquitin-like modifier-activating enzyme 6	UBA6	1504222.055	1526087.312	1583272.622	1786412.413	1521421.6	1660424.492	1862740.935	1669284.71	1830118.522
A0FGR8	Extended synaptotagmin-2	ESYT2	420018.4176	420620.0213	452432.1388	493396.577	444219.1772	445948.4243	412030.377	410524.8839	449422.7346
A0JNW5	UHRF1-binding protein 1-like	UHRF1BP1L	76838.34679	92835.55988	98490.51785	91713.63145	89243.61806	81828.32986	80589.69146	82753.17311	80340.36924
A0MZ66	Shootin-1	SHTN1	222247.0948	240172.1188	231529.3365	227574.1463	222079.1595	189088.9893	199738.3231	208230.3438	196291.5963
A0PJW6	Transmembrane protein 223	TMEM223	27148.58515	38551.20212	22571.20742	26517.2227	45615.93667	23239.16409	35198.45633	40668.53716	24611.93232
A0PJZ3	Glucoside xylosyltransferase 2	GXYLT2	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0PK00	Transmembrane protein 120B	TMEM120B	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1A4S6	Rho GTPase-activating protein 10	ARHGAP10	171025.0571	175361.4531	181214.1891	182041.221	179542.2988	146204.2002	156305.7084	182256.501	183672.5068
A1L0T0	Acetolactate synthase-like protein	ILVBL	1032644.238	1087723.099	819733.3501	846006.855	903808.5655	772441.0415	836886.6	846006.855	915899.2739
A1L157	Tetraspanin-11	TSPAN11	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1L188	NADH dehydrogenase [ubiquinone] 1 alpha subcomplex assembly factor 8	NDUFAF8	132008.7509	141237.367	144729.6935	184930.4746	151341.6131	142467.9572	115708.6773	119140.2275	111112.1263
A1L190	Synaptonemal complex central element protein 3	SYCE3	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1L390	Pleckstrin homology domain-containing family G member 3	PLEKHG3	NA	NA	NA	NA	NA	NA	NA	NA	NA```
```

psiteData.tsv

```
accession	res	phosphoPos	resPos	p_0min_1	p_5min_1	p_30min_1	p_0min_6	p_5min_6	p_30min_6	p_0min_10	p_5min_10	p_30min_10
A0AUZ9	S	714	S714	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0AVT1	S	743	S743	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0FGR8	S	676	S676	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0FGR8	S	755	S755	94373.12342	97026.30243	110925.6985	90732.56506	80781.29388	105189.7109	81291.14389	75084.29861	97007.66739
A0FGR8	S	758	S758	69627.13535	72271.93972	83042.7429	73565.48702	63822.18018	75650.96783	64856.90649	61105.22045	76891.04292
A0FGR8	S	761	S761	54198.4238	57437.10626	64978.40902	54452.3492	51086.13359	57838.47671	52319.67621	50159.36434	63971.92493
A0JNW5	S	414	S414	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0JNW5	S	891	S891	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0JNW5	S	935	S935	65809.40238	62366.92346	79307.50705	50697.43653	57549.96835	55859.11016	65698.05178	60447.64663	66261.31247
A0JNW5	S	953	S953	143342.418	145449.9812	171269.5045	141138.3868	165244.9491	127107.0852	138662.2747	123805.6024	140430.0909
A0JNW5	S	989	S989	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0JNW5	S	1021	S1021	16307.52362	20325.59688	19720.6684	12824.48374	21837.91807	8219.436248	17666.29844	22732.61225	13761.68469
A0MZ66	S	506	S506	30018.71103	20774.59054	26010.29005	23346.59842	20598.16144	25423.04104	24359.40277	23711.62359	37542.77554
A0MZ66	S	534	S534	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1A4S6	S	376	S376	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1A4S6	S	591	S591	80777.86895	95065.19436	94813.62826	85048.89813	78682.87052	73094.56	75847.76494	83476.94901	81184.30477
A1A4S6	S	600	S600	77502.37689	81049.07685	91474.60108	74274.83926	62461.21705	58011.50925	72228.55304	69388.1741	75413.04207
A1L020	S	338	S338	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1L390	S	76	S76	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1L390	S	448	S448	NA	NA	NA	NA	NA	NA	NA	NA	NA
```

Please format your datasets following the above samples. 

## 3 Functions 

### 3.1 Data cleaning and preparation 
```protData_prep``` prepares the protData by adding a fudge factor to the entire data and converting it to log2 scale. 

```{r}

protData_p = protData_prep(protData_filename = "path/to/protData.tsv",
fudge_factor = 0.01)
```

Return value of this function is a data frame ```protData_p```

```
prot_names	p_0min_1	p_5min_1	p_30min_1	p_0min_6	p_5min_6	p_30min_6	p_0min_10	p_5min_10	p_30min_10
A0AVF1_TTC26	15.81274713	15.95441177	16.09428084	16.06316741	15.94688189	15.8884745	15.80265806	15.9715052	16.19891882
A0AVT1_UBA6	20.53161326	20.55227581	20.60495682	20.77792458	20.5478915	20.67311404	20.83790691	20.6807391	20.8125753
A0FGR8_ESYT2	18.71920122	18.72121097	18.82368319	18.94574684	18.79791681	18.80338032	18.69224715	18.68711027	18.81429517
A0JNW5_UHRF1BP1L	16.43142918	16.67144461	16.74756336	16.65585342	16.62092255	16.51066824	16.49140145	16.52488806	16.48749198
A0MZ66_SHTN1	17.83484479	17.94142353	17.89101706	17.86734878	17.8338081	17.61418106	17.68879521	17.74564679	17.66506562
A0PJW6_TMEM223	15.23967485	15.61231171	15.05802001	15.21593806	15.80265263	15.08599626	15.51236758	15.6720385	15.14184346
A0PJZ3_GXYLT2	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0PK00_TMEM120B	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1A4S6_ARHGAP10	17.47806268	17.51192954	17.55641363	17.56259037	17.54384565	17.2672408	17.35678869	17.56419388	17.57469676
A1L0T0_ILVBL	19.99394688	20.06810721	19.66496597	19.70985843	19.80396401	19.58046193	19.69443281	19.70985843	19.82289557
A1L157_TSPAN11	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1L188_NDUFAF8	17.13119594	17.22108511	17.25369198	17.58396373	17.31347717	17.23265906	16.9573077	16.99569748	16.90422991
A1L190_SYCE3	NA	NA	NA	NA	NA	NA	NA	NA	NA
```

```psiteData_prep``` prepares the psiteData by adding a fudge factor to the entire data and converting it to log2 scale. In addition, it normalizes the psite data by subtracting the abundance of mother proteins. So this step has to be done after ```protData_prep```.

```{r}
psiteData_p = psiteData_prep(psiteData_filename = "path/to/psiteData.tsv",
uniprot_gn_filename = "path/to/uniprot_gn.tsv",
fudge_facotr = 0.01,
sub_norm = T,
                          protData_p = protData_p)
```
Return value of this function is a data frame ```psiteData_p```

```
psite_names	p_0min_1	p_5min_1	p_30min_1	p_0min_6	p_5min_6	p_30min_6	p_0min_10	p_5min_10	p_30min_10
A0AVT1_UBA6_S743	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0FGR8_ESYT2_S676	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0FGR8_ESYT2_S755	-2.106807903	-2.071110428	-1.990699629	-2.386751678	-2.395851309	-2.043081626	-2.281712731	-2.383219491	-2.164456067
A0FGR8_ESYT2_S758	-2.516063422	-2.468364267	-2.3844267	-2.66919786	-2.710287156	-2.489419591	-2.583339325	-2.656893792	-2.478539962
A0FGR8_ESYT2_S761	-2.846134578	-2.772319577	-2.71229725	-3.066589113	-3.001674233	-2.84536355	-2.865064552	-2.914557646	-2.723566672
A0JNW5_UHRF1BP1L_S414	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0JNW5_UHRF1BP1L_S891	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0JNW5_UHRF1BP1L_S935	-0.303207769	-0.614282177	-0.370261673	-0.869499344	-0.669459336	-0.598223697	-0.365424571	-0.50891734	-0.350196927
A0JNW5_UHRF1BP1L_S953	0.755074794	0.535301282	0.686538485	0.509173855	0.763244158	0.509584231	0.649110284	0.459079414	0.670565016
A0JNW5_UHRF1BP1L_S989	NA	NA	NA	NA	NA	NA	NA	NA	NA
A0JNW5_UHRF1BP1L_S1021	-1.998013588	-1.997281685	-2.10717494	-2.469556906	-1.865628643	-2.733689605	-1.972000472	-1.723660498	-2.230431546
A0MZ66_SHTN1_S506	-2.705699347	-3.2426938	-2.933000854	-3.035424549	-3.144681898	-2.683032093	-2.807620333	-2.895777724	-2.26096696
A0MZ66_SHTN1_S534	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1A4S6_ARHGAP10_S376	NA	NA	NA	NA	NA	NA	NA	NA	NA
A1A4S6_ARHGAP10_S591	-1.076054244	-0.889605043	-0.937691185	-1.091125203	-1.177169259	-0.999275832	-1.039347165	-1.117905004	-1.165932895
A1A4S6_ARHGAP10_S600	-1.131683166	-1.105409807	-0.98637306	-1.273207046	-1.48468947	-1.305307775	-1.104743767	-1.365632859	-1.264955514
```



### 3.2 1D differential analysis 
With the prepared data we can proceed to differential analysis.


Assume we compare 5 min and 0 min levels in proteome and phosphoproteome. 

Column names for 5 min are
```coln_5min = c("p_5min_1", "p_5min_6", "p_5min_10",...)```

Column names for 0 min are 
```coln_0min = c("p_0min_1", "p_0min_6", "p_0min_10",...)```

Column names for 30 min are 
```coln_30min = c("p_30min_1", "p_30min_6", "p_30min_10",...)```


To make sure the permutation process is trackable, you can design seeds for the permutation in the following way:

```{r}

set.seed(123)
seeds30 = matrix(sample(c(1:30000),30000,replace = F), nrow = 1000, ncol = 30)

```

```{r}

all5_prot = comparison_time_points_1d( 
  s1_col_name = coln_5min,
  s2_col_name = coln_0min,
  d1_data = protData_p,
  nna_cutoff = 6,
  all_seeds = seeds30,
  permute_times = 500,
  working_dir = "/your/working/dir/",
  compare_name = "all_5min_prot")
 
```
The output files of this function are: 

* all_5min_prot_Z_hist.pdf
* all_5min_prot_z_null_hist.pdf
* all_5min_prot_Zz_dens_more.pdf
* all_5min_prot_Z_posterior.pdf
* all_5min_prot_ttest_check.pdf
* all_5min_prot_Z_result.tsv
* all_5min_prot_dens_overlay_prop.pdf


Similarly, in phosphoproteome,

```{r}

all5_psite = comparison_time_points_1d( 
  s1_col_name = coln_5min,
  s2_col_name = coln_0min,
  d1_data = psiteData_p,
  nna_cutoff = 6,
  all_seeds = seeds30,
  permute_times = 500,
  working_dir = "your/working/dir/",
  compare_name = "all_5min_psite")
 
```
The output files of this function are: 

* all_5min_psite_Z_hist.pdf
* all_5min_psite_z_null_hist.pdf
* all_5min_psite_Zz_dens_more.pdf
* all_5min_psite_Z_posterior.pdf
* all_5min_psite_ttest_check.pdf
* all_5min_psite_Z_result.tsv
* all_5min_psite_dens_overlay_prop.pdf


### 3.3 Mapping kinases with substrates
Before we move on to 2D analysis, it is necessary to map protein and phosphorylation data to known kinase-substrate relationships.

```{r}


ks_list = ks_map(protData_filename = "path/to/protData.tsv",
                  psiteData_filename = "path/to/psiteData.tsv",
                  uniprot_gn_filename = "path/to/uniprot_gn.tsv",
                  fudge_factor = 0.01,
                  sub_norm = T,
                  ksNetwork_filename = "path/to/ksNetwork.tsv",
                    working_dir = "your/working/dir/",
                  ks_outputName = "KS_network.tsv")

```

Output of this function is a file named ```KS_network.tsv```, and the return value is a list of two, including the kinase data as the first element, and substrate data as the second element. 


```
name	p_0min_1	p_5min_1	p_30min_1	p_0min_6	p_5min_6	p_30min_6	p_0min_10	p_5min_10	p_30min_10
AAK1_AP2M1_T156_PhosphoSitePlus_prot	19.05130882	19.13135092	19.18658006	19.13823954	18.96507767	18.84372563	19.01603163	19.08365103	18.99200926
AAK1_AP2M1_T156_PhosphoSitePlus_psite	14.61163911	15.55878063	14.64421324	14.59490085	14.5557962	14.45210284	14.12854395	14.38602224	14.48160983
AAK1_AP2M1_T156_PhosphoSitePlus_substrate_prot	22.15692699	22.11692572	21.99599792	22.11216524	22.26042201	22.04791699	21.97906042	21.91812022	21.92638182
AAK1_AP2M1_T156_PhosphoSitePlus_subProt_psite	-7.545287873	-6.558145096	-7.351784681	-7.517264385	-7.70462581	-7.595814155	-7.850516467	-7.532097985	-7.444771989
ABL2_RIN1_Y36_tyrosineKinase_prot	16.47833447	16.49705202	16.36567975	16.14234858	16.05738996	16.16929574	16.17951008	16.31182718	16.2982808
ABL2_RIN1_Y36_tyrosineKinase_psite	NA	NA	NA	NA	NA	NA	NA	NA	NA
ABL2_RIN1_Y36_tyrosineKinase_substrate_prot	18.42161841	18.71502854	18.51235203	18.2204172	18.37937183	18.10725307	18.31578051	18.46716489	18.41527253
ABL2_RIN1_Y36_tyrosineKinase_subProt_psite	NA	NA	NA	NA	NA	NA	NA	NA	NA
ACVR1_ETS1_S251_PhosphoNetworks_prot	17.57801903	17.73966127	17.64331657	17.51420114	17.61194496	17.270375	17.41023305	17.56045551	17.41273569
ACVR1_ETS1_S251_PhosphoNetworks_psite	NA	NA	NA	NA	NA	NA	NA	NA	NA
ACVR1_ETS1_S251_PhosphoNetworks_substrate_prot	17.09798672	17.31724257	17.02214181	17.12735977	17.12735977	16.82173061	17.30203841	16.99061415	17.17685877
ACVR1_ETS1_S251_PhosphoNetworks_subProt_psite	NA	NA	NA	NA	NA	NA	NA	NA	NA
ACVR1_NFATC4_S272_PhosphoNetworks_prot	17.57801903	17.73966127	17.64331657	17.51420114	17.61194496	17.270375	17.41023305	17.56045551	17.41273569
ACVR1_NFATC4_S272_PhosphoNetworks_psite	NA	NA	NA	NA	NA	NA	NA	NA	NA
ACVR1_NFATC4_S272_PhosphoNetworks_substrate_prot	15.46100079	15.9387852	15.66916252	15.63240535	15.97970268	15.67155058	15.3654388	15.73419395	15.82686008
ACVR1_NFATC4_S272_PhosphoNetworks_subProt_psite	NA	NA	NA	NA	NA	NA	NA	NA	NA
AKT1_AKT1S1_T246_PhosphoSitePlus_prot	19.6224442	19.81635347	19.41944128	19.51151207	19.86912028	19.3863779	19.60167954	19.73488315	19.75901951
AKT1_AKT1S1_T246_PhosphoSitePlus_psite	17.71323713	18.83610851	19.70877967	17.96851454	18.19641275	18.36637658	17.56581356	17.90817638	18.04732601
AKT1_AKT1S1_T246_PhosphoSitePlus_substrate_prot	17.23239747	17.42982722	17.19112489	17.56927542	17.38151639	17.22439251	17.15435039	17.07388354	17.00773044
AKT1_AKT1S1_T246_PhosphoSitePlus_subProt_psite	0.480839658	1.406281291	2.517654773	0.399239122	0.814896359	1.141984066	0.41146317	0.834292837	1.039595569
AKT1_ARHGAP22_S16_PhosphoSitePlus_prot	19.6224442	19.81635347	19.41944128	19.51151207	19.86912028	19.3863779	19.60167954	19.73488315	19.75901951
AKT1_ARHGAP22_S16_PhosphoSitePlus_psite	17.87086373	18.01017114	18.0259286	18.21141886	18.20403123	18.16927069	18.0538812	17.61942734	18.0406973
AKT1_ARHGAP22_S16_PhosphoSitePlus_substrate_prot	18.56895248	18.62575146	18.43487308	18.59643362	18.84813564	18.72727272	18.61985413	18.53395605	18.55760479
AKT1_ARHGAP22_S16_PhosphoSitePlus_subProt_psite	-0.698088742	-0.615580316	-0.408944479	-0.385014755	-0.64410441	-0.558002021	-0.565972931	-0.914528709	-0.516907489
AKT1_ATXN1_S775_PhosphoSitePlus_prot	19.6224442	19.81635347	19.41944128	19.51151207	19.86912028	19.3863779	19.60167954	19.73488315	19.75901951
AKT1_ATXN1_S775_PhosphoSitePlus_psite	15.82296078	15.87178627	16.20170683	15.8145091	15.62961163	15.59509085	15.7826247	16.07541037	15.88724731
AKT1_ATXN1_S775_PhosphoSitePlus_substrate_prot	18.11111298	18.28531205	18.03266425	17.86187222	18.00165335	17.75435646	17.98799349	17.9858943	18.2666545
AKT1_ATXN1_S775_PhosphoSitePlus_subProt_psite	-2.288152206	-2.413525774	-1.830957418	-2.047363122	-2.372041717	-2.159265617	-2.205368795	-1.91048393	-2.379407196
AKT1_BABAM1_S29_PhosphoSitePlus_prot	19.6224442	19.81635347	19.41944128	19.51151207	19.86912028	19.3863779	19.60167954	19.73488315	19.75901951
AKT1_BABAM1_S29_PhosphoSitePlus_psite	15.49160485	16.16029459	15.82806277	15.54591313	16.48859232	15.79062292	15.28836727	15.97339172	15.65708376
AKT1_BABAM1_S29_PhosphoSitePlus_substrate_prot	17.84847883	17.95978166	17.66903008	18.02531894	17.72434754	17.74736796	17.98058083	17.89913249	17.87815506
AKT1_BABAM1_S29_PhosphoSitePlus_subProt_psite	-2.356873973	-1.799487069	-1.840967312	-2.479405812	-1.235755221	-1.956745036	-2.692213559	-1.925740775	-2.221071298
```


### 3.4 2D differential analysis 

With prepared kinases and substrates data, suppose we want to identify joint variation in proteome at 5 min and phosphoproteome at 30 min.
We can conduct a 2D analysis using the following function:


```{r}




comparison_time_points_2dZ_5_30(s1_d1_col_name = coln_5min,
  s2_d1_col_name = coln_0min,
  s1_d2_col_name = coln_30min,
  s2_d2_col_name = coln_0min,
  d1_data = ks_list[[1]],
  d2_data = ks_list[[2]],
  nna_cutoff = 6,
  all_seeds = seeds30,
  permute_times = 500,
  working_dir = "your/working/dir/",
  compare_name = "all_5prot_30psite")


```



Output files are:
* all_5prot_30psite_2d_Z_dens_function.pdf
* all_5prot_30psite_2d_z_null_dens_function.pdf
* all_5prot_30psite_2d_Z_posterior.pdf
* all_5prot_30psite_Z_result.tsv
* all_5prot_30psite_2d_f1_function.pdf


## 4 Tips

You can find a sample script of running the entire sequence of functions in this github page "KSA2D_sample_script.R"

