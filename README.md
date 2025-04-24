# PCEGS: <br><small>Identification of protein complexes through integrative gene ontology parent-child semantic similarity and dynamic protein-protein interactions</small>

## Table of contents
* [Features](#Features)
* [Algorithm](#Process)
* [Method](#Method)
* [Experiment](#Experiment)
* [Email](#Email)
* [Acknowledgments](#Acknowledgments)

## Features
PCEGS fuses GO semantic similarity and dynamic networks to identify protein complexes from protein interaction networks.GO semantic similarity can also be used on its own to compute semantic similarity between computed GO terms, or to compute functional similarity of gene products to weight protein interaction networks.

## Algorithm
Despite extensive research on identifying protein complexes from PINs, challenges remain in effectively reducing the impact of noisy data and constructing more meaningful dynamic protein interaction networks(DPINs). To this end, we propose a protein complex recognition algorithm (PCEGS) that innovatively integrates gene expression data, key protein data, and protein-protein interaction data, and sets unique active expression thresholds for each protein, which significantly improves the accuracy and recall of protein complex identification. The article also extends gene ontology semantic similarity and protein semantic similarity, and introduces a network homogeneity index to assign weights to the PIN, which enhances the algorithm's ability to capture complex relationships in organisms. In the DPINs, PCEGS employs a core-attachment structure to recognize protein complexes and uses a cohesion index to decide whether to accept the recognition results. Experimental results show that PCEGS achieves better performance in multiple indicators, such as ACC, Recall and F-score.
## Method

#### Compiling
The PCEGS algorithm is written in Rust code. It can be compiled using the latest version of Rust. The compilation command is
```
cargo build -r
```
#### Runing

```
cargo run --bin -r pcegs
```

## Experiment

The evaluation metrics code, written in python 2.7, is a community-implemented, general-purpose evaluation metric for protein complex identification algorithms.It is placed in the evaluation folder.

The data used were placed in the data folder, including GO annotations, essential proteins, PPIs, gene expression profiles, and protein complexes. The results of the experiments are placed in the results folder by default.

- the number of predicted complexes
- Fscore: the harmonic mean of Recall and Precision
- Precision: the rate of predicted protein complexes that match at least one reference complex 
- Recall: the rate of reference protein complexes that match at least one predicted complex
- ACC: the geometric accuracy
- Sn: the rate of the maximum-sum number of matched proteins to the total number of proteins in the set of reference protein complex.
- PPV: the rate of the maximumsum number of matched proteins to the total matched number of proteins in the set of predicted protein complex. 

To evaluate the "Predicted Complexes.txt" file, it runs as follows:
The match is coded in python. You must use python 2.7 of the interpreter to get the program.
			
python match.py -n ppi_name.txt reference_name.txt result_name.txt
            

The parameters of these methods are set as the recommended values as mentioned in their original papers. For our method, we set the balance index β to 1.5 as to recommended values.

| Datasets         | \#predicated | F\-score | ACC    |
|------------------|--------------|----------|--------|
| Krogan\-core     | 704          | 0\.558   | 0\.528 |
| Krogan\-extended | 778          | 0\.538   | 0\.482 |
| Gavin            | 832          | 0\.668   | 0\.560 |
| Collins          | 794          | 0\.614   | 0\.550 |

The experimental performance proves that the BOPS algorithm can obtain the best results on F-score and ACC when identifying small protein complexes. And the performance of BOPS is better than most of algorithms with respect to the whole protein complexes.

## Email
If you have any question, you can contact with jh.petrichor@gmail.com
## Acknowledgments
Supported by Shandong Province Education Quality Improvement Plan for Graduate Students (Grant No. SDYJG21211) and Undergraduate Teaching Reform Research Program of Shandong Province (Grant No. Z2022036).

At last, thanks very much for using PCEGS.