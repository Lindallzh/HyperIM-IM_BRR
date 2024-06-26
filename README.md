# HyperIM
This project implements the HyperIM and HyperIM_BRR algorithms.
The code is developed on the basis of the project: https://github.com/tangj90/OPIM & https://github.com/qtguo/subsim

Author: Lingling Zhang 

## Compile

```shell
make
```
Before executing, please decompress ‘dataset.rar’
The option -std=c++1x should be included.

## Options

- **-func**=string 

  Specify the function to execute. Possible values are in the following:

  - **format**: format the graph
  - **im**: influence maximization

- **-gname**=string

  Specify the name of the input graph. 

- **-dir**=string

  Specify the directory of the input graphs. The default directory is  **graphInfo** in the current directory.

- **-outpath**=string

  Specify the directory for saving the output. The default directory is **result** in the working directory

**The following options are used with -func=format.**

- **-pdist**=string

  Specify the probability distribution. 

  - **wc** [default]: WC setting (default). The weight of the edge (u, v) is set to p(u, v) = 1/indeg(v).
  - **uniform**: Uniform setting. All the edges are assigned a given weight specified by **-pedge**=float.
  - **skewed**:  Skewed distribution. The weights follow a skewed distribution.

- **-wcvariant**=float

  Specify the weight sum in the wc variant setting. For example, if **-wcvariant**=1.2, it means the weight of the edge (u, v) is set as min{1, 1.2/indegree(v)}.  This option works for **-pdist=wc** only. 

- **-pedge**=float

  Specify the probability for all the edges in the uniform setting. The default value is **0.1**. This option works for **-pdist=uniform** only.

- **-skew**=string

  Specify which skewed distribution is used. This option works for **-pdist=skewed** only.

  -  **exp**[default]: Exponential distribution.
  - **weibull**: Weibull distribution.

**The following options are used with -func=im.**

- **-seedsize**=integer

  Specify the size of the seed set. The default value is **50**.

- **-eps**=float

  Specify the error threshold such that subsim returns **(1-1/e-eps)**-approximate solution. The default value is 0.1.

- **-delta**=float

  Specify the failure probability. It means that the returned solution is valid with at least **1 - delta** probability. The default value is 1/n.

- **-vanilla**=integer

  Specify the RR set generation method. 

  - **0** [default]: the RR set is generated by the subsim method.
  - **1**: the RR set is generated by the vanilla version. 
  
- **-hist**=integer

  Specify whether the HIST algorithm is invoked.

  - **0** [default]: the HIST algorithm is not invoked.
  - **1**: the HIST algorithm is invoked.

## Format the graph

**Before running influence maximization algorithms, please format the graph first.**

Examples:

```shell
//the wc setting and wc variant setting
./subsim -func=format -pdist=wc --dataname  NDC-substances-full --inputpath dataset/
./subsim -func=format -pdist=wc --dataname  tags-ask-ubuntu --inputpath dataset/
./subsim -func=format -pdist=wc --dataname  threads-ask-ubuntu --inputpath dataset/
./subsim -func=format -pdist=wc --dataname  email-Eu-full --inputpath dataset/
./subsim -func=format -pdist=wc --dataname  coauth-MAG-Geology-full --inputpath dataset/
./subsim -func=format -pdist=wc --dataname  walmart --inputpath dataset/

./subsim -func=format -pdist=wc --dataname  tags-math-sx --inputpath dataset/

tags-math-sx

./subsim -func=format -pdist=wc --dataname  email-Enron-full --inputpath dataset/
./subsim -func=format -pdist=wc --dataname  stack-overflow --inputpath dataset/ 
$ ./subsim -func=format -gname=facebook -pdist=wc
$ ./subsim -func=format -gname=facebook -pdist=wc -wcvariant=1.2

//the uniform setting
$ ./subsim -func=format -gname=NDC-substances-full -pdist=uniform -pedge=0.01

./subsim -func=format -gname=facebook -pdist=uniform -pedge=0.01

//the skewed setting with exponential or weibull distribution
$ ./subsim -func=format -gname=facebook -pdist=skewed -skew=exp
$ ./subsim -func=format -gname=facebook -pdist=skewed -skew=weibull

./subsim -func=format -pdist=uniform -pedge=0.01 --dataname  NDC-substances-full --inputpath dataset/
```

## Influence maximization

Examples:

```shell
//use the subsim method to sample RR sets
$ ./subsim -func=im -gname=facebook -seedsize=100 -eps=0.01

./subsim -func=im -seedsize=300 -eps=0.01 --inputpath dataset/ --dataname  NDC-substances-full 
./subsim -func=im -seedsize=100 -eps=0.01 --inputpath dataset/ --dataname  contact-high-school
./subsim -func=im -seedsize=100 -eps=0.01 --inputpath dataset/ --dataname  coauth-MAG-Geology-full
./subsim -func=im -seedsize=100 -eps=0.01 --inputpath dataset/ --dataname  threads-ask-ubuntu
./subsim -func=im -seedsize=100 -eps=0.01 --inputpath dataset/ --dataname  email-Eu-full
./subsim -func=im -seedsize=200 -eps=0.01 --inputpath dataset/ --dataname  tags-ask-ubuntu
./subsim -func=im -seedsize=100 -eps=0.01 --inputpath dataset/ --dataname  walmart
./subsim -func=format -pdist=wc --dataname  walmart --inputpath dataset/
./subsim -func=im -seedsize=400 -eps=0.01 --inputpath dataset/ --dataname  threads-ask-ubuntu

./subsim -func=im -seedsize=300 -eps=0.01 --inputpath dataset/ --dataname  coauth-MAG-Geology-full

./subsim -func=format -pdist=wc --dataname  email-Eu-full --inputpath dataset/
./subsim -func=im -seedsize=100 -eps=0.01 --inputpath dataset/ --dataname  stack-overflow

./subsim -func=format -pdist=wc --dataname  catCooking_he --inputpath dataset/
./subsim -func=im -seedsize=100 -eps=0.01 --inputpath dataset/ --dataname  catCooking_he
#./subsim -func=im -pdist=wc --dataname  NDC-substances-full --inputpath dataset/ stack-overflow

//use the vanilla method to sample RR sets
$ ./subsim -func=im -gname=facebook -seedsize=100 -eps=0.01 -vanilla=1
```

