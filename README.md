[![Build Status](https://travis-ci.org/mplex/multigraph.svg?branch=master)](https://travis-ci.org/mplex/multigraph)
[![CRAN version](http://www.r-pkg.org/badges/version/multigraph)](https://cran.r-project.org/package=multigraph)


<br />

### **`multigraph`**: Plot and Manipulate Multigraphs in R
#### Author: Antonio Rivero Ostoic (@mplex)

<br />

Install **`multigraph`** 

```{r }
### from CRAN
install.packages("multigraph")
```

or

```{r }
### from Github
devtools::install_github("mplex/multigraph")
```

<br />


When you load the package then **`multiplex`** is automatically invoked.

```{r }
library("multigraph")
# Loading required package: multiplex
```

<br />

### Multigraph: Florentine Families data set

We work with Padgett's Florentine Families data set, which is publicly available as a Ucinet DL file format. We use function `read.dl` of the **`multiplex`** package to retrieve this data with the **R** console.

<br />


```{r }
### Read the Padgett Florentine Families data set as a Ucinet DL file
### from a public repository and storage it as an object

floflies <- read.dl(file = "http://moreno.ss.uci.edu/padgett.dat")


### take a look at this data
floflies

, , PADGM

          ACCIAIUOL ALBIZZI BARBADORI BISCHERI CASTELLAN GINORI GUADAGNI LAMBERTES MEDICI PAZZI PERUZZI PUCCI RIDOLFI SALVIATI STROZZI TORNABUON
ACCIAIUOL         0       0         0        0         0      0        0         0      1     0       0     0       0        0       0         0
ALBIZZI           0       0         0        0         0      1        1         0      1     0       0     0       0        0       0         0
BARBADORI         0       0         0        0         1      0        0         0      1     0       0     0       0        0       0         0
BISCHERI          0       0         0        0         0      0        1         0      0     0       1     0       0        0       1         0
CASTELLAN         0       0         1        0         0      0        0         0      0     0       1     0       0        0       1         0
GINORI            0       1         0        0         0      0        0         0      0     0       0     0       0        0       0         0
GUADAGNI          0       1         0        1         0      0        0         1      0     0       0     0       0        0       0         1
LAMBERTES         0       0         0        0         0      0        1         0      0     0       0     0       0        0       0         0
MEDICI            1       1         1        0         0      0        0         0      0     0       0     0       1        1       0         1
PAZZI             0       0         0        0         0      0        0         0      0     0       0     0       0        1       0         0
PERUZZI           0       0         0        1         1      0        0         0      0     0       0     0       0        0       1         0
PUCCI             0       0         0        0         0      0        0         0      0     0       0     0       0        0       0         0
RIDOLFI           0       0         0        0         0      0        0         0      1     0       0     0       0        0       1         1
SALVIATI          0       0         0        0         0      0        0         0      1     1       0     0       0        0       0         0
STROZZI           0       0         0        1         1      0        0         0      0     0       1     0       1        0       0         0
TORNABUON         0       0         0        0         0      0        1         0      1     0       0     0       1        0       0         0

, , PADGB

          ACCIAIUOL ALBIZZI BARBADORI BISCHERI CASTELLAN GINORI GUADAGNI LAMBERTES MEDICI PAZZI PERUZZI PUCCI RIDOLFI SALVIATI STROZZI TORNABUON
ACCIAIUOL         0       0         0        0         0      0        0         0      0     0       0     0       0        0       0         0
ALBIZZI           0       0         0        0         0      0        0         0      0     0       0     0       0        0       0         0
BARBADORI         0       0         0        0         1      1        0         0      1     0       1     0       0        0       0         0
BISCHERI          0       0         0        0         0      0        1         1      0     0       1     0       0        0       0         0
CASTELLAN         0       0         1        0         0      0        0         1      0     0       1     0       0        0       0         0
GINORI            0       0         1        0         0      0        0         0      1     0       0     0       0        0       0         0
GUADAGNI          0       0         0        1         0      0        0         1      0     0       0     0       0        0       0         0
LAMBERTES         0       0         0        1         1      0        1         0      0     0       1     0       0        0       0         0
MEDICI            0       0         1        0         0      1        0         0      0     1       0     0       0        1       0         1
PAZZI             0       0         0        0         0      0        0         0      1     0       0     0       0        0       0         0
PERUZZI           0       0         1        1         1      0        0         1      0     0       0     0       0        0       0         0
PUCCI             0       0         0        0         0      0        0         0      0     0       0     0       0        0       0         0
RIDOLFI           0       0         0        0         0      0        0         0      0     0       0     0       0        0       0         0
SALVIATI          0       0         0        0         0      0        0         0      1     0       0     0       0        0       0         0
STROZZI           0       0         0        0         0      0        0         0      0     0       0     0       0        0       0         0
TORNABUON         0       0         0        0         0      0        0         0      1     0       0     0       0        0       0         0
```

<br /> 

Object `floflies` represents this network where `"PADGM"` are marriage relations, and `"PADGB"` correspond to business ties among 16 Florentine families.


<br /> 


### Ploting the Multigraph

We plot this network with the `multigraph` function:

```{r }
multigraph(floflies)
```
![Default layout of `multigraph`](figs/floflies.png)


<br /> 

The network is symmetric and the default layout is circular. Check also out the [vector image](figs/floflies.pdf) of this multigraph.

<br /> 


#### Force directed layout

Besides the circular layout, another possibility is to apply a force-directed layout for the visualization, and below the Florentine Families data set is depicted with the force directed algorithm with a number of arguments. 

<br /> 


```{r }
multigraph(floflies, directed = FALSE, layout = "force", seed = 2, cex = 6, tcex = .7, pos = 0, vcol = 8,
+  ecol = 1, lwd = 2, bwd = .5)
```
![Force directed layout of `multigraph`](figs/floflies-force.png)
[vector image](figs/floflies-force.pdf)

<br />
<br /> 


```{r }
multigraph(floflies, directed = FALSE, layout = "force", seed = 2, cex = 6, tcex = .7, pos = 0, vcol = 8,
+  ecol = 1, lwd = 2, bwd = .5, lty = 2:1, pch = 13)
```
![Force directed layout of `multigraph` different shapes](figs/floflies-force2.png)
[vector image](figs/floflies-force2.pdf)

<br />

Note that when the graph is depicted as *undirected* then the reciprocal ties by default are collapsed, and you can prevent this to happen by setting the argument `collRecip`  to  `FALSE`. Some arguments such as `cex`, `lwd`, `lty`, `pch` are graphical parameters of the **`graphics`** package to set the shape of both the vertices and the edges, whereas other arguments like `bwd` to specify the width of the bundle type, `tcex` for the size of the node labels, or `ecol` and `vcol` for the color of edges and vertices respectively are complementary in **`multigraph`**. Moreover, by setting the `pos` argument to `0`, the actor labels are placed in the middle of the nodes.


<br /> 
<br /> 


## Ploting with Actor Attributes


In a similar way than before, we obtain some actor attributes of the Florentine Families network with the `read.dl` function from this repository.


```{r }
flofliesatt <- read.dl(file = "http://moreno.ss.uci.edu/padgw.dat")
```

and we take a look at the ` flofliesatt` that storages this type of information


```{r }
flofliesatt

          WEALTH #PRIORS #TIES
ACCIAIUOL     10      53     2
ALBIZZI       36      65     3
RIDOLFI       27      38     4
STROZZI      146      74    29
BARBADORI     55       0    14
BISCHERI      44      12     9
CASTELLAN     20      22    18
GUADAGNI       8      21    14
LAMBERTES     42       0    14
MEDICI       103      53    54
PAZZI         48       0     7
PERUZZI       49      42    32
SALVIATI      10      35     5
TORNABUON     48       0     7
GINORI        32       0     9
PUCCI          3       0     1
```


However, in order to depict the multigraph of `floflies` with the information contained in `flofliesatt`, we need to be sure that the order of the actors matches in both objects.

```{r }
flofliesatt <- flofliesatt[order(rownames(flofliesatt)), ]

flofliesatt

          WEALTH #PRIORS #TIES
ACCIAIUOL     10      53     2
ALBIZZI       36      65     3
BARBADORI     55       0    14
BISCHERI      44      12     9
CASTELLAN     20      22    18
GINORI        32       0     9
GUADAGNI       8      21    14
LAMBERTES     42       0    14
MEDICI       103      53    54
PAZZI         48       0     7
PERUZZI       49      42    32
PUCCI          3       0     1
RIDOLFI       27      38     4
SALVIATI      10      35     5
STROZZI      146      74    29
TORNABUON     48       0     7
```

That's much better now.


<br /> 
<br /> 

The following code serves to depict this network in a way that the size of the vertices reflects the wealth of the actors.


```{r }
multigraph(floflies, directed = FALSE, layout = "force", seed = 1, cex = flofliesatt[,1],
+  tcex = .8, pos = 0, lwd = 2, ecol = 1, vcol = 5)
```
![Force directed layout layout of `multigraph` with attributes](figs/flofliesatt-force.png)
[vector image](figs/flofliesatt-force.pdf)


<br /> 

And with the `clu` argument we establish the clustering of the network with three classes of actors differentiated by the colors of the vertices.


```{r }
multigraph(floflies, directed = FALSE, layout = "force", seed = 1, cex = flofliesatt[,1], tcex = .8, 
+  pos = 0, lwd = 2, ecol = "white", vcol = c("orange", "blue", "white"), alpha = c(.5, 1, .2), 
+  clu = c(1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 3, 1, 1, 2, 1), bg = 1)
```
![Force directed layout of `multigraph` with clustering](figs/flofliesatt-force2.png)
[vector image](figs/flofliesatt-force2.pdf)

Hence, colors can be established in different ways, and the `alpha` vector argument serves to set the transparency of vertices, edges, and background colors respectively.




<br /> 
===
<br /> 


### Bipartite Graph: Southern Women data set

Support for the visualization of two-mode networks is also given by **`multigraph`** and we work with the Southern Women classic data set to illustrate some of the layout options with this package.

```{r }
### Read the Ucinet DL file of Davis, Gardner, Gardner Southern Women
### data set from a public repository and storage it as an object

swomen <- read.dl(file = "http://moreno.ss.uci.edu/davis.dat")

### take a look...
swomen

          E1 E2 E3 E4 E5 E6 E7 E8 E9 E10 E11 E12 E13 E14
EVELYN     1  1  1  1  1  1  0  1  1   0   0   0   0   0
LAURA      1  1  1  0  1  1  1  1  0   0   0   0   0   0
THERESA    0  1  1  1  1  1  1  1  1   0   0   0   0   0
BRENDA     1  0  1  1  1  1  1  1  0   0   0   0   0   0
CHARLOTTE  0  0  1  1  1  0  1  0  0   0   0   0   0   0
FRANCES    0  0  1  0  1  1  0  1  0   0   0   0   0   0
ELEANOR    0  0  0  0  1  1  1  1  0   0   0   0   0   0
PEARL      0  0  0  0  0  1  0  1  1   0   0   0   0   0
RUTH       0  0  0  0  1  0  1  1  1   0   0   0   0   0
VERNE      0  0  0  0  0  0  1  1  1   0   0   1   0   0
MYRA       0  0  0  0  0  0  0  1  1   1   0   1   0   0
KATHERINE  0  0  0  0  0  0  0  1  1   1   0   1   1   1
SYLVIA     0  0  0  0  0  0  1  1  1   1   0   1   1   1
NORA       0  0  0  0  0  1  1  0  1   1   1   1   1   1
HELEN      0  0  0  0  0  0  1  1  0   1   1   1   1   1
DOROTHY    0  0  0  0  0  0  0  1  1   1   0   1   0   0
OLIVIA     0  0  0  0  0  0  0  0  1   0   1   0   0   0
FLORA      0  0  0  0  0  0  0  0  1   0   1   0   0   0
```

In this case the information can be contained in a data frame or an array as before.


<br/>


#### Plotting two-mode data

Function `bmgraph` serves to plot two-mode data or an affiliation network as a bipartite graph.

```{r }
bmgraph(swomen)
```
![Default layout of `bmgraph`](figs/swomen.png)
[vector image](figs/swomen.pdf)


In this case actor and events have different shape by default.


<br />

Similarly to `multigraph` the color and shape of edges and vertices can be modified by equal arguments, and we can mirror the *x* axis of the plot.

```{r }
bmgraph(swomen, cex = 3, tcex = .8, pch = c(19, 15), lwd = 1.5, vcol = 2:3, mirrorX = TRUE)
```
![Mirror X of `bmgraph`](figs/swomen2.png)
[vector image](figs/swomen2.pdf)


<br />

Option `bip3` splits the actors in two columns, whereas `bip3e` will split the events.


```{r }
bmgraph(swomen, layout = "bip3", cex = 3, tcex = .8, pch = c(19, 15), lwd = 1.5, vcol = 2:3)
```
![Mirror X of `bmgraph`](figs/swomen3.png)
[vector image](figs/swomen3.pdf)


<br />

The binomial projection of a two-mode data set allows obtaining a force directed layout that in this case the image is clockwise rotated 65 degrees.

```{r }
bmgraph(swomen, layout = "force", seed = 1, cex = 3, tcex = .8, pch = c(19, 15), lwd = 2, vcol = 2:3, 
+  ecol = 8, rot = 65)
```
![Force directed layout of `bmgraph`](figs/swomen-force.png)
[vector image](figs/swomen-force.pdf)


<br />
<br />

Finally, function `bmgraph` stands for a bipartite *multigraph*, and this is because the actors can be affiliated at different levels.


```{r }
bmgraph(floflies, ecol = 1)
```
![bipartite graph of `floflies`](figs/floflies-bmgraph.png)
[vector image](figs/floflies-bmgraph.pdf)


<br />
<br />


*Note that the rendering of vector images may vary according to the device used.*

