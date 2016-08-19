

### **`multigraph`**: Plot and Manipulate Multigraphs in R

<br />

When you load the **`multigraph`** package then **`multiplex`** is automatically invoked

```{r }
library("multigraph")
Loading required package: multiplex
```

<br />

### Example: Florentine Families data set

We work with Padgett's Florentine Families data set, which is publicly available as a Ucinet DL file format. We use function `read.dl` of the **`multiplex`** package to retrieve this data with the **R** console.

<br />


```{r }
### Read the Padgett Florentine Families data set as a Ucinet DL file
### from a public repository and storage it as an object
floflies <- read.dl(file="http://moreno.ss.uci.edu/padgett.dat")


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


`floflies` represents this network where `"PADGM"` are marriage relations, and `"PADGB"` correspond to business ties among 16 Florentine families.


<br /> 


## Ploting the Multigraph

Now we plot this network with the `multigraph` function:

```{r }
multigraph(floflies)
```
![Default layout of `floflies`](figs/floflies.png)


<br /> 

The network is symmetric and the default layout is circular. Check also out the [vector image](figs/floflies.pdf) of this multigraph.

<br /> 


#### Stress majorization layout

Besides the circular layout, another possibility is to apply a force-directed layout for the visualization, and below the Florentine Families data set is depicted with the stress majorization algorithm with a number of arguments. 

<br /> 


```{r }
multigraph(floflies, directed = FALSE, layout = "stress", seed = 2, cex = 6, tcex = .7, pos = 0, vcol = 8,
           ecol = 1, lwd = 2, bwd = .5)
```
![stress layout of `floflies`](figs/floflies-stress.png)
[vector image](figs/floflies-stress.pdf)

<br />
<br /> 


```{r }
multigraph(floflies, directed = FALSE, layout = "stress", seed = 2, cex = 6, tcex = .7, pos = 0, vcol = 8,
           ecol = 1, lwd = 2, bwd = .5, lty = 2:1, pch = 13)
```
![stress layout of `floflies` different shapes](figs/floflies-stress2.png)
[vector image](figs/floflies-stress2.pdf)

<br />

Note that when the graph is depicted as *undirected* then the reciprocal ties by default are collapsed, and you can prevent this to happen by setting the argument `collRecip`  to  `FALSE`. Some arguments such as `cex`, `lwd`, `lty`, `pch` are graphical parameters of the **`graphics`** package, and other arguments like `bwd` to specify the width of the bundle type, `tcex` for the size of the vertices labels, or `ecol` and `vcol` for the color of edges and vertices to set the shape and color of both the vertices and the edges in the multigraph are complementary. Moreover, by setting the `pos` argument to `0`, the actor labels are placed in the middle of the vertices.

<br /> 
<br /> 


## Ploting with Actor Attributes


In a similar way than before, we obtain some actor attributes of the Florentine Families network with the `read.dl` function from this repository.


```{r }
flofliesatt <- read.dl(file="http://moreno.ss.uci.edu/padgw.dat")
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
multigraph(floflies, directed = FALSE, layout = "stress", seed = 1, cex = flofliesatt[,1],
           tcex = .8, pos = 0, lwd = 2, ecol = 1, vcol = 5)
```
![stress layout of `floflies`](figs/flofliesatt-stress.png)
[vector image](figs/flofliesatt-stress.pdf)


<br /> 

And with the `clu` argument we establish the clustering of the network with three classes of actors differentiated by the colors of the vertices.


```{r }
multigraph(floflies, directed = FALSE, layout = "stress", seed = 1, cex = flofliesatt[,1], tcex = .8, 
           pos = 0, lwd = 2, ecol = "white", vcol = c("orange", "blue", "white"), alpha = c(.5, 1, .2), 
           clu = c(1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 2, 3, 1, 1, 2, 1), bg = 1)
```
![stress layout of `floflies` with clustering](figs/flofliesatt-stress2.png)
[vector image](figs/flofliesatt-stress2.pdf)


*Beware that alpha transparency renders differently according to the device used by* **`grDevices`**



<br /> 
<br /> 


### Bipartite Graphs

TBD




