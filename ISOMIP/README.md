
# Isomip demonstration case
Here a description of the test case + link toward the src and notebooks. 
<br>
We here provide a physical description of this experiment and additional details as to how to run this experiment within NEMO. This experiment is **created and tested** for NEMO **code at revision 8097**. 

A **ipython notebook is also provided** as a demonstration of possible analysis. If you have already run the NEMO experiment and want to analyse the resulting output, you can directly look at the notebook : **[here](....)**.

## Objectives
The ..... experiment is......


## Physical description
NEMO .... demonstration case follows the specifications of .... et al. (2012): 
......  <br>


### Exemple of run
In this exemple we assess the behaviour of the ....<br>

* The **Reference Simulation** : **FCT2** is the first simulation, ....

```
cd TEST_CASES/........
ln -sf namelist_....._cfg namelist_cfg
```
choice of..... is done in namtra_.... block of namelist: 

~~~fortran
!-----------------------------------------------------------------------

~~~

Run the executable : (if you haven't compiled NEMO see [here](https://github.com/sflavoni/NEMO-test-cases) )

``` 
 mpirun -np 1 ./opa 
```
Output files are: <br>

~~~

~~~

* The **Sensibility Simulation** : .


```
ln -sf namelist_xxxx_cfg namelist_cfg
```

Run the executable again : 

``` 
 mpirun -np 1 ./opa 
```

Output files : <br>

~~~

~~~

* You can change output file name  in variable @expname@ in file\_def\_nemo-opa.xml

~~~xml
<file_definition type="multiple_file" name="@expname@" sync_freq="10d" min_digits="4">
~~~

* Available notebook python is **[here]()**.

## References
