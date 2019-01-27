
# NEMO Demonstration Cases
<!---
COMMENT ON STYLE: 
this is useful sintax to have nice grey box
[comment]: <

 ```
cd TEST\_CASES\_NEMO
mkdir MY\_TEST 
cd MY\_TEST 
svn --username 'mylogin' co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM NEMOGCM_r8097 -r 
8097
```
>
-->

<!--
equivalent commands:
<pre> </pre>
or
````

-->

<!--
comment:
can choise:
    type of font 
    dimension in points
    and color

<code>
<span style="color:red; font-family: 'Courier'; font-size: 18pt;"> memo for style
 : prova</span>
</code>
-->

This repository contains informations on :

1. How to use demonstration cases within NEMO and how to analyse their ouputs :
<br> demonstration cases are ready to be used, all files needed are availables within NEMO.

2. How to contribute and add new demontration case :
<br> you can follow the roadmap to add your own demontration case 

Demonstration cases available in NEMO repository  are :

- LOCK_EXCHANGE
- OVERFLOW
- ISOMIP
- WAD (Wetting & Dry)
- SAS_BIPER
- VORTEX
- CANAL

### How to use demonstration cases in NEMO 

* Download & compile **XIOS** code : 

<pre>
mkdir ~/XIOS; cd ~/XIOS
svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0
cd xios-2.0
./make_xios --arch your-compiler --jobs 8

( ./make_xios --help  (to choice your compiler) )
</pre>

*  Download & compile **NEMO** code

<pre>
mkdir my_TEST 
cd my_TEST 
svn --username 'mylogin' co http://forge.ipsl.jussieu.fr/nemo/svn/trunk/NEMOGCM NEMOGCM
cd NEMOGCM/CONFIG
./makenemo -a TEST_CASES -n <i>name_of_test_case</i> -m <i>your_ARCH_FILE</i>

( example ./makenemo -a TEST_CASES -n OVERFLOW -m X64_ADA )
</pre>

### Set of runs
If you want to run one of these test cases you can read README instructions here :


[LOCK_EXCHANGE](https://github.com/sflavoni/NEMO-test-cases/blob/master/lock-exchange/README.md)

[OVERFLOW](https://github.com/sflavoni/NEMO-test-cases/blob/master/overflow/README.md)

[ISOMIP](https://github.com/sflavoni/NEMO-test-cases/blob/master/isomip/README.md)

[WAD (Wetting & Dry)](https://github.com/sflavoni/NEMO-test-cases/blob/master/wetting_and_drying/README.md)

[SAS_BIPER](https://github.com/sflavoni/NEMO-test-cases/blob/master/sas_biperiodic/README.md)

[VORTEX](https://github.com/sflavoni/NEMO-test-cases/blob/master/vortex/README.md)

[CANAL](https://github.com/sflavoni/NEMO-test-cases/blob/master/canal/README.md)


### How to contribute and add new demontration case 

TO BE DONE......
