Above of all, thanks for your grateful help to run my codes in the computer cluster! I will tell you some simple but necessary items to complete jobs successfully.

-   There should be [3]{style="color:red"} RMD files in HCP codes, including "data_process.Rmd", "activation_brain.Rmd" and "summary_brain.Rmd". Please [override]{style="color:red"} all RMD files!

-   Meanwhile, you will receive the HCP data. There is a compressed file called "WM_contrasts" in the directory. Please unzip it initially.

-   Before running codes, please modify the corresponding file path to data in [all]{style="color:red"} RMD files. The correct file path should be the data directory, not the root directory, with an example ".../HCP data". They are in the top block and you can see them easily.

-   If you want to change the loop in different subjects and submit multiple tasks, you can modify the variable value *nsub* and make the following modifications:

```         
nsub_1 = ?; nsub_2 = ?
for idx in c(1:nsub)   ------>   for idx in c(nsub_1:nsub_2)
```

-   The order to run codes is: "data_process" ---\> "activation_brain" ---\> "summary_brain". [Before]{style="color:red"} running "summary_brain.Rmd", be sure that all result files have been collected if you submit multiple tasks.

-   Finally, just return the summary csv files to me. They are in the "HCP data/tables" sub-dir and called "summary_brain.csv".
