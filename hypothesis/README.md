## Running all Scripts
The core functions are implemented in the hypothesis module, which can be run with the following command from the
 root of the repository:
```commandline
(test_env) C:\\...\\VHL-Scripts>py -m hypothesis
```

Without any commandline arguments, the scripts will fetch all annotations from the  VHL Hypothesis Annotation group 
(id: dKymJJpZ) on hypothes.is and extract summary statistics. For the scripts to run, the user must be part of this 
group and have a hyothes.is developer token. The token must be copied and pasted into the directory,

```
hypothesis\files\input\secret_token.txt
```
which is saved locally.